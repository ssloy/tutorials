#define F_CPU 16000000L

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/atomic.h>
#include <util/delay.h>

#include <stdlib.h>
#include <stdio.h>
#include <avr/pgmspace.h> // PSTR

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PENDULUM_ANGLE_MOD2PI 0
#define PENDULUM_USE_INDEX_SIGNAL 0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define  INPUT2(port,pin)   DDR ## port &= ~_BV(pin) 
#define OUTPUT2(port,pin)   DDR ## port |=  _BV(pin) 
#define  CLEAR2(port,pin)  PORT ## port &= ~_BV(pin) 
#define    SET2(port,pin)  PORT ## port |=  _BV(pin) 
#define   READ2(port,pin) ((PIN ## port & _BV(pin))?1:0)

#define  INPUT(x)  INPUT2(x) 
#define OUTPUT(x) OUTPUT2(x)
#define  CLEAR(x)  CLEAR2(x)
#define    SET(x)    SET2(x)
#define   READ(x)   READ2(x)
#define  WRITE(x,b) ((b)?(SET2(x)):(CLEAR2(x)))

#define SHIFT_PL_PIN      D,2
#define SHIFT_CP_PIN      D,3
#define SHIFT_Q7_PIN      C,3

#define DECODER_CHIX_PIN  D,7 
#define DECODER_OE_PIN    D,5
#define DECODER_XY_PIN    C,2
#define DECODER_SEL1_PIN  D,6
#define DECODER_SEL2_PIN  C,1
#define DECODER_RST_PIN   D,4

#define AMPLIFIER_RST_PIN B,1

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volatile uint8_t timer1_overflows; // global variable to count the number of overflows of timer1

volatile int32_t pendulum_offset = 4000; // index pulse position in ticks
volatile int32_t pendulum_ticks  = 0;    // current angle measurement, in ticks

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void spi_master_init() {
    DDRB  |=  _BV(PB3) | _BV(PB5) | _BV(PB2);   // set MOSI, SCK and SS as output
    DDRB  &= ~_BV(PB4);                         // set MISO as input
    SPCR   =  _BV(SPE) | _BV(MSTR) | _BV(SPR0); // enable master SPI at clock rate Fck/4
}

//uint8_t spi_tranceiver(uint8_t data) {
//    PORTB &= ~_BV(PB2);          // select slave
//    SPDR = data;                 // send data
//    while (!(SPSR & _BV(SPIF))); // wait for transmition complete
//    PORTB |= _BV(PB2);           // release slave
//    return SPDR;
//}

// this functions sends (and receives) a trash byte to initiate communication
// then len bytes are exchanged
volatile uint8_t trash;
void spi_transfer_sync(uint8_t *dataout, uint8_t *datain, uint8_t len) {
    PORTB &= ~_BV(PB2);
    SPDR = 0;
    while (!(SPSR & _BV(SPIF)));
    trash = SPDR;
    
    _delay_us(1);
    for (uint8_t i = 0; i<len; i++) {
        SPDR = dataout[i];
        while (!(SPSR & _BV(SPIF)));
        datain[i] = SPDR;
        _delay_us(1);
    }
    PORTB |= _BV(PB2);
}
        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void decoder_init() {
    OUTPUT(DECODER_OE_PIN);
    OUTPUT(DECODER_XY_PIN);
    OUTPUT(DECODER_SEL1_PIN);
    OUTPUT(DECODER_SEL2_PIN);
    OUTPUT(DECODER_RST_PIN);
    SET(DECODER_RST_PIN); // reset is active low 

    INPUT(DECODER_CHIX_PIN);
    PCICR  |= _BV(PCIE2);  // interrupt will be fired on any change on pin DECODER_CHIX_PIN
    PCMSK2 |= _BV(7);      // DECODER_CHIX_PIN = D7
    sei();
}

void shift_register_init() {
    OUTPUT(SHIFT_PL_PIN);
    OUTPUT(SHIFT_CP_PIN);
    INPUT(SHIFT_Q7_PIN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const uint8_t avago_byte_select[4][2] = {{1,0},{0,0},{1,1},{0,1}};
inline void decoder_latch_byte(uint8_t i) {
    WRITE(DECODER_SEL1_PIN, avago_byte_select[i][0]);
    WRITE(DECODER_SEL2_PIN, avago_byte_select[i][1]);
}

uint8_t read_shift_register() {
    uint8_t result = 0;
    CLEAR(SHIFT_PL_PIN);
    _delay_us(10);
    SET(SHIFT_PL_PIN);
    for (uint8_t i=8; i--;) {
        result |= (READ(SHIFT_Q7_PIN)<<(i));
        SET(SHIFT_CP_PIN);
        _delay_us(10);
        CLEAR(SHIFT_CP_PIN);
        _delay_us(10);
    }
    return result;
}

int32_t read_decoder(uint8_t encoder) {
    int32_t result = 0;
    WRITE(DECODER_XY_PIN, encoder);
    CLEAR(DECODER_OE_PIN);
    _delay_us(10); // The OE/, SEL1, and SEL2 inputs are sampled by the internal inhibit logic on the falling edge of the clock
    for (uint8_t i=4; i--;) {
        decoder_latch_byte(i);
        result += read_shift_register();
        if (i) result <<= 8;
    }
    SET(DECODER_OE_PIN);
    return result;
}

void decoder_reset() {
    CLEAR(DECODER_RST_PIN);
    _delay_ms(1);
    SET(DECODER_RST_PIN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void uart_write(char x) {
    while ((UCSR0A & (1<<UDRE0))==0); // wait for empty receive buffer
    UDR0 = x; // send
}

uint8_t uart_char_is_waiting() { // returns 1 if a character is waiting, 0 otherwise
    return (UCSR0A & (1<<RXC0));
}

char uart_read() {
    while (!uart_char_is_waiting());
    char x = UDR0;
    return x;
}

int uart_putchar(char c, FILE *stream __attribute__((unused))) {
    uart_write(c);
    return 0;
}

int uart_getchar(FILE *stream __attribute__((unused))) {
    return uart_read();
}

void uart_init() {
    UBRR0H = 0;        // For divisors see table 19-12 in the atmega328p datasheet.
    UBRR0L = 16;       // U2X0, 16 -> 115.2k baud @ 16MHz. 
    UCSR0A = 1<<U2X0;  // U2X0, 207 -> 9600 baud @ 16Mhz.
    UCSR0B = 1<<TXEN0; // Enable  the transmitter. Reciever is disabled.
    UCSR0C = (1<<UDORD0) | (1<<UCPHA0);
    fdevopen(&uart_putchar, &uart_getchar);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void timer_init() {
    TCCR1B = (1 << CS10) | (1<<CS12); // set up timer with prescaler = 1024
    TCNT1 = 0;                        // initialize counter
    TIMSK1 |= (1 << TOIE1);           // enable overflow interrupt
    sei();                            // enable global interrupts
    timer1_overflows = 0;             // initialize overflow counter variable
}

ISR(TIMER1_OVF_vect) { // TIMER1 overflow interrupt service routine called whenever TCNT1 overflows
    timer1_overflows++;    // keep a track of number of overflows
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void reset_current_time() {
    ATOMIC_BLOCK(ATOMIC_FORCEON) {
        TCNT1 = 0;
        timer1_overflows = 0;
    }
}

inline int32_t get_current_time() {
    int32_t useconds;
    ATOMIC_BLOCK(ATOMIC_FORCEON) {
        useconds = (TCNT1 + ((uint32_t)timer1_overflows<<16)) << 6;
    }
    return useconds;
}

inline int32_t get_cart_position() {
    return read_decoder(0)*18L; // microns, #ticks*(36 teeth * 0.002m pitch)/(1000 ticks/rev * 4x) gives meters, multiplying it by 10^6 we get microns
}

ISR (PCINT2_vect) { // D9 has changed, thus index pulse fired, updating its position
#if PENDULUM_USE_INDEX_SIGNAL
    pendulum_offset = pendulum_ticks;
#endif
}

inline int32_t get_pendulum_angle() {
    cli();
    pendulum_ticks = read_decoder(1);
    long tticks = pendulum_ticks - pendulum_offset; // if index pulse position off the balance, correct it here
    sei();

#if PENDULUM_ANGLE_MOD2PI
    while (tticks> 4000) tticks -= 8000; // mod 2pi
    while (tticks<-4000) tticks += 8000;
#endif

    return tticks*45L; // millidegrees, #ticks * 360 / (2000 ticks/rev * 4x) gives degrees
}

void speed_closed_loop(int32_t reference_speed, uint32_t duration, FILE *uart_stream_ptr) {
    uint32_t time_of_start, useconds, useconds_prev;
    useconds_prev = useconds = time_of_start = get_current_time();

    int32_t speed = 0; // microns per second
    int32_t microns, microns_prev;
    microns = microns_prev = get_cart_position();

    int16_t current = 0, ack = 0;
    int32_t e = 0, x = 0, dt = 0;
    while (useconds < time_of_start+duration) {
        useconds_prev = useconds;
        useconds = get_current_time();
        dt = useconds - useconds_prev;

        microns_prev = microns;
        microns = get_cart_position();
        speed = (microns-microns_prev)*(1000000L/dt);

        e = reference_speed - speed;        // current error
        x = x + (((94L*e)/1000L)*dt)/1000L; // error integral, multiplied by Ki = 94.3 ; x = x + Ki*e*dt/10^6
        
        if (x >  3000000L) x =  3000000L; // three amps saturation
        if (x < -3000000L) x = -3000000L;

        current = (7L*e + x)/1000L; // mA, Kp = 6.87

        if (current >  3000L) current =  3000L; // three amps saturation
        if (current < -3000L) current = -3000L;

        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);

        if (uart_stream_ptr)
            fprintf_P(uart_stream_ptr, PSTR("%ld,%ld,%d,%d\r\n"), useconds, speed, current, ack);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

void rewind(uint8_t right) {
    speed_closed_loop(300L*1000L*(right?1L:-1L), 3000000L, NULL);
    _delay_ms(500);
}

void center_cart() {
    int32_t center = ((900L-45L)/2L)*1000L;
    uint8_t to_right = center>get_cart_position();
    while (to_right==(center>get_cart_position())) {
        speed_closed_loop(200L*1000L*(to_right ? 1L : -1L), 100L*1000L, NULL);
    }
}


const uint8_t sinewave_halfperiod[] PROGMEM = {
0x0, 0x3, 0x6, 0x9, 0xc, 0xf, 0x12,0x15,0x18,0x1c,0x1f,0x22,0x25,0x28,0x2b,0x2e,
0x31,0x34,0x37,0x3a,0x3d,0x40,0x44,0x47,0x4a,0x4d,0x4f,0x52,0x55,0x58,0x5b,0x5e,
0x61,0x64,0x67,0x6a,0x6d,0x6f,0x72,0x75,0x78,0x7a,0x7d,0x80,0x83,0x85,0x88,0x8b,
0x8d,0x90,0x92,0x95,0x97,0x9a,0x9c,0x9f,0xa1,0xa4,0xa6,0xa8,0xab,0xad,0xaf,0xb2,
0xb4,0xb6,0xb8,0xba,0xbc,0xbf,0xc1,0xc3,0xc5,0xc7,0xc9,0xca,0xcc,0xce,0xd0,0xd2,
0xd4,0xd5,0xd7,0xd9,0xda,0xdc,0xdd,0xdf,0xe0,0xe2,0xe3,0xe5,0xe6,0xe7,0xe9,0xea,
0xeb,0xec,0xed,0xef,0xf0,0xf1,0xf2,0xf3,0xf4,0xf4,0xf5,0xf6,0xf7,0xf8,0xf8,0xf9,
0xfa,0xfa,0xfb,0xfb,0xfc,0xfc,0xfd,0xfd,0xfd,0xfe,0xfe,0xfe,0xfe,0xfe,0xfe,0xfe,
0xff,0xfe,0xfe,0xfe,0xfe,0xfe,0xfe,0xfe,0xfd,0xfd,0xfd,0xfc,0xfc,0xfb,0xfb,0xfa,
0xfa,0xf9,0xf8,0xf8,0xf7,0xf6,0xf5,0xf4,0xf4,0xf3,0xf2,0xf1,0xf0,0xef,0xed,0xec,
0xeb,0xea,0xe9,0xe7,0xe6,0xe5,0xe3,0xe2,0xe0,0xdf,0xdd,0xdc,0xda,0xd9,0xd7,0xd5,
0xd4,0xd2,0xd0,0xce,0xcc,0xca,0xc9,0xc7,0xc5,0xc3,0xc1,0xbf,0xbc,0xba,0xb8,0xb6,
0xb4,0xb2,0xaf,0xad,0xab,0xa8,0xa6,0xa4,0xa1,0x9f,0x9c,0x9a,0x97,0x95,0x92,0x90,
0x8d,0x8b,0x88,0x85,0x83,0x80,0x7d,0x7a,0x78,0x75,0x72,0x6f,0x6d,0x6a,0x67,0x64,
0x61,0x5e,0x5b,0x58,0x55,0x52,0x4f,0x4d,0x4a,0x47,0x44,0x40,0x3d,0x3a,0x37,0x34,
0x31,0x2e,0x2b,0x28,0x25,0x22,0x1f,0x1c,0x18,0x15,0x12,0xf, 0xc, 0x9, 0x6, 0x3};

inline int32_t sine(int32_t idx) { // one period equals 512 idx units
    if (idx & 0x100) {
        return -pgm_read_byte(&sinewave_halfperiod[idx & 0xFF]);
    } else {
        return  pgm_read_byte(&sinewave_halfperiod[idx & 0xFF]);
    }
}


inline int32_t cosine(int32_t idx) { // one period equals 512 idx units
    return sine(idx+128);
}

void constant_current_open_loop(int32_t ref, int32_t duration_ms, FILE *uart_stream_ptr) {
    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
    int32_t time_of_start, useconds;
    useconds = time_of_start = get_current_time();

    int32_t cart_position, cart_position_prev, pendulum_angle;
    cart_position = cart_position_prev = get_cart_position();

    int16_t current = ref, ack;
    while (useconds < time_of_start+duration_ms*1000L) {
        useconds = get_current_time();

        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);

        cart_position = get_cart_position();
        pendulum_angle = get_pendulum_angle();

        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds, current, ack, cart_position, pendulum_angle);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

// frequencies are given in mHz
void chirp_current_open_loop(int32_t start_freq, int32_t end_freq, int32_t duration_ms, FILE *uart_stream_ptr) {
    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
    int32_t time_of_start, useconds;
    useconds = time_of_start = get_current_time();

    int32_t cart_position, cart_position_prev, pendulum_angle;
    cart_position = cart_position_prev = get_cart_position();

    int16_t current, ack;
    while (useconds < time_of_start+duration_ms*1000L) {
        useconds = get_current_time();
        int32_t cur_ms = (useconds-time_of_start)/1000L;
        
        int32_t phi = ((start_freq + ((end_freq-start_freq)*cur_ms)/(2L*duration_ms))*cur_ms)/1000L; // период синусоиды равен тысяче единиц phi

        int32_t sine_idx = ((phi*512L)/1000L) & 0x1FF;
        current = (sine(sine_idx)*2500L)/255L;

        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
        cart_position = get_cart_position();
        pendulum_angle = get_pendulum_angle();

        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds, current, ack, cart_position, pendulum_angle);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

// frequencies are given in mHz
void sine_wave(int32_t freq, int32_t amplitude, int32_t duration_ms, FILE *uart_stream_ptr) {
    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
    int32_t time_of_start, useconds;
    useconds = time_of_start = get_current_time();

    int32_t cart_position, pendulum_angle;
    cart_position =  get_cart_position();

    int16_t current, ack;
    while (useconds < time_of_start+duration_ms*1000L) {
        useconds = get_current_time();
        int32_t cur_ms = (useconds-time_of_start)/1000L;


       
        int32_t sine_idx = ((((cur_ms*freq)/1000L)*512L)/1000L);
        current = (sine(sine_idx)*amplitude)/255L;
               
        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
        cart_position  = get_cart_position();
        pendulum_angle = get_pendulum_angle();

        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds, current, ack, cart_position, pendulum_angle);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

// frequencies are given in mHz
void square_wave(int32_t freq, int32_t amplitude, int32_t duration_ms, FILE *uart_stream_ptr) {
    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
    int32_t time_of_start, useconds;
    useconds = time_of_start = get_current_time();

    int32_t cart_position, pendulum_angle;
    cart_position =  get_cart_position();

    int16_t current, ack;
    while (useconds < time_of_start+duration_ms*1000L) {
        useconds = get_current_time();
        int32_t cur_ms = (useconds-time_of_start)/1000L;
        if (cur_ms%(1000000L/freq) > (1000000L/2L)/freq) {
            current = -amplitude;
        } else {
            current =  amplitude;
        }
        
        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
        cart_position  = get_cart_position();
        pendulum_angle = get_pendulum_angle();

        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds, current, ack, cart_position, pendulum_angle);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

void random_current(int32_t duration_ms, FILE *uart_stream_ptr) {
    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
    int32_t time_of_start, useconds, useconds_prev;
    useconds_prev = useconds = time_of_start = get_current_time();

    int32_t pendulum_angle;

    const int8_t memlen = 6;
    int32_t X[memlen], cnt=0;
    for (uint8_t i=memlen; i--; X[i] = 0);

    const uint8_t rndlen = 16;
    int16_t rnd[rndlen];
    for (uint8_t i=rndlen; i--; rnd[i] = 0);

    int16_t ridx = 0;
    int16_t current, ack;
    uint8_t right = 1;


    while (useconds < time_of_start+duration_ms*1000L) {
        useconds_prev = useconds;
        useconds = get_current_time();


        X [cnt%memlen] = get_cart_position();
        cnt++;
        pendulum_angle = get_pendulum_angle();
        int32_t v = ((X[(cnt+memlen-1)%memlen]-X[cnt%memlen])*1000L)/(memlen*(useconds-useconds_prev)); //   mm/sec

        rnd[ridx%rndlen] = (rand()&1 ? 1 : -1)*rand()%1800;
        ridx++;
        current = 0;
        for (uint8_t c=rndlen; c--; current += rnd[c]);
        if (right) {
            current = 300 + current/rndlen;
            if (v<250 && current < 400) current = 400;
            if (v>750) current = -400;
            if (X[(cnt+memlen-1)%memlen] > 600L*1000L) right = 0;
        } else {
            current = -300 + current/8;
            if (v>-250 && current>-400) current = -400;
            if (v<-750) current = 400;
            if (X[(cnt+memlen-1)%memlen] < 200L*1000L) right = 1;
        }
        
        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds, current, ack, X[(cnt+memlen-1)%memlen], pendulum_angle);
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}


//#define labs(x)  ( ( (x) < 0) ? -(x) : (x) )

void stabilize(FILE *uart_stream_ptr) {
    int32_t useconds, useconds_prev, time_of_start;
    time_of_start = useconds = useconds_prev = get_current_time();

    int32_t cart_position, cart_position_prev, pendulum_angle, pendulum_angle_prev;
    cart_position  = cart_position_prev  = get_cart_position();
    pendulum_angle = pendulum_angle_prev = get_pendulum_angle();

    const int8_t memlen = 5;
    int32_t X[memlen], TH[memlen], cnt=0;


    int8_t state = 0;
    int16_t current, ack;
    int32_t duration_ms = 30000L;
    while (1) {
        useconds_prev          = useconds;
        cart_position_prev     = cart_position;
//        pendulum_angle_prev = pendulum_angle;
        useconds       = get_current_time();
        cart_position  = get_cart_position();
        pendulum_angle = get_pendulum_angle();
        X [cnt%memlen] = cart_position;
        TH[cnt%memlen] = pendulum_angle;
        cnt++;

        if (0==state) {
            fprintf_P(uart_stream_ptr, PSTR("%ld %ld\r\n"), cart_position, pendulum_angle);
            if (labs(pendulum_angle)<5000L) {
                state = 1;
                time_of_start = useconds;
            }
        }
        if (1==state) {
            int32_t cur_ms = (useconds-time_of_start)/1000L;
            int32_t freq = 50L; 
            int32_t sine_idx = ((((cur_ms*freq)/1000L)*512L)/1000L);
//            cart_position += ((sine(sine_idx)*1000L)/255L)*1000L/4;


            if (useconds - time_of_start > duration_ms*1000L) break;
            if (labs(cart_position)>300000L || labs(pendulum_angle)>30000L) break;

            int32_t v = ((cart_position  - cart_position_prev)*1000L)/(useconds-useconds_prev); //   mm/sec
            int32_t w = ((pendulum_angle - TH[cnt%memlen])*1000000L)/(memlen*(useconds-useconds_prev)); // mdeg/sec


//            int32_t v = ((cart_position  -  cart_position_prev)*1000L   )/(useconds-useconds_prev); //   mm/sec
//            int32_t w = ((pendulum_angle - pendulum_angle_prev)*1000000L)/(useconds-useconds_prev); // mdeg/sec
            int32_t K[] = {30L, 23L, 101L, 19L};

            int32_t f = K[0]*(cart_position/1000L) + K[1]*v + (K[2]*314L*pendulum_angle)/18000L + K[3]*((314L*w)/18000L); // mN

            current = (int16_t)((f + (v==0L ? (f>0L?4690L:-3350L) : (v>0L?4690L:-3350L)))/13L);
            if (current >  5000L) current =  5000L;
            if (current < -5000L) current = -5000L;
            spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);

        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds-time_of_start, current, ack, cart_position, pendulum_angle);
//            int32_t w2 = ((pendulum_angle - TH[(cnt+memlen-2)%memlen])*1000000L)/(useconds-useconds_prev); // mdeg/sec
//            fprintf_P(uart_stream_ptr, PSTR("%ld\t%ld\t%ld\t%ld\r\n"), v, w, w2, pendulum_angle);
//            fprintf_P(uart_stream_ptr, PSTR("%ld %ld %ld %ld\t%ld\t%d %d\r\n"), cart_position, v, pendulum_angle, w, f, current, ack);
        }
    }
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
}

int main(void) {
    // reset the current amplifier
    OUTPUT(AMPLIFIER_RST_PIN);
    CLEAR(AMPLIFIER_RST_PIN);
    _delay_ms(1);
    SET(AMPLIFIER_RST_PIN);
    _delay_ms(3000); // wait for the current amplifier to boot and calibrate its zero

    timer_init();
    spi_master_init();
    decoder_init();
    shift_register_init();

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    decoder_reset();
    /*
    rewind(0);
    decoder_reset();
    center_cart();
    */
    reset_current_time();
//    speed_closed_loop(-350L*1000L,  60L*1000L*1000L, &uart_stream);

//    square_wave(2000L, 1000L, 5L*1000L, &uart_stream);

    chirp_current_open_loop(2000L, 9000L, 5L*1000L, &uart_stream);
    if (0) {
        stabilize(&uart_stream);
        while (1) {
            fprintf_P(&uart_stream, PSTR("%ld %ld\r\n"), get_pendulum_angle(), pendulum_offset);
            _delay_ms(100);
        }
    }

//    constant_current_open_loop(700L, 3L*1000L, &uart_stream);
//    square_wave(10000L, 1500L, 5L*1000L, &uart_stream);
//    sine_wave(1000L, 300L, 5L*1000L, &uart_stream);


//    square_wave(3000L, 1200L, 5L*1000L, &uart_stream);
//  srand(642);
// rewind(0);


//    square_wave(2000L, 1000L, 5L*1000L, &uart_stream);
//    chirp_current_open_loop(1000L, 5000L, 5L*1000L, &uart_stream);
//  random_current(5L*1000L, &uart_stream);



#if 0 
    fprintf_P(&uart_stream, PSTR("time(us),speed(um/sec),reference current(mA),current(mA)\r\n"));
    speed_closed_loop(  50L*1000L, 18L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( -50L*1000L, 18L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( 100L*1000L,  9L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop(-100L*1000L,  9L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( 150L*1000L,  6L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop(-150L*1000L,  6L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( 200L*1000L,  5L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop(-200L*1000L,  5L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( 250L*1000L,  4L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop(-250L*1000L,  4L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop( 300L*1000L,  3L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
    speed_closed_loop(-300L*1000L,  3L*1000L*1000L, &uart_stream);
    _delay_ms(1000);
#endif

    return 0;
}

