#define F_CPU 16000000L

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/atomic.h>
#include <util/delay.h>

#include <stdio.h>
#include <avr/pgmspace.h> // PSTR

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

#define SHIFT_PL_PIN                D,2
#define SHIFT_CP_PIN                D,3
#define SHIFT_Q7_PIN                C,3

#define DECODER_CHIX_PIN            D,7 
#define DECODER_OE_PIN              D,5
#define DECODER_XY_PIN              C,2
#define DECODER_SEL1_PIN            D,6
#define DECODER_SEL2_PIN            C,1
#define DECODER_RST_PIN             D,4

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volatile uint8_t tot_overflow;     // global variable to count the number of overflows of timer1

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void spi_master_init() {
    DDRB  |=  _BV(PB3) | _BV(PB5) | _BV(PB2);   // set MOSI, SCK and SS as output
    DDRB  &= ~_BV(PB4);                         // set MISO as input
    SPCR   =  _BV(SPE) | _BV(MSTR) | _BV(SPR0); // enable master SPI at clock rate Fck/4
}

uint8_t spi_tranceiver(uint8_t data) {
    PORTB &= ~_BV(PB2);
    SPDR = data;                 // send data
    while (!(SPSR & _BV(SPIF))); // wait for transmition complete
    PORTB |= _BV(PB2);
    return SPDR;
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

ISR (PCINT2_vect) { // D9 has changed, thus index pulse fired, updating its position
}


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
    _delay_us(10);
    SET(DECODER_RST_PIN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void uart_write(char x) {
    while ((UCSR0A & (1<<UDRE0))==0); // wait for empty receive buffer
    UDR0 = x; // send
}

uint8_t uart_char_is_waiting() {
    // returns 1 if a character is waiting
    // returns 0 if not
    return (UCSR0A & (1<<RXC0));
}

char uart_read() {
    while (!uart_char_is_waiting());
    char x = UDR0;
    return x;
}

int uart_putchar(char c, FILE *stream) {
    uart_write(c);
    return 0;
}

int uart_getchar(FILE *stream) {
    return uart_read();
}

void uart_init() {
    // For divisors see table 19-12 in the atmega328p datasheet.
    // U2X0, 16 -> 115.2k baud @ 16MHz. 
    // U2X0, 207 -> 9600 baud @ 16Mhz.
    UBRR0H = 0;
    UBRR0L = 16;
    UCSR0A = 1<<U2X0;
    // Enable  the transmitter. Reciever is disabled.
    UCSR0B = 1<<TXEN0;
    UCSR0C = (1<<UDORD0) | (1<<UCPHA0);  //(3 << UCSZ00);  

    fdevopen(&uart_putchar, &uart_getchar);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void timer_init() {
    TCCR1B = (1 << CS10) | (1<<CS12); // set up timer with prescaler = 1024
    TCNT1 = 0;                        // initialize counter
    TIMSK1 |= (1 << TOIE1);           // enable overflow interrupt
    sei();                            // enable global interrupts
    tot_overflow = 0;                 // initialize overflow counter variable
}

ISR(TIMER1_OVF_vect) { // TIMER1 overflow interrupt service routine called whenever TCNT1 overflows
    tot_overflow++;    // keep a track of number of overflows
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const uint8_t sinewave_data[] PROGMEM = {
0,2,4,6,8,10,12,14,16,19,21,23,25,28,30,32,
34,36,38,40,42,44,47,49,51,53,55,57,59,61,63,65,
67,70,72,74,76,77,79,81,84,85,87,89,91,93,95,97,
98,100,102,104,105,107,109,111,112,114,116,117,119,121,122,124,
125,127,128,130,131,133,135,136,137,139,140,141,142,144,145,147,
148,149,150,151,152,154,154,156,156,158,158,160,161,161,163,163,
164,165,165,167,168,168,169,170,170,170,171,172,172,173,173,174,
175,175,175,175,176,176,177,177,177,177,177,177,177,177,177,177,
178,177,177,177,177,177,177,177,177,177,177,176,176,175,175,175,
175,174,173,173,172,172,171,170,170,170,169,168,168,167,165,165,
164,163,163,161,161,160,158,158,156,156,154,154,152,151,150,149,
148,147,145,144,142,141,140,139,137,136,135,133,131,130,128,127,
125,124,122,121,119,117,116,114,112,111,109,107,105,104,102,100,
98,97,95,93,91,89,87,85,84,81,79,77,76,74,72,70,
67,65,63,61,59,57,55,53,51,49,47,44,42,40,38,36,
34,32,30,28,25,23,21,19,16,14,12,10,8,6,4,2};


int main(void) {
    timer_init();
    spi_master_init();
    decoder_init();
    decoder_reset();
    shift_register_init();

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    uint32_t micros=0;
    while(1) {
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            micros = (TCNT1 + ((uint32_t)tot_overflow<<16)) << 6;
        }
        uint16_t sine_idx = (micros>>5) & 0x1FF; // sine 61.04 Hz : 32uS per increment, one sine period in 32uS*512 samples, thus the frequency is 1/(32*512/10^6) = 61.0351
        uint8_t g = 0;
        if (sine_idx & 0x100) {
            g =      pgm_read_byte(&sinewave_data[sine_idx & 0xFF])>>1;
        } else {
            g = 127+(pgm_read_byte(&sinewave_data[sine_idx & 0xFF])>>1);
        }

        uint8_t ack = spi_tranceiver(g);
        fprintf_P(&uart_stream, PSTR("%ld %d %d\r\n"), micros, g, ack);

//      int32_t ticks0 = read_decoder(0);
//      int32_t ticks1 = read_decoder(1);

//        fprintf_P(&uart_stream, PSTR("X: %ld \t Y: %ld \t ack: %d %d\r\n"), ticks0, ticks1, ack, cnt);
//        _delay_ms(200);
    }
    return 0;
}

