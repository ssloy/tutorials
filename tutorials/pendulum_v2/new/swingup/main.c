#define F_CPU 16000000L

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/atomic.h>
#include <util/delay.h>

#include <stdlib.h>
#include <stdio.h>
#include <avr/pgmspace.h> // PSTR

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PENDULUM_ANGLE_MOD2PI     0
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
#define DECODER_OE_PIN    B,1
#define DECODER_XY_PIN    C,2
#define DECODER_SEL1_PIN  B,0
#define DECODER_SEL2_PIN  C,1
#define DECODER_RST_PIN   D,4

#define SK6812_DATA_PIN   B,3

#define AMPLIFIER_ENABLE_PIN B,2

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volatile uint16_t adc0; // adc readings
volatile int32_t adc0_zero = 51200;
volatile uint8_t timer1_overflows; // global variable to count the number of overflows of timer1

volatile int32_t pendulum_offset = 4000; // index pulse position in ticks
volatile int32_t pendulum_ticks  = 0;    // current angle measurement, in ticks

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int16_t current(int16_t reference) { // mA
    if (reference> 5000) reference =  5000;
    if (reference<-5000) reference = -5000;

    reference = ((int32_t)(reference)*15620L)/10000L; // на диффвыходе ардуины и входе в усилитель получается делитель с резисторами 4.7k и 10k

    int16_t pwm = ((int32_t)(reference)*255L)/5000L;

    if (pwm>=0) {
        OCR0A = 255;
        OCR0B = 255-pwm;
    } else {
        OCR0A = 255+pwm;
        OCR0B = 255;
    }

    int16_t reality = 0;
    ATOMIC_BLOCK(ATOMIC_FORCEON) {
        reality = adc0;
    }

    // АЦП adc0_val бегает от 0 до 1023, при этом ноль ампер соответствует 1023/2 АЦП
    // поэтому сначала преобразуем adc0_val в значение от -512 до 512, с нулём ампер в нуле: adc0_val - 512
    // чтобы получить амперы из выражения (adc0_val - 512), надо его умножить на (1/512)*(2.5V/.185V/A)
    reality = ((int32_t)(reality*100L - adc0_zero)*1000L)/4733L; // измеряется в миллиамперах

    return reality;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void adc_init() {
    // Bits 7:6    01 - Set REFS1..0 in ADMUX to AREF, to change reference voltage to the proper source
    // Bits 7:6    00 - Set REFS1..0 in ADMUX to AREF, internal Vref turned off
    // Bit  5       0 - clear ADLAR in ADMUX (0x7C) to right-adjust the result ADCL will contain lower 8 bits, ADCH upper 2 (in last two bits)
    // Bit  4       0 - unused
    // Bits 3:0  0000 - ADC0 (Table 24-4)
    ADMUX = 0b01000000;

    // Bit  7      1 - ADC Enable
    // Bit  6      0 - ADC Start conversion (well, no start here)
    // Bit  5      1 - ADC Auto trigger enable
    // Bit  4      0 - ADC Interrupt flag (this bit is set when an ADC conversion completes and the Data Registers are updated)
    // Bit  3      1 - ADC Interrupt enable (without this, the internal interrupt will not trigger)
    // Bits 2:0  111 - Set the prescaler to 128 (16000KHz/128 = 125KHz)
    ADCSRA = 0b10101111;

    // Bits 7, 5:3 - Reserved, these bits must be written to zero when ADCSRB is written.
    // Bit  6      - ACME: Analog Comparator Multiplexer Enable
    // Bits 2:0    - ADC Auto Trigger Source (free running here, this means that as soon as an ADC has finished, the next will be immediately started)
    ADCSRB = 0b00000000;

    sei();

    // Set ADSC in ADCSRA (0x7A) to start the ADC conversion
    ADCSRA |= 0b01000000;
}

ISR(ADC_vect) { // Interrupt service routine for the ADC completion
    adc0 = ADCL | (ADCH<<8);  // Must read low first
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pwm_init() {
    DDRD |= (1 << DDD6) | (1 << DDD5); // PD6 and PD5 are now output

    OCR0A = OCR0B = 255; // set PWM for 0% duty cycle

    TCCR0A  = (1 << COM0A1) | (1 << COM0A0) | (1 << COM0B1) | (1 << COM0B0); // set inverting mode  to avoid spikes for zero duty cycles
    TCCR0A |= (1 << WGM01) | (1 << WGM00);   // set fast PWM Mode
    TCCR0B |= (1 << CS00);                   // no prescaling
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int32_t get_cart_position() {
    int32_t tmp = read_decoder(0)*18L; // microns, #ticks*(36 teeth * 0.002m pitch)/(1000 ticks/rev * 4x) gives meters, multiplying it by 10^6 we get microns
    return tmp;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define ASM_STRIP_PIN2(port,pin)  "I" (_SFR_IO_ADDR(PORT ## port)), "I" (pin)
#define ASM_STRIP_PIN(x) ASM_STRIP_PIN2(x)  

/** led_strip_write sends a series of colors to the LED strip, updating the LEDs.
 The colors parameter should point to an array of rgb_color structs that hold the colors to send.
 The count parameter is the number of colors to send.
 This function takes about 1.1 ms to update 30 LEDs.
 Interrupts must be disabled during that time, so any interrupt-based library
 can be negatively affected by this function.
 Timing details at 20 MHz (the numbers slightly different at 16 MHz and 8MHz):
  0 pulse  = 400 ns
  1 pulse  = 850 ns
  "period" = 1300 ns
 */
void __attribute__((noinline)) led_strip_write(volatile uint8_t *colors, uint16_t count) {
    cli();   // Disable interrupts temporarily because we don't want our pulse timing to be messed up.
    while (count--) {
        // Send a color to the LED strip.
        // The assembly below also increments the 'colors' pointer,
        // it will be pointing to the next color at the end of this loop.
        asm volatile(
                "ld __tmp_reg__, %a0+\n"
                "rcall led_strip_send_byte%=\n"  // Send green component.
                "ld __tmp_reg__, %a0+\n"
                "rcall led_strip_send_byte%=\n"  // Send red component.
                "ld __tmp_reg__, %a0+\n"
                "rcall led_strip_send_byte%=\n"  // Send blue component.
                "ld __tmp_reg__, %a0+\n"
                "rcall led_strip_send_byte%=\n"  // Send white component.
                "rjmp led_strip_asm_end%=\n"     // Jump past the assembly subroutines.

                // led_strip_send_byte subroutine:  Sends a byte to the LED strip.
                "led_strip_send_byte%=:\n"
                "rcall led_strip_send_bit%=\n"  // Send most-significant bit (bit 7).
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"
                "rcall led_strip_send_bit%=\n"  // Send least-significant bit (bit 0).
                "ret\n"

                // led_strip_send_bit subroutine:  Sends single bit to the LED strip by driving the data line
                // high for some time.  The amount of time the line is high depends on whether the bit is 0 or 1,
                // but this function always takes the same time (2 us).
                "led_strip_send_bit%=:\n"
                "sbi %2, %3\n"                           // Drive the line high.
                "rol __tmp_reg__\n"                      // Rotate left through carry.
                "nop\n" "nop\n"
                "brcs .+2\n" "cbi %2, %3\n"              // If the bit to send is 0, drive the line low now.
                "nop\n" "nop\n" "nop\n" "nop\n" "nop\n"
                "brcc .+2\n" "cbi %2, %3\n"              // If the bit to send is 1, drive the line low now.
                "ret\n"
                "led_strip_asm_end%=: "
                : "=b" (colors)
                : "0" (colors),         // %a0 points to the next color to display
                ASM_STRIP_PIN(SK6812_DATA_PIN) // %2 is the port register (e.g. PORTB), %3 is the pin number (0-8)
                       );
    }
    sei();          // Re-enable interrupts now that we are done.
    _delay_us(80);  // Send the reset signal.
}


#define LED_COUNT 60
uint8_t red[]    = {0,128,0,0};
uint8_t green[]  = {128,0,0,0};
uint8_t gray[]   = {0,0,0,128};
uint8_t blue[]   = {0,0,128,0};
uint8_t yellow[] = {64,255,0,0};
uint8_t max[] = {255, 255, 255, 255};

volatile uint8_t colors[LED_COUNT*4];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
void swingup(FILE *uart_stream_ptr) {
    int32_t useconds, useconds_prev, time_of_start;
    time_of_start = useconds = useconds_prev = get_current_time();

    int32_t cart_position, cart_position_prev, pendulum_angle, pendulum_angle_prev;
    cart_position  = cart_position_prev  = get_cart_position();
    pendulum_angle = pendulum_angle_prev = get_pendulum_angle();


    float hatQ = pendulum_angle*3.14159/180.000/1000.; // radians
    float hatX = cart_position/1000000.;               // meters
    float hatLPX = 0;
    float hatLPQ = 0;

    int16_t current_ref, current_measured = 0;
    int32_t duration_ms = 30000L;

    uint32_t iteration = 0;
    while (1) {

        useconds_prev      = useconds;
        cart_position_prev = cart_position;
        useconds           = get_current_time();
        cart_position  = get_cart_position();
        pendulum_angle = get_pendulum_angle();

        if (useconds - time_of_start > duration_ms*1000L) break;


        int32_t v = ((cart_position  - cart_position_prev)*1000L)/(useconds-useconds_prev); //   mm/sec

        float Ki = 9.268;
        float fricpos =  3.986*.8;
        float fricneg = -2.875*.8;

        float u = current_measured*Ki/1000.;
        if (labs(v)>1) {
            if (v>0)  
                u -= fricpos;
            else
                u -= fricneg;
        } else {
            if (u>0) {
                if (u>fricpos) u -= fricpos;
                else u = 0;
            } else {
                if (u<fricneg) u -= fricneg;
                else u = 0;
            }
        }

        float dt = (useconds-useconds_prev)/1000000.; // seconds
        float mesQ = pendulum_angle*3.14159/180000.;  // radians
        float mesX = cart_position/1000000.;          // meters
        float cosQ = cos(mesQ);
        float sinQ = sin(mesQ);

        while (pendulum_angle> 180000) pendulum_angle -= 360000; // mod 2pi
        while (pendulum_angle<-180000) pendulum_angle += 360000;

        if (labs(pendulum_angle)<3000L) break;

        float A = 1./sqrt(731.775 - 152.361*cosQ*cosQ);
        float hatDQ = (123.018*hatLPQ) * A;
        float hatDX = 0.933*hatLPX - (11.512*hatLPQ*cosQ) * A;
        float diffX = hatDX + 50*(mesX - hatX);
        float diffQ = 80*(mesQ - hatQ) + hatDQ;
        float diffLPX = 373.066*mesX - 373.066*hatX + .933*u + (4604.916*cosQ*(hatQ - mesQ)) * A;
        float diffLPQ = (18452.7315*(mesQ - hatQ) + 129.832*sinQ - 11.512*u*cosQ) * A;

        hatX   = hatX   + dt*diffX;
        hatQ   = hatQ   + dt*diffQ;
        hatLPX = hatLPX + dt*diffLPX;
        hatLPQ = hatLPQ + dt*diffLPQ;

        

//        float E = .5*.3*(pow(hatDX+ cosQ*.38*hatDQ,2)+pow(hatDQ*.38*sinQ,2)) + .5*.5*hatDX*hatDX    + .3*9.8*.5*(cosQ-1.)  - .01;
        float E = .5*.048*hatDQ*hatDQ + .3*9.8*.5*(cosQ-1.)  - .15;

        current_ref = (int16_t)( (cosQ*hatDQ<0 ? -1 : 1 )*E*700. );

        if (current_ref >  700L) current_ref =  700L;
        if (current_ref < -700L) current_ref = -700L;

        if (current_ref > 0 && cart_position >  200000L) current_ref = 0;
        if (current_ref < 0 && cart_position < -200000L) current_ref = 0;

        current_measured = current(current_ref);
        fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld,%ld\r\n"), useconds-time_of_start, current_ref, current_measured, cart_position, pendulum_angle, (int32_t)(E*1000L));
        if (iteration%20==0){
        for (uint8_t j=0; j<LED_COUNT; j++) {
            float t = ((float)j/(float)(LED_COUNT)) * .3*9.8*.5*(-2.);
            
            for (uint8_t k=0; k<4; k++) {
                if (t>E)  {
                    colors[j*4+k] = 0;
                } else {
                    colors[j*4+k] = yellow[k]/32;
                }
            }
        }
        led_strip_write(colors, LED_COUNT);
        }
        iteration++;
    }
    current(0);
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  UNSTABLE EQUILIBRIUM  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stabilize(FILE *uart_stream_ptr) {
    int32_t useconds, useconds_prev, time_of_start;
    time_of_start = useconds = useconds_prev = get_current_time();

    int32_t cart_position, cart_position_prev, pendulum_angle, pendulum_angle_prev;
    cart_position  = cart_position_prev  = get_cart_position();
    pendulum_angle = pendulum_angle_prev = get_pendulum_angle();

    float hatX = 0;
    float hatQ = 0;
    float hatLPX = 0;
    float hatLPQ = 0;

    int8_t state = 0;
    int16_t current_ref, current_measured = 0;
    int32_t duration_ms = 30000L;
    int32_t target = 0, dtarget = 400;
    uint8_t pb4_pin_toggle = 1; 
    OUTPUT2(B,4);
    uint32_t iteration = 0;
    while (1) {
        pb4_pin_toggle = 1 - pb4_pin_toggle;
        if (pb4_pin_toggle) {
            SET2(B,4);
        } else {
            CLEAR2(B,4);
        }

        useconds_prev      = useconds;
        cart_position_prev = cart_position;
        useconds           = get_current_time();
        cart_position  = get_cart_position();
        pendulum_angle = get_pendulum_angle();
        while (pendulum_angle> 180000) pendulum_angle -= 360000; // mod 2pi
        while (pendulum_angle<-180000) pendulum_angle += 360000;



        if (2==state) {
            fprintf_P(uart_stream_ptr, PSTR("%ld %ld\r\n"), cart_position, pendulum_angle);
        }
        if (0==state) {
            hatQ = pendulum_angle*3.14159/180.000/1000.; // radians
            hatX = cart_position/1000000.;               // meters
            fprintf_P(uart_stream_ptr, PSTR("%ld %ld\r\n"), cart_position, pendulum_angle);

            if (labs(pendulum_angle)<5000L) {
                state = 1;
                time_of_start = useconds;
                if (dtarget) {
                    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg),target(um)\r\n"));
                } else {
                    fprintf_P(uart_stream_ptr, PSTR("time(us),reference current(mA),current(mA),cart position(um),pendulum angle(mdeg)\r\n"));
                }
            }
        }
        if (1==state) {
            if (useconds - time_of_start > duration_ms*1000L) break;
            if (labs(cart_position)>300000L || labs(pendulum_angle)>30000L) {
                state = 2;
                current(0);
                continue;
            }

            int32_t v = ((cart_position  - cart_position_prev)*1000L)/(useconds-useconds_prev); //   mm/sec

            float Ki = 9.268;
            float fricpos =  3.986*.8;
            float fricneg = -2.875*.8;

            float u = current_measured*Ki/1000.;
            if (labs(v)>1) {
                if (v>0)  
                    u -= fricpos;
                else
                    u -= fricneg;
            } else {
                if (u>0) {
                    if (u>fricpos) u -= fricpos;
                    else u = 0;
                } else {
                    if (u<fricneg) u -= fricneg;
                    else u = 0;
                }
            }

            float dt = (useconds-useconds_prev)/1000000.; // seconds
            float mesQ = pendulum_angle*3.14159/180000.;  // radians
            float mesX = cart_position/1000000.;          // meters
            float cosQ = cos(mesQ);
            float sinQ = sin(mesQ);

            float A = 1./sqrt(731.775 - 152.361*cosQ*cosQ);
            float hatDQ = (123.018*hatLPQ) * A;
            float hatDX = 0.933*hatLPX - (11.512*hatLPQ*cosQ) * A;
            float diffX = hatDX + 50*(mesX - hatX);
            float diffQ = 80*(mesQ - hatQ) + hatDQ;
            float diffLPX = 373.066*mesX - 373.066*hatX + .933*u + (4604.916*cosQ*(hatQ - mesQ)) * A;
            float diffLPQ = (18452.7315*(mesQ - hatQ) + 129.832*sinQ - 11.512*u*cosQ) * A;

            hatX   = hatX   + dt*diffX;
            hatQ   = hatQ   + dt*diffQ;
            hatLPX = hatLPX + dt*diffLPX;
            hatLPQ = hatLPQ + dt*diffLPQ;

            int32_t K[] = {21L, 23L, 104L, 20L};

            if (target > 200000L || target<-200000L) dtarget = -dtarget;
            target += dtarget;
            int32_t f = K[0]*((cart_position-target)/1000L) + K[1]*(int32_t)(hatDX*1000 - dtarget/5) + (K[2]*314L*pendulum_angle)/18000L + K[3]*(int32_t)(hatDQ*1000); // mN

            int32_t antifriction = 0;
            if (labs(v)>1) {
                if (v>0) antifriction = (int32_t)(fricpos*1000);
                else antifriction = (int32_t)(fricneg*1000);
            } else {
                if (labs(f)>10) {
                    if (f>0) antifriction = (int32_t)(fricpos*1000);
                    else antifriction = (int32_t)(fricneg*1000);
                } else f = 0;
            }
            current_ref = (int16_t)((f + antifriction)/Ki);

            if (current_ref >  5000L) current_ref =  5000L;
            if (current_ref < -5000L) current_ref = -5000L;
            current_measured = current(current_ref);
            if (dtarget) {
                fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld,%ld\r\n"), useconds-time_of_start, current_ref, current_measured, cart_position, pendulum_angle, target);
            } else {
                fprintf_P(uart_stream_ptr, PSTR("%ld,%d,%d,%ld,%ld\r\n"), useconds-time_of_start, current_ref, current_measured, cart_position, pendulum_angle);
            }


            if (iteration%20==0) {
                for (uint8_t j=0; j<LED_COUNT; j++) {
                    int32_t t  = ((-(int32_t)j +(LED_COUNT>>1))*1000L)/LED_COUNT;
                    for (uint8_t k=0; k<4; k++) {
                        colors[j*4+k] = 0;
                        if (labs(t)>300L)  {
                            colors[j*4+k] |= red[k]/16;
                        } 
                    }
                }

                int32_t idx1 = LED_COUNT/2L-(cart_position*(LED_COUNT/10L))/100000L;
                int32_t idx2 = LED_COUNT/2L-(       target*(LED_COUNT/10L))/100000L;
                idx1 = (idx1<0?0:(idx1>LED_COUNT-1?LED_COUNT-1:idx1));
                idx2 = (idx2<0?0:(idx2>LED_COUNT-1?LED_COUNT-1:idx2));
                for (uint8_t k=0; k<4; k++) {
                    colors[idx1*4+k] |=  blue[k];
                    colors[idx2*4+k] |= green[k];
                }

                led_strip_write(colors, LED_COUNT);
            }
            iteration++;

        }
    }
    current(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  main loop  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(void) {
    OUTPUT(SK6812_DATA_PIN);
    CLEAR(SK6812_DATA_PIN);

    OUTPUT(AMPLIFIER_ENABLE_PIN);
    CLEAR(AMPLIFIER_ENABLE_PIN);

    adc_init();
    pwm_init();
    timer_init();
    pwm_init();
    decoder_init();
    shift_register_init();

    SET(AMPLIFIER_ENABLE_PIN);
    _delay_ms(1000);
    {
        int32_t tmp = 0;
        for (int16_t i=0; i<5000; i++) {
            _delay_ms(1);
            ATOMIC_BLOCK(ATOMIC_FORCEON) {
                tmp += adc0;
            }
        }
        adc0_zero = tmp/50L;
    }

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    decoder_reset();
    reset_current_time();

    stabilize(&uart_stream);
/*
    for (uint8_t i=0; i<5; i++) {
        swingup(&uart_stream);
        stabilize(&uart_stream);
        _delay_ms(1000);
    }
*/

    CLEAR(AMPLIFIER_ENABLE_PIN);
    return 0;
}

