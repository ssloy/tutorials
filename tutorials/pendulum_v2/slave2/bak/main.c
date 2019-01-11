#define F_CPU 16000000L

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/atomic.h>
#include <util/delay.h>

#include <stdlib.h>
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volatile uint8_t timer1_overflows; // global variable to count the number of overflows of timer1

volatile uint16_t adc0; // adc readings

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

int main(void) {
    adc_init();
    pwm_init();

    OCR0A = OCR0B = 255; // set PWM for 0% duty cycle

    DDRD |= _BV(PD7); // пин PD7 будет мигать, это нужно для контроля длительности одной итерации цикла (осциллограф в помощь)
    DDRD |= _BV(PD4); // пин PD4 будет гореть, означая, что контроллер готов к работе
    uint8_t pd7_pin_toggle = 1; 
    PORTD |= _BV(PD4);

    timer_init();

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    reset_current_time();
    


    int32_t time_of_start, useconds;
    useconds = time_of_start = get_current_time();

    int16_t pwm = 100;
    while (useconds<10000000L) {
        pd7_pin_toggle = 1 - pd7_pin_toggle;
        if (pd7_pin_toggle) {
            PORTD |=  _BV(PD7);
        } else {
            PORTD &= ~_BV(PD7);
        }

        useconds = get_current_time();

        int32_t sine_idx = (useconds>>5) & 0x1FF;  // sine 61.04 Hz : 32uS per increment, one sine period in 32uS*512 samples, thus the frequency is 1/(32*512/10^6) = 61.0351

        pwm = sine(sine_idx);

        if (pwm>=0) {
            OCR0A = 255;
            OCR0B = 255-pwm;
        } else {
            OCR0A = 255+pwm;
            OCR0B = 255;
        }


        uint16_t y = 0;
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            y = adc0;
        }

        //        fprintf_P(&uart_stream, PSTR("%ld,%d,%d\r\n"), useconds, y, pwm);
        //        _delay_us(100);
    }
    PORTD &= ~_BV(PD4);

    return 0;
}

