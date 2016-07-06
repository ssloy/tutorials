//#define F_CPU 14745600L
#define F_CPU 16000000L

#include <avr/io.h>
#include <util/atomic.h>
#include <avr/interrupt.h>
//#include <util/delay.h>
#include <stdio.h>
#include <avr/pgmspace.h> // PSTR



FILE mystream;

//volatile int analog_val;
volatile unsigned char analog_val = 0; // care should be taken about non-atomic access


void uart_write(char x) {
    // wait for empty receive buffer
    while ((UCSR0A & (1<<UDRE0))==0);
    // send
    UDR0 = x;
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
    // For devisors see table 19-12 in the atmega328p datasheet.
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


// global variable to count the number of overflows
volatile unsigned char tot_overflow;

const unsigned char sinewave_data[] PROGMEM = {
0x0,0x3,0x6,0x9,0xc,0xf,0x12,0x15,0x18,0x1c,0x1f,0x22,0x25,0x28,0x2b,0x2e,
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
0x31,0x2e,0x2b,0x28,0x25,0x22,0x1f,0x1c,0x18,0x15,0x12,0xf,0xc,0x9,0x6,0x3};



#define logsize 1024
unsigned char readings[logsize];


int main(void) {
    for (int i=0; i<logsize; i++) readings[i] = 255;

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    DDRB = 0b00000000; // sets all the port B as input


    { // setup PWM: https://sites.google.com/site/qeewiki/books/avr-guide/pwm-on-the-atmega328
        // TODO inverted PWM in order to avoid spikes for zero duty cycles
        DDRD |= (1 << DDD6) | (1 << DDD5); // PD6 and PD5 are now output

        OCR0A = OCR0B = 0; // set PWM for 0% duty cycle

        TCCR0A |= (1 << COM0A1) | (1 << COM0A0) | (1 << COM0B1) | (1 << COM0B0); // set inverting mode
        //        TCCR0A |= (1 << COM0A1) | (1 << COM0B1); // set none-inverting mode
        TCCR0A |= (1 << WGM01) | (1 << WGM00);   // set fast PWM Mode
        TCCR0B |= (1 << CS01);                   // set prescaler to 8 and starts PWM
    }
    { // setup ADC
        // Bits 7:6    01 - Set REFS1..0 in ADMUX to AREF, to change reference voltage to the proper source
        // Bits 7:6    00 - Set REFS1..0 in ADMUX to AREF, internal Vref turned off
        // Bit  5       0 - clear ADLAR in ADMUX (0x7C) to right-adjust the result ADCL will contain lower 8 bits, ADCH upper 2 (in last two bits)
        // Bit  4       0 - unuseed
        // Bits 3:0  0000 - ADC0 (Table 24-4)
        ADMUX = 0b01000000;

        // Bit  7      1 - ADC Enable
        // Bit  6      0 - ADC Start conversion (well, no start here)
        // Bit  5      1 - ADC Auto trigger enable
        // Bit  4      0 - ADC Interrupt flag (this bit is set when an ADC conversion completes and the Data Registers are updated)
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
    { //timer 1 setup
        // set up timer with prescaler = 1024
        TCCR1B = (1 << CS10) | (1<<CS12);

        // initialize counter
        TCNT1 = 0;

        // enable overflow interrupt
        TIMSK1 |= (1 << TOIE1);

        // enable global interrupts
        sei();

        // initialize overflow counter variable
        tot_overflow = 0;

    }

    //      DDRC |= (1<<PC4);
    while (1) {


        /*
        //            PORTC |= (1<<PC4);

        //        fprintf_P(&uart_stream, PSTR("voltage: %.2f\n"), analog_val*5./255.);
        fprintf_P(&uart_stream, PSTR("voltage: %.2f %d\n"), analog_val*5./1023., analog_val);
        //        PORTC &= ~(1<<PC4);
         */


        unsigned long micros;
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            micros = (TCNT1 + (((unsigned long)tot_overflow)<<16))<<6;
        }

        int idx = ((int)(micros*(30000/1000000.)))&0x1FF;

        int voltage =  pgm_read_byte(&sinewave_data[idx&0xFF])*(idx>0xFF?1:-1);
        if (voltage>0) {
            OCR0A = 255-voltage;
            OCR0B = 255-0;
        } else {
            OCR0A = 255-0;
            OCR0B = 255+voltage;
        }


    }
}


// Interrupt service routine for the ADC completion
ISR(ADC_vect) {
//    analog_val = ADCL | (ADCH << 8); // Must read low first
    analog_val = (ADCL>>2) | (ADCH<<6);
}


// TIMER1 overflow interrupt service routine
// called whenever TCNT1 overflows
ISR(TIMER1_OVF_vect) {
    // keep a track of number of overflows
    tot_overflow++;
}

