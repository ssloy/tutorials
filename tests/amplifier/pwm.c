#define F_CPU 14745600

#include <avr/io.h>
//#include <util/atomic.h>
#include <avr/interrupt.h>
#include <util/delay.h>
#include <stdio.h>
#include <avr/pgmspace.h> // PSTR


FILE mystream;

volatile int analog_val;
//volatile unsigned char analog_val; // care should be taken about non-atomic access


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
    // set baud rate
    UBRR0H = 0;
    UBRR0L = 7;	// for 115200bps with 14.7456MHz clock
    // enable uart RX and TX
    UCSR0B = (1<<RXEN0)|(1<<TXEN0);
    // set 8N1 frame format
    UCSR0C = (1<<UCSZ01)|(1<<UCSZ00);
    // set up STDIO handlers so you can use printf, etc
    fdevopen(&uart_putchar, &uart_getchar);
}



int main(void) {
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

    //      DDRC |= (1<<PC4);
    while (1) {
        OCR0B = 255;
        OCR0A = 0;
        //            PORTC |= (1<<PC4);

//        fprintf_P(&uart_stream, PSTR("voltage: %.2f\n"), analog_val*5./255.);
        fprintf_P(&uart_stream, PSTR("voltage: %.2f %d\n"), analog_val*5./1023., analog_val);
        //        PORTC &= ~(1<<PC4);
    }
}


// Interrupt service routine for the ADC completion
ISR(ADC_vect) {
    analog_val = ADCL | (ADCH << 8); // Must read low first
//    analog_val = (ADCL>>2) | (ADCH << 6);
}

