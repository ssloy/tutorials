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

volatile uint16_t adc0; // adc readings
volatile  int16_t spi_in_msg  = 0;   // задание силы тока
volatile  int16_t spi_out_msg = 0;   // реально протекающий ток
volatile uint8_t trash = 0;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void spi_slave_init() {
    DDRB &= ~(_BV(PB3) | _BV(PB5)  | _BV(PB2)); // set MOSI, SCK and SS as input
    DDRB |=   _BV(PB4);                         // set MISO as output
    SPCR  =   _BV(SPE) | _BV(SPR0) | _BV(SPIE); // enable slave SPI at clock rate Fck/4, turn on interrupts
    SPDR = 0;
}

// transfer sessions consist of three bytes: 
// first a trash byte is exchanged
// then int16_t milliamps values come and go
ISR(SPI_STC_vect) {
    trash = SPDR;
    for (uint8_t i=0; i<2; i++) {
        SPDR = ((uint8_t *)&spi_out_msg)[i];
        while (!(SPSR & _BV(SPIF)));
        ((uint8_t *)&spi_in_msg)[i] = SPDR;
    }
    SPDR = trash;
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

int main(void) {
    spi_slave_init();
    adc_init();
    pwm_init();

    OCR0A = OCR0B = 255; // set PWM for 0% duty cycle
    int32_t tmp = 0;
    for (int16_t i=0; i<1024; i++) {
        _delay_ms(1);
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            tmp += adc0;
        }
    }
    int16_t acs714_zero = tmp >> 10;


    DDRD |= _BV(PD7); // пин PD7 будет мигать, это нужно для контроля длительности одной итерации цикла (осциллограф в помощь)
    DDRD |= _BV(PD4); // пин PD4 будет гореть, означая, что контроллер готов к работе
    uint8_t pd7_pin_toggle = 1; 
    PORTD |= _BV(PD4);

    int16_t g = 0, pwm = 0, y = 0;
    while (1) {
        pd7_pin_toggle = 1 - pd7_pin_toggle;
        if (pd7_pin_toggle) {
            PORTD |=  _BV(PD7);
        } else {
            PORTD &= ~_BV(PD7);
        }

        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            y = adc0;
        }

        // АЦП adc0_val бегает от 0 до 1023, при этом ноль ампер соответствует 1023/2 АЦП
        // поэтому сначала преобразуем adc0_val в значение от -512 до 512, с нулём ампер в нуле: adc0_val - 512
        // чтобы получить амперы из выражения (adc0_val - 512), надо его умножить на (1/512)*(2.5V/.185V/A)
        y = (y - acs714_zero)*26; // y измеряется в миллиамперах, 2.5 / (.000185 * 512) = 26.4

        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            g = spi_in_msg; // g тоже должен измеряться в миллиамперах
            spi_out_msg = y;
        }

        // насытим до 97% выходного максимума
        if (g> 5000L) g =  5000L;
        if (g<-5000L) g = -5000L;

        pwm = (g*255L)/5000L;

        if (pwm>=0) {
            OCR0A = 255;
            OCR0B = 255-pwm;
        } else {
            OCR0A = 255+pwm;
            OCR0B = 255;
        }
    }

    PORTD &= ~_BV(PD4);
    return 0;
}

