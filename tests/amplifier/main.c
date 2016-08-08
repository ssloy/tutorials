#define F_CPU 16000000L

#include <avr/io.h>
#include <util/atomic.h>
#include <util/delay.h>
#include <avr/pgmspace.h>

volatile uint8_t analog_val = 127; // adc readings
volatile uint8_t tot_overflow;     // global variable to count the number of overflows of timer1

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void spi_master_init(void) {
    DDRB  |= (1<<PB3)|(1<<PB5)|(1<<PB2);   // set MOSI, SCK and SS as output
    DDRB  &= ~(1<<PB4);                    // set MISO as input
    PORTB |= (1<<PB2);                     // set SS to high
    SPCR   = (1<<SPE)|(1<<MSTR)|(1<<SPR0); // enable master SPI at clock rate Fck/8
    SPSR   = (1<<SPI2X);
}

uint8_t spi_tranceiver(uint8_t data) {
    PORTB &= ~(1<<PB2); // select slave
    SPDR = data;        // send data
    while (!(SPSR &(1<<SPIF))); // wait for transmition complete
    PORTB |= (1<<PB2);  // release slave
    return SPDR;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void adc_init() {
    // Bits 7:6    01 - Set REFS1..0 in ADMUX to AREF, to change reference voltage to the proper source
    // Bits 7:6    00 - Set REFS1..0 in ADMUX to AREF, internal Vref turned off
    // Bit  5       0 - clear ADLAR in ADMUX (0x7C) to right-adjust the result ADCL will contain lower 8 bits, ADCH upper 2 (in last two bits)
    // Bit  4       0 - unuseed
    // Bits 3:0  0000 - ADC0 (Table 24-4)
    ADMUX = 0b00000000;

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

ISR(ADC_vect) { // Interrupt service routine for the ADC completion
    analog_val = (ADCL>>2) | (ADCH<<6);  // Must read low first
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void timer_init() {
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

ISR(TIMER1_OVF_vect) { // TIMER1 overflow interrupt service routine called whenever TCNT1 overflows
    tot_overflow++; // keep a track of number of overflows
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pwm_init() {
    DDRD |= (1 << DDD6) | (1 << DDD5); // PD6 and PD5 are now output

    OCR0A = OCR0B = 255; // set PWM for 0% duty cycle

    TCCR0A  = (1 << COM0A1) | (1 << COM0A0) | (1 << COM0B1) | (1 << COM0B0); // set inverting mode  to avoid spikes for zero duty cycles
    TCCR0A |= (1 << WGM01) | (1 << WGM00);   // set fast PWM Mode
    TCCR0B |= (1 << CS01);                   // set prescaler to 8
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
    spi_master_init();
    adc_init();
    timer_init();
    pwm_init();

    int32_t tmp = 0;
    for (uint16_t i=0; i<512; i++) {
        tmp += analog_val;
        _delay_ms(1);
    }
    int16_t acs714_zero = (int16_t)(tmp*2./512.); // calibrate current sensor zero


    for (uint8_t measurement=0; measurement<11; measurement++) {
        uint32_t micros;
        uint16_t milbis;
        int16_t  g, y;

        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            TCNT1 = 0;
            tot_overflow = 0;
        }
        while (1) {
            ATOMIC_BLOCK(ATOMIC_FORCEON) {
                micros = (TCNT1 + ((unsigned long)tot_overflow<<16)) << 6;
            }
            milbis = micros >> 10;
            if (milbis>=500) break;

            uint16_t sine_idx = 0;
            switch (measurement) {
                case 0: sine_idx = (micros>>5) & 0x1FF; break; // sine 61.04 Hz : 32uS per increment, one sine period in 32uS*512 samples, thus the frequency is 1/(32*512/10^6) = 61.0351
                case 1: sine_idx = (micros>>6) & 0x1FF; break; // sine 30.52 Hz
                case 2: sine_idx = (micros>>7) & 0x1FF; break; // sine 15.26 Hz
            }
            if (measurement<3) {
                if (sine_idx & 0x100) {
                    g = -pgm_read_byte(&sinewave_data[sine_idx & 0xFF]);
                } else {
                    g =  pgm_read_byte(&sinewave_data[sine_idx & 0xFF]);
                }
            } else {
                switch (measurement) {
                    case  3: g = milbis & 64 ? 0 :  204; break;
                    case  4: g = milbis & 64 ? 0 : -204; break;
                    case  5: g = milbis & 64 ? 0 :  153; break;
                    case  6: g = milbis & 64 ? 0 : -153; break;
                    case  7: g = milbis & 64 ? 0 :   77; break;
                    case  8: g = milbis & 64 ? 0 :  -77; break;
                    case  9: g = milbis & 64 ? 0 :   26; break;
                    case 10: g = milbis & 64 ? 0 :  -26; break;
                    case 11: g = 0; break;
                }
            }

            if (g>=0) {
                OCR0A = 255-g;
                OCR0B = 255;
            } else {
                OCR0A = 255;
                OCR0B = 255+g;
            }

            spi_tranceiver(*((uint8_t *)(void *)(&g)+1));
            spi_tranceiver(*((uint8_t *)(void *)(&g)+0));

            y = ((((int16_t)(analog_val))<<1) - acs714_zero);
            spi_tranceiver(*((uint8_t *)(void *)(&y)+1));
            spi_tranceiver(*((uint8_t *)(void *)(&y)+0));
        }
        OCR0B = OCR0A = 255;
        _delay_ms(100);
    }

    return 0;
}



