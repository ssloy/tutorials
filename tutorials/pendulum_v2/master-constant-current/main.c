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

#define AMPLIFIER_RST_PIN           B,1

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

int uart_putchar(char c, FILE *stream) {
    uart_write(c);
    return 0;
}

int uart_getchar(FILE *stream) {
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
    tot_overflow = 0;                 // initialize overflow counter variable
}

ISR(TIMER1_OVF_vect) { // TIMER1 overflow interrupt service routine called whenever TCNT1 overflows
    tot_overflow++;    // keep a track of number of overflows
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(void) {
    OUTPUT(AMPLIFIER_RST_PIN);
    CLEAR(AMPLIFIER_RST_PIN);
    _delay_ms(1);
    SET(AMPLIFIER_RST_PIN);

    _delay_ms(3000);
    timer_init();
    spi_master_init();
    decoder_init();
    shift_register_init();

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    int16_t current = 300; // mA
    int16_t ack = 0;
    /*
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
    _delay_ms(4000);
    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
    _delay_ms(1000);
*/
    decoder_reset();

    ATOMIC_BLOCK(ATOMIC_FORCEON) {
        TCNT1 = 0;
        tot_overflow = 0;
    }
    uint32_t micros=0;

    
    while (1) {
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            micros = (TCNT1 + ((uint32_t)tot_overflow<<16)) << 6;
        }
        if (micros<300000L) {
            current = +700;
        } else {
            current = 0;
        }
        spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
        if (micros>1000000L) {
            break;
        }
        int32_t ticks1 = read_decoder(0); // (36 teeth * 2mm pitch) / (1000 ticks * 4x)
        fprintf_P(&uart_stream, PSTR("%ld %ld %d %d\r\n"), micros, ticks1, current, ack);
    }

    current = 0;
    spi_transfer_sync((uint8_t *)&current, (uint8_t *)&ack, 2);
    return 0;
}

