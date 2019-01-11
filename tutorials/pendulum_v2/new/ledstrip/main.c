#define F_CPU 16000000L

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/atomic.h>
#include <util/delay.h>

#include <stdlib.h>
#include <stdio.h>
#include <avr/pgmspace.h> // PSTR

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

#define SK6812_DATA_PIN        B,0

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
//  main loop  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
void __attribute__((noinline)) led_strip_write(uint8_t *colors, uint16_t count) {
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
uint8_t red[]   = {0,128,0,0};
uint8_t green[] = {128,0,0,0};
uint8_t gray[]  = {0,0,0,128};
uint8_t blue[] = {0,0,128,0};
uint8_t max[] = {255, 255, 255, 255};

uint8_t colors[LED_COUNT*4];


int main(void) {
    OUTPUT(SK6812_DATA_PIN);
    CLEAR(SK6812_DATA_PIN);

    uart_init();
    FILE uart_stream = FDEV_SETUP_STREAM(uart_putchar, uart_getchar, _FDEV_SETUP_RW);
    stdin = stdout = &uart_stream;

    while(1) {
        for (uint8_t i=0; i<4; i++) {
            for (uint8_t j=0; j<LED_COUNT; j++) {
                uint8_t *c = max;
                if (j%5==0) {
                    c = red;
                } else if (j%5==1) {
                    c = green;
                } else if (j%5==2) {
                    c = blue;
                } else if (j%5==3) {
                    c = gray;
                } else {
                    c = max;
                }
                for (uint8_t k=0; k<4; k++) {
                    colors[((i+j)*4+k)%(LED_COUNT*4)] = c[k];
                }
            }
            led_strip_write(colors, LED_COUNT);
            _delay_ms(100);
        }
    }

    return 0;
}

