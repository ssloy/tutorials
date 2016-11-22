#define F_CPU 16000000L

#include <avr/io.h>
#include <util/atomic.h>
#include <util/delay.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volatile uint16_t adc0_val  = 512;   // adc readings
volatile uint8_t  spi_in_msg  = 128; // задание силы тока
volatile uint8_t  spi_out_msg = 0;   // реально протекающий ток
volatile uint8_t tot_overflow;       // global variable to count the number of overflows of timer1

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void spi_slave_init() {
    DDRB &= ~(_BV(PB3) | _BV(PB5)  | _BV(PB2)); // set MOSI, SCK and SS as input
    DDRB |=   _BV(PB4);                         // set MISO as output
    SPCR  =   _BV(SPE) | _BV(SPR0) | _BV(SPIE); // enable slave SPI at clock rate Fck/4, turn on interrupts
}

ISR(SPI_STC_vect) {
    spi_in_msg = SPDR;
    SPDR = spi_out_msg;
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
    adc0_val = ADCL | (ADCH<<8);  // Must read low first
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
    spi_slave_init();
    adc_init();
    pwm_init();

    OCR0A = OCR0B = 255; // set PWM for 0% duty cycle
    int32_t tmp = 0;
    for (int16_t i=0; i<1024; i++) {
        _delay_ms(1);
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            tmp += adc0_val;
        }
    }
    int16_t acs714_zero = tmp >> 10;

    int16_t g=0, y=0;
    int32_t e=0, x=0, u=0;
    int16_t pwm=0;
    uint32_t micros=0, oldmicros=0, dt=0;

    // Вход: y (измерение силы тока), g (задание по току)
    // Внутренняя переменная x.
    // Выход: u (напряжение)
    //
    // определим ошибку как e( k ):
    // e ( k ) = g( k ) - y( k );
    //
    // интегрируем ошибку, ограничивая интеграл на каждом шаге суммирования:
    // x ( k ) = min(x_max, max(x_min, x ( k-1 ) + Ki*dt*e( k )));
    //
    // тогда выход это взвешенная сумма ошибки и интеграла ошибки:
    // u ( k ) = min(u_max, max(u_min, x ( k ) + Kp*e( k )));

    // АЦП adc0_val бегает от 0 до 1023, при этом ноль ампер соответствует 1023/2 АЦП
    // поэтому сначала преобразуем adc0_val в значение от -512 до 512, с нулём ампер в нуле:  adc0_val - 512
    // затем y измеряется в миллиамперах : [-512,512]/512 * (2.5V / .000185V/mA)

    DDRD |= _BV(PD7);        // пин PD7 будет мигать, это нужно для контроля длительности одной итерации цикла (осциллограф в помощь)
    int pd7_pin_toggle = 1; 

    timer_init();
    while (1) {
        pd7_pin_toggle = 1 - pd7_pin_toggle;
        if (pd7_pin_toggle) {
            PORTD |=  _BV(PD7);
        } else {
            PORTD &= ~_BV(PD7);
        }

        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            y = adc0_val;
        }
        spi_out_msg = y >> 2;
        y = (y - acs714_zero)*26;

        g = ((int32_t)spi_in_msg<<1) - 256;
        // g тоже должен измеряться в миллиамперах, поэтому домножим его на 5000/255 = 19.6, округляем до 20
        g *= 20L;

        // ошибка - это разница между желаемым и реальным током
        e = g - y;

        // x - это интеграл ошибки (с насыщением), помноженной на коэффициент Ki (у нас = 2268)
        // x = x + Ki * dt * e
        // давайте предположим 130 микросекунд на цикл интегратора, тогда dt = .000130 (в реальности надо измерить длину цикла)
        // тогда x = x + e * 2268 * .000130
        // чтобы не делать вычислений с плавающей точкой, давайте делать вычисления с фиксированной точкой, помножим этот x ещё и на 1000
        // тогда x = x + e * (2268*.130), максимальная ошибка у нас ограничена +- десятью тысячами миллиампер (+-5А максимальный ток * 2)
        // e*2268*130 вылезает за int32_t, поэтому вместо (e*2268*dt)/1000 напишем (e*(2268/4)*dt)/250
        oldmicros = micros;
        ATOMIC_BLOCK(ATOMIC_FORCEON) {
            micros = (TCNT1 + ((uint32_t)tot_overflow<<16)) << 6;
        }
        dt = micros - oldmicros; // TODO uint32_t overflow (через час работы контроллер может сойти с ума)
        x = x + ((e*567L*(int32_t)dt)/250L);

        // теперь насытим это дело до 80% от выходного значения. Максимум выходного значения +-24V, тогда 80% это +- 19.2V
        // x, по факту, измеряется в микровольтах
        if (x> 19200000L) x =  19200000L;
        if (x<-19200000L) x = -19200000L;

        // пропорциональный коэффициент Kp=3, но фиксированная точка требует ещё множителя 1000, поэтому 3000L*e
        u = 3000L*e + x; // микровольты

        // насытим до 97% выходного максимума
        if (u> 23280000L) u =  23280000L;
        if (u<-23280000L) u = -23280000L;

        // итого пересчитаем скважность шим из микровольт: pwm = u*255/24000000
        // если напрямую умножать на 255, то будет переполнение int32_t, поэтому сократим дробь до 255/24000000 51/4800000
        // максимальный u это 24*10^6, если умножить на 51, то оно ещё (со скрипом) влезает в int32_t

        pwm = (u*51L)/4800000L;

        if (pwm>=0) {
            OCR0A = 255-pwm;
            OCR0B = 255;
        } else {
            OCR0A = 255;
            OCR0B = 255+pwm;
        }
    }
    return 0;
}

