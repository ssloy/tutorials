GCCFLAGS=-g -Os -Wall -mmcu=atmega328p -std=gnu99
LINKFLAGS=-Wl,-u,vfprintf -lprintf_flt -Wl,-u,vfscanf -lscanf_flt -lm
AVRDUDEFLAGS= -c usbasp -P usb -p m328p -B 5

all:    pwm-upload

pwm.hex:  pwm.c
	avr-gcc ${GCCFLAGS} ${LINKFLAGS} -o pwm.o pwm.c
	avr-objcopy -j .text -O ihex pwm.o pwm.hex

pwm.ass:  pwm.hex
	avr-objdump -S -d pwm.o > pwm.ass

pwm-upload:       pwm.hex
	sudo avrdude ${AVRDUDEFLAGS} -U flash:w:pwm.hex:a

