GCCFLAGS=-g -O3 -mmcu=atmega328p -std=gnu99 -Wall -Wextra -pedantic
LINKFLAGS=-Wl,-u,vfprintf -lprintf_flt -Wl,-u,vfscanf -lscanf_flt -lm
AVRDUDEFLAGS=  -F -V -c arduino -p ATMEGA328P -P /dev/ttyUSB1 -b 57600

all:    main-upload

main.hex:  main.c
	avr-gcc ${GCCFLAGS} ${LINKFLAGS} -o main.o main.c
	avr-objcopy -O ihex -R .eeprom main.o main.hex

main.ass:  main.hex
	avr-objdump -S -d main.o > main.ass

main-upload:       main.hex
	sudo avrdude ${AVRDUDEFLAGS} -U flash:w:main.hex:a

