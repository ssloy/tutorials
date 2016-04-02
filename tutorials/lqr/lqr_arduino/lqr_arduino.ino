volatile long xi = 0, xi_1 = 0;
int vi = 0;

volatile char ABprev = 0;
const int increment[16] = {0,-1,1,0, 1,0,0,-1, -1,0,0,1, 0,1,-1,0};

ISR (PCINT0_vect) { // D8 or D9 has changed
    char AB = PINB & 3;
    xi += increment[AB+ABprev*4];
    ABprev = AB;
}

void set_speed(int speed) {
    if (speed<0) {
      analogWrite(3, 0);
      analogWrite(5, speed<-255?-255:-speed);
    } else if (speed>0) {
      analogWrite(3, speed>255?255:speed);
      analogWrite(5, 0);
    } else {
      analogWrite(3, 0);
      analogWrite(5, 0);
    }
}

void setup() {
    pinMode(8, INPUT);  // A
    pinMode(9, INPUT);  // B
    pinMode(3, OUTPUT);  // LPWM
    pinMode(5, OUTPUT);  // RPWM

    set_speed(64); // move the cart slowly to the rightmost position
    delay(4000);
    set_speed(0);
    xi_1 = xi = 24500;
    vi = 0;
    
    ABprev = 0;
    PCICR |= (1 << PCIE0);  // interrupt will be fired on any change on pins d8 and d9
    PCMSK0 |= 3;
}

void loop() {
    vi = xi-xi_1;
    int ui = 255.*(-0.000973669*xi -0.0563218*vi)/24.;
    xi_1 = xi;
    set_speed(ui);
    delay(2);
}

