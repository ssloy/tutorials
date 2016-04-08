volatile long xi = 0, xi_1 = 0; // encoder ticks
int vi = 0;                     // encoder ticks per 2ms
double zi = 0;

volatile char ABprev = 0;
const int increment[16] = {0,-1,1,0, 1,0,0,-1, -1,0,0,1, 0,1,-1,0};

int cnt;
unsigned long start_time;
long xi_log[256];
int  vi_log[256];

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
    Serial.begin(115200);

    pinMode(8, INPUT);  // A
    pinMode(9, INPUT);  // B
    pinMode(3, OUTPUT);  // LPWM
    pinMode(5, OUTPUT);  // RPWM

    set_speed(32); // move the cart slowly to the rightmost position
    delay(6000);
    set_speed(-32);
    delay(100);
    set_speed(0);
    delay(500);

    xi_1 = xi = 24500;
    vi = 0;

    ABprev = 0;
    PCICR  |= (1 << PCIE0);  // interrupt will be fired on any change
    PCMSK0 |= 3;             // on pins d8 and d9

    cnt = 0;
    start_time = millis();
}

void loop() {
    vi = xi-xi_1;
    zi += xi*.001*.01*.002;
    int ui = 255.*(-0.000973669*xi -0.0563218*vi)/24.;

    if (1024>cnt) {
        if (0==cnt&3) {
            xi_log[cnt>>2] = xi;
            vi_log[cnt>>2] = vi;
        }
        cnt++;
        set_speed(ui);
    } else if (1024==cnt) {
        set_speed(0);
        cnt++;
        Serial.print(millis()-start_time);
        Serial.println("ms ");

        for (int i=0; i<256; i++) {
            Serial.print(xi_log[i]*.001*.01); // meters
            Serial.print(" ");
        }
        Serial.println("");

        for (int i=0; i<256; i++) {
            Serial.print(vi_log[i]*.005); // meters/second
            Serial.print(" ");
        }
        Serial.println("");
    }

    xi_1 = xi;
    delay(2);
}

