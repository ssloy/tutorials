volatile long cnt = 0;
volatile char ABprev = 0;
const int increment[16] = {0,-1,1,0, 1,0,0,-1, -1,0,0,1, 0,1,-1,0};

int nb_print;
int prev;
int result[600];

ISR (PCINT0_vect) { // one of pins D8 to D10 has changed
    char AB = PINB & 3;
    cnt += increment[AB+ABprev*4];
    ABprev = AB;
}


void go_right(int speed) {
    analogWrite(3, speed);
    analogWrite(5, 0);
}
void go_left(int speed) {
    analogWrite(3, 0);
    analogWrite(5, speed);
}
void brake() {
    analogWrite(3, 0);
    analogWrite(5, 0);
}

unsigned long start;



void setup() {
    Serial.begin(115200);
    pinMode(8, INPUT);  // A
    pinMode(9, INPUT);  // B
    pinMode(10, INPUT); // Z
    pinMode(3, OUTPUT);  // LPWM
    pinMode(5, OUTPUT);  // RPWM

    go_right(64);
    delay(4000);
    brake();

    prev = cnt = 24500;
    ABprev = 0;

    nb_print = -1;


    PCICR |= (1 << PCIE0);
    PCMSK0 |= 7;
    go_left(96);

    for (int i=0; i<600; i++) {
        result[i] = 0;
    }

    start = millis();
}

void loop() {
    if (cnt!=prev && nb_print==-1) nb_print = 0;

    if (nb_print!=-1 && nb_print<600) {
        result[nb_print] = cnt-prev;
    }
    if (nb_print<2000) {
        nb_print++;
        int speed = 255.*cnt/20000;
        if (speed>255) speed = 255;
        go_left(speed);
    } else {
        brake();
        Serial.println(nb_print);
        Serial.println(cnt);
        Serial.println(millis()-start);
        nb_print = -1;
    }
    
    /*
    if (nb_print==600 || cnt>40000) {
        brake();
        for (int i=0; i<nb_print; i++) {
            Serial.print(" ");
            Serial.print(result[i]);
        }
        cnt=0;
        nb_print = 601;
        Serial.println(" ");
    }
    */
    prev = cnt;
    delay(2);
}

