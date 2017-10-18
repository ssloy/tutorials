volatile long xi = 0;
volatile char ABprev = 0;
const int increment[16] = {0,-1,1,0, 1,0,0,-1, -1,0,0,1, 0,1,-1,0};

ISR (PCINT0_vect) { // D8 or D9 has changed
  char AB = PINB & 3;
  xi += increment[AB+ABprev*4];
  ABprev = AB;
}

void setup() {
  pinMode(8, INPUT);  // A
  pinMode(9, INPUT);  // B
  PCICR |= (1 << PCIE0);  // interrupt will be fired on any change on pins d8 and d9
  PCMSK0 |= 3;
  Serial.begin(115200);
}

void loop() {
  Serial.print("ATmega328p: ");
  Serial.println(xi);
  delay(100);
}

