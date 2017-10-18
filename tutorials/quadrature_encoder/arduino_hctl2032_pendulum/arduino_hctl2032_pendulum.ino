const int SHIFT_PL_PIN = 11;
const int SHIFT_CP_PIN = 12;
const int SHIFT_Q7_PIN = 13;

const int DECODER_CHIX_PIN  = 9;

const int DECODER_OE_PIN   = 8;
const int DECODER_XY_PIN   = 2;
const int DECODER_SEL1_PIN = 7;
const int DECODER_SEL2_PIN = 3;
const int DECODER_RST_PIN  = 4;

const int MOTOR_PWM = 5;
const int MOTOR_DIR = 6;

void set_speed(int speed) {
   analogWrite(MOTOR_PWM, abs(speed) & 0xFF); // & 0xFF is a cheap clamp 0-255
   digitalWrite(MOTOR_DIR, speed<0);
}

const unsigned char avago_byte_select[4][2] = {{1,0},{0,0},{1,1},{0,1}};
inline void decoder_latch_byte(unsigned char i) {
  digitalWrite(DECODER_SEL1_PIN, avago_byte_select[i][0]);
  digitalWrite(DECODER_SEL2_PIN, avago_byte_select[i][1]);
}

unsigned char read_shift_register() {
  unsigned char result = 0;
  digitalWrite(SHIFT_PL_PIN, 0);
  digitalWrite(SHIFT_PL_PIN, 1);
  for (unsigned char i=8; i--;) {
    result |= (digitalRead(SHIFT_Q7_PIN)<<(i));
    digitalWrite(SHIFT_CP_PIN, 1);
    digitalWrite(SHIFT_CP_PIN, 0);
  }
  return result;
}

long read_decoder(unsigned char encoder) {
  long result = 0;
  digitalWrite(DECODER_XY_PIN, encoder);
  digitalWrite(DECODER_OE_PIN, 0);
  delayMicroseconds(10); // The OE/, SEL1, and SEL2 inputs are sampled by the internal inhibit logic on the falling edge of the clock
  for (unsigned char i=4; i--;) {
    decoder_latch_byte(i);
    result += read_shift_register();
    if (i) result <<= 8;
  }
  digitalWrite(DECODER_OE_PIN, 1);
  return result;
}

void decoder_reset() {
  digitalWrite(DECODER_RST_PIN, 0);
  delayMicroseconds(10);
  digitalWrite(DECODER_RST_PIN, 1);
}

void setup() {
  pinMode(MOTOR_PWM, OUTPUT);
  pinMode(MOTOR_DIR, OUTPUT);

  pinMode(SHIFT_PL_PIN, OUTPUT);
  digitalWrite(SHIFT_PL_PIN,   1);

  pinMode(SHIFT_CP_PIN, OUTPUT);
  digitalWrite(SHIFT_CP_PIN,   0);

  pinMode(SHIFT_Q7_PIN, INPUT);

  pinMode(DECODER_OE_PIN, OUTPUT);
  digitalWrite(DECODER_OE_PIN, 1);

  pinMode(DECODER_RST_PIN, OUTPUT);
  digitalWrite(DECODER_RST_PIN, 1);

  pinMode(DECODER_XY_PIN, OUTPUT);
  digitalWrite(DECODER_XY_PIN, 0); // X encoder

  pinMode(DECODER_SEL1_PIN, OUTPUT);
  pinMode(DECODER_SEL2_PIN, OUTPUT);
  decoder_latch_byte(0);
  
  pinMode(DECODER_CHIX_PIN, INPUT);  
//  PCICR  |= (1 << PCIE0);  // interrupt will be fired on any change on pin d9
//  PCMSK0 |= 2;

  Serial.begin(115200);
  decoder_reset();
}

volatile long offset = 4000; // index pulse position in ticks
volatile long ticks  = 0;    // current angle measurement, in ticks

inline double read_cart_position() { // the encoder is supposed to be reset at the zero (center) position
    return read_decoder(1)*(.08/4000.); // .08 m/rev[40 GT2 teeth] / 1000 pulses/rev
}

ISR (PCINT0_vect) { // D9 has changed, thus index pulse fired, updating its position
    offset = ticks;
}

inline double read_pendulum_angle() {
    cli();
    ticks = read_decoder(0);
    long tticks = ticks - offset; // if index pulse position off the balance, correct it here
    sei();
    
    while (tticks> 4000) tticks -= 8000; // mod 2pi
    while (tticks<-4000) tticks += 8000;
    
    return tticks*(3.1416*2./8000.); // 2pi / 2000 pulses/rev
}

int state = 0;
long time_prev;

long cnt = 0;
const int memlen = 2;
double X[memlen], TH[memlen];

long tticks,tticks_prev;
void loop() {
  tticks = read_decoder(0);
  long dt = micros() - time_prev;
  time_prev = micros();

  Serial.print("HCTL2032: ");
  Serial.print(tticks);
  Serial.print(", ");
  Serial.print((tticks-tticks_prev)/8000.*1000000./dt);
  Serial.println(" RPM");
  tticks_prev = tticks;

  delay(100);
  return;
  
  if (0==state) {
    double x  = read_cart_position();
    double th = read_pendulum_angle();

    Serial.print(offset);
    Serial.print("\t");
    Serial.print(read_decoder(0));
    Serial.print("\t");
    Serial.print(th/3.1416*180);
    Serial.print("\t");
    Serial.println(x);

    if (fabs(th)<.1) {
       state = 1;
       return;
    }
  }
  if (1==state) {
    long dt = micros() - time_prev;
    time_prev = micros();
    double x  = read_cart_position();
    double th = read_pendulum_angle();
    double v = (x  -  X[cnt%memlen])/(.02*memlen);
    double w = (th - TH[cnt%memlen])/(.02*memlen);
    X [cnt%memlen] = x;
    TH[cnt%memlen] = th;
    cnt++;

    if ((fabs(x)>.22 && x*v>0) || fabs(th)>.3) {
        set_speed(0);
        state = 2;
        return;
    }
    
    double c[] = {34, 25, 89, 19};

    double f = c[0]*x + c[1]*v + c[2]*th + c[3]*w;
    double u = ((1-.404)*v + .02*f/.49 + copysignf(.01, v))/.055;
    if (cnt>memlen) set_speed(u*255./24.);

    Serial.print(offset);
    Serial.print("\t");
    Serial.print(x);
    Serial.print("\t");
    Serial.print(v);
    Serial.print("\t");
    Serial.print(th/3.1416*180);
    Serial.print("\t");
    Serial.print(w/3.1416*180);
    Serial.print("\t");
    Serial.print(u);
    Serial.print("\t");
    Serial.println(dt);
  }
  delay(20-4); // 3ms for printing
  if (2==state) {
//    Serial.print(read_cart_position());
//    Serial.print("\t");
//    Serial.println(read_pendulum_angle()/3.1416*180);
    set_speed(copysignf(20,-read_cart_position()));
    while (fabs(read_cart_position())>.01) delay(10);
    set_speed(0);
    state = 3;
  }
}

