const int SHIFT_PL_PIN = 11;
const int SHIFT_CP_PIN = 12;
const int SHIFT_Q7_PIN = 13;

const int DECODER_OE_PIN   = 8;
const int DECODER_XY_PIN   = 2;
const int DECODER_SEL1_PIN = 7;
const int DECODER_SEL2_PIN = 3;
const int DECODER_RST_PIN  = 4;

const int MOTOR_PWM = 5;
const int MOTOR_DIR = 6;


void set_speed(int speed) {
   analogWrite(MOTOR_PWM, abs(speed));
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

double cart_x_prev = 0;
long start_time, time_prev;

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

  Serial.begin(115200);
}

int voltage = -255;
int state = 0;

void loop() {
  if (0==state) { // retour chariot
    Serial.print("float data_");
    Serial.print(voltage);
    Serial.print("[] = {");
    Serial.print(voltage/255.*24.);
    Serial.print(",");
    state = 1;

    set_speed(voltage<0 ? 20 : -20);
    delay(6000);
    set_speed(voltage<0 ? -20 : 20);
    delay(100);
    set_speed(0);
    delay(500);
    decoder_reset();

    set_speed(voltage);

    cart_x_prev = 0;
    start_time = micros();
    time_prev = 0;
  }
  if (1==state) {
    long t  = micros() - start_time;
    long dt = t - time_prev;
    time_prev = t;
    double cart_x = read_decoder(1)*.00008; // *1000 pulses/rev * .08 m/rev[40 GT2 teeth]
    cart_x_prev = cart_x;

    if (fabs(cart_x)<.3 && t<2000000) {
      Serial.print(cart_x*1000);
      Serial.print(",");
      Serial.print(dt);
      Serial.print(",");
    } else {
      set_speed(0);
      state = 2;
    }
  }
  if (2==state) {
    if (voltage+10<256) {
      Serial.println("};");
      voltage += 10;
      state = 0;
    }
  }
}

