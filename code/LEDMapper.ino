/*
  LEDMapper

  Controls an array of LEDs in common-row cathode configuration

*/
#define BLINK_DELAY_MS 1000
#define NLEDs 36
#define Nrows 6
#define Ncols 6
#define Npins 12
#define Pmax 256

const int chgPin = 12;
const int measPin = 7;
const int PinMapping[] = {2,3,4,5,6,7,8,9,10,11,12,44};

void getrowcol(const int index, int rowcol[]){
  rowcol[0] = index / Ncols;
  rowcol[1] = index % Ncols;
}

// the setup function runs once when you press reset or power the board
void setup() {
  // Need 12 digital pins for output. Initialize all of them for output and ground them
  for (int pin = 0; pin < Npins; ++pin){
    pinMode(pin, OUTPUT);
    digitalWrite(pin, LOW);
  }
}

// the loop function runs over and over again forever
void loop() {
  for (int power = 0; power < Pmax; ++power){
    for (int led = 0; led < NLEDs; ++led){
        int rowcol[] = {0, 0};
        getrowcol(led, rowcol);
        printf("Row %d Col %d.", rowcol[0], rowcol[1]);
        analogWrite(PinMapping[rowcol[1]], power);
        delay(BLINK_DELAY_MS);
        analogWrite(PinMapping[rowcol[1]], 0);
    }
  }
  for (int power = Pmax; power >= 0; --power){
    for (int led = 0; led < NLEDs; ++led){
        int rowcol[] = {0, 0};
        getrowcol(led, rowcol);
        printf("Row %d Col %d.", rowcol[0], rowcol[1]);
        analogWrite(PinMapping[rowcol[1]], power);
        delay(BLINK_DELAY_MS);
        analogWrite(PinMapping[rowcol[1]], 0);
    }
  }

}
