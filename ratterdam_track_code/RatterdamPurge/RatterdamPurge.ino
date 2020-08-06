#include <Wire.h>

int targetAlley = 24;
void setup() {
  Wire.begin();
  Serial.begin(9600);

}

void loop() {
  for(int i=21; i<25; i++){
    Serial.println(i);
    for(int r=0; r < 5; r++){
      Serial.println(r);
  Wire.beginTransmission(i);
  Wire.write(0);
  int o = Wire.endTransmission();
  //Serial.println(o);
  delay(20000);
    }
    delay(10000);
  }


//  Wire.beginTransmission(targetAlley);
//  Wire.write(0);
//  int o = Wire.endTransmission();
//  Serial.println(o);
//  delay(60000);
}
