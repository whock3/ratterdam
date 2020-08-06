#include <Wire.h>


int input;
int outcome;

void setup() {
  Wire.begin();
	Serial.begin(115200);
	Serial.print("Ready");

}

void loop() {

//  Wire.beginTransmission(19);
//	Wire.write(0);
//	outcome = Wire.endTransmission();
//  delay(1000);
//  Wire.requestFrom(10, 1);
//	input = Wire.read();

	for(int i=8; i < 25;i++){
	Serial.println(i);
	Wire.beginTransmission(i);
	Wire.write(0);
	outcome = Wire.endTransmission();
  Serial.println(outcome);
  delay(500);
	
	}
		

	

}
