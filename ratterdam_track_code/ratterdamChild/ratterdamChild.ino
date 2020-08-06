  #include <Wire.h>

int childAddress = 18;
int idx = childAddress - 8;
int sensorAvalues[] = {780, 500,800, 1023, 1023,0,1023,1023,0,1023,450,0,0,0,0,1023,500};
int sensorBvalues[] = {1023, 900, 1023, 900, 670,0,1023,980,0,800,950,0,0,0,0,1023,600};
int pd1 = A1;
int pd2 = A0;
int motor = 3;
int led = 4;
int buzzer = 8;
int purgePulses = 39;
int BeamState = 0; // status of PDs - 0 None; 1 A; 2 B; 3 Both
int RewardState = 0; //active status of alley - 0 not; 1 yes rewarded 

int pdASensor;
int pdBSensor;

int pdABreak = 0;
int pdBBreak = 0;


void setup() {
  //Serial.begin(115200);
  Wire.begin(childAddress);
  pinMode(pd1, INPUT);
  pinMode(pd2, INPUT);

	pinMode(buzzer, OUTPUT);
  pinMode(motor, OUTPUT);
  pinMode(led, OUTPUT);

		 Wire.onReceive(processInstruction);
		 Wire.onRequest(sendInformation);
}

void loop() {
foragingV2(); 
//testParts();
// testLogic();
//testI2C();
}

void sendInformation(){
  Wire.write(BeamState);
}


void processInstruction(){
  if(Wire.available()){
  int c = Wire.read();
  RewardState = c;
  }
  }

void foragingV2(){
  pdASensor = analogRead(pd1);
  pdASensor = map(pdASensor, 0,  sensorAvalues[idx], 0, 1023);
  if(pdASensor < 512){
    pdASensor = 0;
  }
  else if(pdASensor >= 512){
    pdASensor = 1;
  }

  pdBSensor = analogRead(pd2);
  pdBSensor = map(pdBSensor, 0, sensorBvalues[idx], 0, 1023);
  if(pdBSensor < 512){
    pdBSensor = 0;
  }
  else if(pdBSensor >= 512){
    pdBSensor = 1;
  }

// Neither broken
	if(pdASensor == 0 && pdBSensor == 0){
		BeamState = 0;
      if(RewardState == 1){
      reward();
      RewardState = 0;
    }
	}

	// A broken
	else if(pdASensor == 0 && pdBSensor == 1){
		BeamState = 1;
		if(RewardState == 1){
			reward();
			RewardState = 0;
		}
	}

	
	else if(pdASensor == 1 && pdBSensor == 0){

		 BeamState = 2;
		 if(RewardState == 1){
			reward();
			RewardState = 0;
		}
	}

	
	else if(pdASensor == 1 && pdBSensor == 1){

		BeamState = 3;

	}

 if(RewardState == 5){
  digitalWrite(motor, HIGH);
  delay(10);
 }

 if(RewardState == 9){
  digitalWrite(motor, LOW);
  delay(10);
 }

 if(RewardState == 7){
  if(pdASensor == 0){
    for(int i = 0; i<6;i++){
      digitalWrite(led, HIGH);
      delay(50);
      digitalWrite(led,LOW);
    }
  }
  if(pdBSensor == 0){
   for(int i = 0; i<6;i++){
      digitalWrite(led, HIGH);
      delay(200);
      digitalWrite(led,LOW);
    }
  }
  
 }

 if(RewardState == 6){
  digitalWrite(motor,HIGH);
  delay(500);
  digitalWrite(motor,LOW);
 }

}

void reward(){
	digitalWrite(motor,HIGH);
	delay(500);
	digitalWrite(motor,LOW);
	tone(buzzer,500,750); 
  digitalWrite(led,HIGH);
  delay(2000);
  digitalWrite(led,LOW);
  }

void error(){
	tone(buzzer, 1000, 250);
  for(int i=0; i<5; i++){
    digitalWrite(led,HIGH);
    delay(100);
    digitalWrite(led,LOW);
  }
}

void purge(){
   for(int x = 0; x < purgePulses; x++){
	 digitalWrite(motor,HIGH);
	 delay(500);
	 digitalWrite(motor,LOW);
   }
 }

void testParts(){

  for(int i=0; i<11; i++){
	pdASensor = analogRead(pd1);
	pdBSensor = analogRead(pd2);
	Serial.print("PD A: ");
	Serial.print(pdASensor);
	Serial.print("\n");
  Serial.print("PD B: ");
  Serial.print(pdBSensor);
  Serial.print("\n");
  delay(100);
  }

  delay(100);
  for(int i=0; i<4; i++){
  digitalWrite(led, HIGH);
  delay(100);
  digitalWrite(led, LOW);
  }
  delay(100);
  for(int i =0; i<5; i++){
  digitalWrite(motor,HIGH);
  delay(100);
  digitalWrite(motor,LOW);
  delay(500);
  }
  delay(100);

  tone(buzzer,500,300);
}


void testI2C(){
if(RewardState == 1){
  digitalWrite(led,HIGH);
  delay(75);
  digitalWrite(led,LOW);

  digitalWrite(motor, HIGH);
  delay(150);
  digitalWrite(motor, LOW);

  tone(buzzer, 500, 100);
  RewardState = 0;
}
  
}

