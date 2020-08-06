#include <Wire.h>
#include <Keypad.h>


// Keypad Initialization - Start
const byte numRows= 4; //number of rows on the keypad
const byte numCols= 4; //number of columns on the keypad
//keymap defines the key pressed according to the row and columns just as appears on the keypad
char keymap[numRows][numCols]= 
{
{'R', '1', '1', 'X'}, 
{'P', '1', '1', 'T'},   
{'1', '1', '1', '1'},
{'Q', '1', '1', 'S'}
};
//Code that shows the the keypad connections to the arduino terminals
byte rowPins[numRows] = {46,48,50,52}; 
byte colPins[numCols]= {36,38,40,42}; 
//initializes an instance of the Keypad class
Keypad myKeypad= Keypad(makeKeymap(keymap), rowPins, colPins, numRows, numCols);
// Keypad Initialization - End


//Serial comm managment decl start
const byte numBytes = 3;
char data[numBytes];
int v;
int alleyCode;
// 3 byte command from python over Serial:
int dir; //
int alley; // 
int command; //
int input;
int buzzer = 10;
//Serial comm management delc end

void setup() {
  Wire.begin();
  Serial.begin(115200);
  Serial.print("Ready");

	//pinMode(buzzer,OUTPUT);
}

void loop() {

    //Read instructions from master arduino
    v = Serial.readBytes(data, numBytes);
    dir = data[0]; // direction I/O of comm.
    alley = data[1]; // alley code
    command = data[2]; // command. NB X means nothing, used for receiving data not sending


    if(dir == 0){
    Wire.requestFrom(alley, 1);
    input = Wire.read();
    Serial.println(input);
    }

    if(dir == 3){
      char keypressed = myKeypad.getKey();
    if(keypressed == 'S'){
      Serial.print('S');
    }
    else if(keypressed == 'P'){
      Serial.print('P');
    }
		else if(keypressed == 'R'){
			Serial.print('R');
		}
   else if(keypressed == 'X'){
    Serial.print('X');
   }
   else if(keypressed == 'T'){
    Serial.print('T');
   }
   else if(keypressed == 'Q'){
    Serial.print('Q');
   }
    else{
      Serial.println(3);
    }  
    }
    
    if(dir == 1){
       Wire.beginTransmission(alley);
       Wire.write(command);
       Wire.endTransmission(); 
       Serial.println(5); //throwaway -py expects something
    }


}
