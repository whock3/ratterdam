/* 
Ratterdam Project
Integration with Crown IR Head Position and Bearing Tracking
WH August 1 2019


This script serves as a function generator which sends a square wave pulse
of a set duration and frequency. This pulse is sent to two machines simultaneously
1) A neuralynx cheetah recording system
2) A FlyCap FLIR camera

The pulse to the FLIR camera triggers a frame to be recorded (FLIR cam is set to operate via TTL)
The pulse to Cheeta generates a timestamp on cheetah's end to serve as synchronization record of when
in terms of the Neuralynx clock those frames were captured.

An optional step is to have the onboard LED of the Uno blink in time with the pulses.
This allows a rough sanity check of the ongoing frequency generation. If it's found to be
disruptive (take too long, interfere with pulses) then simply comment out LINES 31, 37, 40
*/

int pulsePin = 5; // This is the Arduino Uno pin which is used for the pulse.
                  // It is connected to a cable which is split with one wire going
                  // to a Neuralynx Cheetah system and another going to a FlyCap FLIR Cam
int pulse_duration_ms = 10; // This represents the duration of the pulse in milliseconds
int pulse_frequency_Hz = 30; // This represents the frequency in Hz of the pulses and therefore FLIR framerate
int inter_pulse_delay = 20;
int pinVal;
//round(((1.0/pulse_frequency_Hz)*1000)-pulse_duration_ms); // The delay between pulses.
                                                                              // It 1/frequency less the duration of the pulse

void setup() {
  Serial.begin(115200);
  pinMode(pulsePin, OUTPUT); // This sets the pulsePin to an output mode

}

void loop() {
  digitalWrite(pulsePin, HIGH); // This action begins the pulse
  delay(pulse_duration_ms); // This delay represents duration of pulse
  digitalWrite(pulsePin, LOW); // This action ends the pulse
  
  delay(inter_pulse_delay);
  

}
