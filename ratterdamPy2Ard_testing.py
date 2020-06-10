
import serial, random, datetime, numpy as np, time
from sys import argv


l = []
ser = serial.Serial('COM5',115200)
stillRunning = 1
ready = 1
while not ready:
    if "Ready" in str(ser.readline()):
        ready = 1
while stillRunning:
    ser.write(bytes([np.int8(3)]))
    ser.write(bytes([np.int8(10)]))
    ser.write(bytes([np.int8(3)]))
    response  = ser.read(1)
    l.append(response)
    print(response)

    
