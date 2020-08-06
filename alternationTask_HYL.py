# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 18:49:41 2020

@author: whockei1
"""

import numpy as np, time, datetime, serial

iti_len = 60 # time in start box
num_r = 6
iri = 5 # gap between rewards
nTrials = 100
tepsilon = 0.1
alley = 20
rat = '999'
ser = serial.Serial('COM5', 115200, timeout = 1)

def generateRewardTimes():
    sampleFound = False
    while not sampleFound:
        sample = np.random.randint(0,iti_len,num_r)
        gaps = np.diff(sorted(sample))
        if min(gaps) >= iri:
            sampleFound = True
    return sample

def transceive(dirVal, alleyVal, commandVal):
    direction, alley, command  = map(np.int8, [dirVal, alleyVal, commandVal])
    for info in [direction, alley, command]:
        ser.write(bytes([info]))
    incomingData = ser.readline()
    return str(incomingData)

def padDate(x):
    if len(str(x)) == 1:
        return f"0{x}"
    else:
        return str(x)

def initTimestamp():
    """init timestamp so you can reopen a file on autosave and not make
    nee file per lap"""
    dateTimestamp = datetime.datetime.now()
    year = str(dateTimestamp.year)[2:]
    month = padDate(dateTimestamp.month)
    day = padDate(dateTimestamp.day)
    hour = padDate(dateTimestamp.hour)
    minute = padDate(dateTimestamp.minute)
    
    fname = f"{rat}_{year}{month}{day}_{hour}-{minute}_BehavioralData.csv"
    return fname

    
def communicate(self):
    allSamples = []
    deliveredTrials = []
    trialStartTimes = []
    for n in range(nTrials):
        allSamples.append(generateRewardTimes())      
    

    stillRunning = 1
    ready = 0
    
    while not ready:
        if "Ready" in str(self.ser.readline()):
            ready = 1
            
    while stillRunning:
        
        for i, trial in enumerate(allSamples):
            trialStart = False
            
            while not trialStart:
                incoming = self.transceive(0, alley, 9)
                if incoming == '1':
                    trialStart = True
                
            trialStartTime = time.time()
            trialStartTimes.append(trialStartTime)
            deliveredTrials.append([i,trialStartTime,trial])
            
            for rewardTime in trial:
                if rewardTime-(time.time()- trialStartTime) <= tepsilon:
                    print(f"Reward Delivered at {rewardTime}s")
                    _ = self.transceive(1, alley, 1)
        
        
        
    

