"""
Ratterdam Stimuli Selection
William Hockeimer
V1 June 11 2018

Purpose is to pseudorandomly select textured stimuli for each alley. 
This script selects which textures should be placed at what alley
however the user inputs when this is actually done. 
For now, V1, this script does not have any information about
what alleys have actually been visited and relies on the user
to input what alleys have had their stimuli updated. Data is saved
at the end of a session.  
"""

import serial, random, datetime, numpy as np, time, datetime, json, pickle
from sys import argv
import ratterdamGUI as gui
import matplotlib.pyplot as plt


def initializeStimuli():
    good = False
    while not good:
        pattern = {i:np.random.choice(['A','B','C'],1,replace=True) for i in range(1,18)}
        l = [pattern[i] for i in range(1,18)]        
        if l.count('A') > 6 or l.count('B') > 6 or l.count('C') > 6:
            pass
        else:
            good = True
    return pattern

def printPattern(pattern, upcomingStimuli, textureRecord):
    print(f" Replacement Number: {len(updateRecord)}")
    for a,t in pattern.items():
        print(f"Alley {a}: Texture {t}, next Stim: {upcomingStimuli[a]}, Switches: {len(textureRecord[a])-1}")
        
def selectNewStimuli(textureRecord, alley):
    x = ['A','B','C']
    x.remove(textureRecord[alley][-1])
    newStim = np.random.choice(x,1)
    return newStim


def updateStimState(alley, upcomingStimuli):
    textureRecord[alley].append(upcomingStimuli[alley])
    upcomingStimuli[alley] = selectNewStimuli(textureRecord, alley)
    updateRecord.append(((datetime.datetime.now()-timeZero).total_seconds(),alley))
    return upcomingStimuli

def padDate(x):
    if len(str(x)) == 1:
        return f"0{x}"
    else:
        return str(x)


def save():
    dateTimestamp = datetime.datetime.now()
    year = str(dateTimestamp.year)[2:]
    month = padDate(dateTimestamp.month)
    day = padDate(dateTimestamp.day)
    hour = padDate(dateTimestamp.hour)
    minute = padDate(dateTimestamp.minute)
    data = {'stimuli':textureRecord, 'order':updateRecord}
    with open(f"R765_{year}{month}{day}_{hour}-{minute}_Stimuli.p", 'wb') as outfile:
            pickle.dump(data, outfile)
    print("SAVED")


pattern = initializeStimuli()
textureRecord = {i:[pattern[i]] for i in range(1,18)}
updateRecord = []
upcomingStimuli = {i:selectNewStimuli(textureRecord,i) for i in range(1,18)}
running = 1
term = 0
timeZero = datetime.datetime.now()

while running:
    if term == 1:
        save()
        running = 0
    else:
        printPattern(pattern, upcomingStimuli, textureRecord)
        userInput = input("Enter input: ")
        if userInput == 'q':
            term = 1
        elif userInput.isdigit():
            upcomingStimuli = updateStimState(int(userInput), upcomingStimuli)
            pattern = {i:textureRecord[i][-1] for i in range(1,18)}

