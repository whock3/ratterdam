"""
18-03-07
William Hockeimer
Knierim Lab
Ratterdam Project

Communication script - python to arduino mega in Ratterdam V1
Control logic: Each trial is a string, 17 characters long (for 17 alleys)
Each character byte is an instruction for the alley number corresponding to that idx
i.e 5th entry for alley 5. Code is: A (A+,B-), B(B+,A-), E(A-,B-), N(A0,B0), O (A+,B+)
This is sent to master arduino which in turn sends each char to the child ards.
String passed back to this python script is a 17 char long script with statuses at time
of trial outcome. Statuses - N - waiting, A - A is broken, A - B is broken
"""

import serial, random, datetime, numpy as np, time, datetime, json
from sys import argv
import ratterdamGUI as gui
import matplotlib.pyplot as plt

class CommChannel():

    def __init__(self, rat, TaskObj, commName):
        self.rat = rat
        self.ser = serial.Serial(commName, 115200)
        self.TaskObj = TaskObj
        self.x = 0

    def transceive(self, dirVal, alleyVal, commandVal):
        direction, alley, command  = map(np.int8, [dirVal, alleyVal, commandVal])
        for info in [direction, alley, command]:
            self.ser.write(bytes([info]))
        incomingData = self.ser.readline()
        return str(incomingData)

        
    def communicate(self):
        stillRunning = 1
        ready = 0
        while not ready:
            if "Ready" in str(self.ser.readline()):
                ready = 1
        while stillRunning:
            profile.stamp(0)
            # print(self.TaskObj.runStatus)
            # print("-------------------------")
            # print(self.x)
            # print(self.TaskObj.alleyRules)
            # print(self.TaskObj.currComm)

            alley = np.int8(self.TaskObj.currComm + 8) #i2c starts at 8
            # #check to see if there are any user-generated diagnostic commands
            #profile.stamp(1)
            incomingData = self.transceive(3,alley,5) # command 5 is a null throwaway command, same bytes/trans
            #profile.look(1)
            self.TaskObj.processUserInput(incomingData)

            # Readfrom alley. Use a toggled while loop to check for a diagnostic state
            #profile.stamp(2)
            incomingData = self.transceive(0, alley, 5)
            #profile.look(2)
#            print(incomingData)
            # engage any master-controlled reward contingencies. Buzzer is only one for now. 
            reward = self.TaskObj.updateRules(incomingData)
            # if reward == 1:
            #  _ =  self.transceive(2, alley, 3)

             #Send to alley
            command = np.int8(self.TaskObj.nextState())
 #           print(command)
            #profile.stamp(3)
            _ = self.transceive(1, alley, command)
            #profile.look(3)
            
            self.x +=1
            self.TaskObj.currCommIncrement()
           # profile.look(0)

class Profiler():

    def __init__(self):
        self.stamps = {}

    def stamp(self, key):
        self.stamps[key] = datetime.datetime.now()

    def look(self,key):       
        #print((datetime.datetime.now() - self.stamps[key]).total_seconds() * 1000)
        pass

    
class Foraging():

    def __init__(self, nActive, nChoices, nReplace, nAlleys = 17):
        ''' 
        nActive - number of active (rewarded) alleys at a given time. Either side accessible
        nChoices - number of choices animal makes before textures are swapped out
        nReplace - number of textures to swap out when a replacement occurs
        '''
        self.nActive = nActive
        self.numAlleys = nAlleys # change to smaller num for easier debugging. remember to change arduino code to match (for loop over the alleys ln 33 abouts)
        self.nChoices = nChoices
        self.choiceNum = 0 
        self.nReplace = nReplace
        self.currComm =  0
        self.runStatus = 'P' # R - running, P - paused, T - terminate
        self.textures = ['A','B','C','D']
        self.dataStates = [0,1]
        self.rewardNums = 0
        self.alleyRulesLookup = {0:0, 1:1, 3:1}
        self.dataRecord = []

        # Initialize active alleys and texture pattern
        self.alleyRules = self.initializeRules()
        self.texturePattern = np.random.choice(self.textures, self.numAlleys, replace=True)
        self.printStatus()


    def initializeRules(self):
        activeAlleys = np.random.choice(range(self.numAlleys), self.nActive, replace=False)
        alleyRules = np.array([0]*self.numAlleys)
        alleyRules[activeAlleys] = 1
        return alleyRules
                               

    def getps(self, t):
        return [1./float(len(self.textures)-1) if i != self.textures.index(t) else 0 for i in range(len(self.textures))]
    
    def replaceTextures(self):
        '''Select which alleys to replace textures and what textures to replace with'''
        replacementAlleys = np.random.choice(range(self.numAlleys), self.nReplace, replace=False)
        for r in replacementAlleys:
            self.texturePattern[r] = np.random.choice(self.textures, 1, p=self.getps(self.texturePattern[r]))[0]
        print(self.texturePattern)
        self.dataRecord.append(str(self.texturePattern))
        
    def nextState(self):
        if self.runStatus == 'P':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]
            nextState = self.alleyRulesLookup[rule]
            if rule == 3:
                self.alleyRules[self.currComm] = 1 #update that youve now sent a reward command once
        return nextState

    
    def updateChoices(self, currentState):
        self.choiceNum += 1
        if self.choiceNum >= self.nChoices:
            self.runStatus = 'P'
            self.choiceNum = 0
            self.replaceTextures()

    def  printStatus(self):
        #print(self.alleyRules)
        t = self.texturePattern
        r = ['_' if i == 1 or i == 3 else '#' for i in self.alleyRules]
        print("-----------------------------------------------")
        print(f"     {t[0]}      {t[4]}      {t[6]}     \n    #{r[0]}#    #{r[4]}#    #{r[6]}#\n  {t[2]} {r[2]}##  {t[3]} {r[3]}##  {t[5]} {r[5]}#{r[7]} {t[7]}\n    #{r[1]}#    #{r[12]}#    #{r[8]}#\n     {t[1]}      {t[12]}     {t[8]}\n    ###    ###    ###\n {t[16]}  {r[16]}## {t[14]}  {r[14]}##  {t[11]} {r[11]}#{r[9]} {t[9]}\n    #{r[15]}#    #{r[13]}#    #{r[10]}#\n     {t[15]}     {t[13]}     {t[10]}")
        
    
    def updateRules(self, currentState):
        rule = self.alleyRules[self.currComm]
        reward = 0
         #this is only for O-rule. would need to change for other policy 
        if '7' in currentState or '0' in currentState:
            self.dataRecord.append(f"{self.currComm}-{currentState}")
        if rule == 1 and '1' in currentState:
            pass
        elif rule == 1 and '0' in currentState:
            self.alleyRules[self.currComm] = 0
            self.updateDataRecord(0)
            self.rewardNums += 1
            newAlley = self.getNewAlley()
            self.alleyRules[newAlley] = 3
            self.dataRecord.append(str(self.alleyRules))
            self.printStatus()
            reward = 1
        elif rule == 0 and '7' in currentState:
            self.updateDataRecord(7)
    
        return reward

    def updateDataRecord(self,state):
        pass #self.dataRecord.append(f"{self.currComm}-{state}")

    def processUserInput(self, data):
        if '0' in data:
            self.runStatus = 'P'
            print("PAUSED")
        elif '1' in data:
            self.runStatus = 'R'
            print("RUNNING")
        elif '2' in data:
            self.runStatus = 'P'
            self.replaceTextures()
        elif '4' in data:
            self.save()
        elif '3' in data:
            pass
    
    def currCommIncrement(self):
        if self.currComm == self.numAlleys - 1:
            self.currComm = 0
        else:
            self.currComm += 1
            
    def getNewAlley(self):
        inactiveAlleys = [i for i,x in enumerate(self.alleyRules) if x == 0 and i != self.currComm]
        newAlley = np.random.choice(inactiveAlleys, 1)
        return newAlley

    def save(self):
        with open('R765_180426.txt', 'w') as outfile:
            json.dump(self.dataRecord, outfile)
        print("SAVED")
        
 
if __name__ == '__main__':
    profile = Profiler()
    forg = Foraging(12,1,10,17)    
    comm = CommChannel('999',forg,'COM5')
    comm.communicate()



    
