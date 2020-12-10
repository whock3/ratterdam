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
            incomingData = self.transceive(0, alley, 9)
            #profile.look(2)
#            print(incomingData)
            # engage any master-controlled reward contingencies. Buzzer is only one for now. 
            self.TaskObj.updateRules(incomingData)
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
        self.runStatus = 'R' # R - running, P - paused, T - terminate
        self.StimRProbCodes = {'A':2, 'B': 3, 'C': 4} # These codes correspond, in alley Ard code, to 0.75, 0.5, 0.2, respectively
        self.dataStates = [0,1]
        self.rewardNums = 0
        self.dataRecord = []
        self.textures = ['A', 'B', 'C']
        self.alleyRulesLookup = {1:1, 0:0, 3:1}

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

    def printTextures(self):
        for i,t in enumerate(self.texturePattern):
            print(f"Alley {i+1}: {t}")
    
    def replaceTextures(self, i):
        '''Select which alleys to replace textures and what textures to replace with'''
        pass

    def nextState(self):
        if self.runStatus == 'P' or self.runStatus == 'L':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]
            nextState = self.alleyRulesLookup[rule]
            if rule == 3:
                self.alleyRules[self.currComm] = 1 #update that youve now sent a reward command once
        return nextState


    def nextState_DF(self):
        if self.runStatus == 'P' or self.runStatus == 'L':
            nextState = 0
        elif self.runStatus == 'R':
            pass
        return nextState

    
    def  printStatus(self):
        #print(self.alleyRules)
        t = self.texturePattern
        r = ['_' if i == 1 or i == 3 else '#' for i in self.alleyRules]
        print("-----------------------------------------------")
        print(f"     {t[0]}      {t[4]}      {t[6]}     \n    #{r[0]}#    #{r[4]}#    #{r[6]}#\n  {t[2]} {r[2]}##  {t[3]} {r[3]}##  {t[5]} {r[5]}#{r[7]} {t[7]}\n    #{r[1]}#    #{r[12]}#    #{r[8]}#\n     {t[1]}      {t[12]}     {t[8]}\n    ###    ###    ###\n {t[16]}  {r[16]}## {t[14]}  {r[14]}##  {t[11]} {r[11]}#{r[9]} {t[9]}\n    #{r[15]}#    #{r[13]}#    #{r[10]}#\n     {t[15]}     {t[13]}     {t[10]}")
        

    
    def updateRules(self, currentState):
        if self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]

            # 0 - neither broken
            if '0' in currentState:
                pass

            # 1 - A is broken
            if '1' in currentState:
                maybeAlley = tracker.inputNewPosition(self.currComm, 'A', rule)
                if maybeAlley:
                    self.alleyRules[self.currComm] = 0
                    self.alleyRules[maybeAlley] = 3
                    #self.runStatus = 'L'
                    self.printStatus()
            # 2 - B is broken
            if '2' in currentState:
                maybeAlley = tracker.inputNewPosition(self.currComm, 'B', rule)
                if maybeAlley:
                    self.alleyRules[self.currComm] = 0
                    self.alleyRules[maybeAlley] = 3
                    #self.runStatus = 'L'
                    self.printStatus()

            # 3 - both are broken.
            if '3' in currentState:
                pass
        
        # if self.runStatus == 'L':
        #     if '1' in currentState:
        #         maybeResume = tracker.inputNewPosition(self.currComm, 'A', 'lockout')
        #     if '2' in currentState:
        #         maybeResume = tracker.inputNewPosition(self.currComm, 'B', 'lockout')
        #     if '0' in currentState:
        #         pass
        #     if maybeResume == '1':
        #         self.runStatus = 'R'


    def processUserInput(self, data):
        if 'P' in data:
            self.runStatus = 'P'
            print("PAUSED")
        elif 'R' in data:
            self.runStatus = 'R'
            print("RUNNING")
        elif '2' in data:
            self.runStatus = 'P'
        elif 'S' in data:
            tracker.save()
        elif 'M' in data:
            tracker.manual_sync()
        else:
            pass
    
    def processUserInput_DF(self, data):
        for i in range(1,18):
            if str(i) in data:
                self.replaceTextures(i)

    def currCommIncrement(self):
        if self.currComm == self.numAlleys - 1:
            self.currComm = 0
        else:
            self.currComm += 1
            

class PositionTracker():

    def __init__(self,rat, task='foraging'):
        self.task = task
        self.rat = rat
        self.dataRecord = {'alleys':[], 'sides':[], 'events':[], 'ts':[]}
        self.timeZero = datetime.datetime.now()
        self.alleyVisitRec= {i:{'A':0,'B':0,'currentVisits':0, 'totalVisits':0} for i in range(1,18)}
        self.lockoutZone = {1:[2, 6, 7, 1, 15, 17],2:[1, 5, 6, 9, 12, 14, 16],3:[4, 5, 13, 15, 16],4:[3, 6, 7, 9, 12, 14, 16, 17],5:[2, 3, 8, 9, 12, 13, 15], 6:[1,2, 4, 8, 10, 11, 14, 15],7:[1, 4, 9, 10, 12, 13],8:[5, 6, 11, 12, 13],9:[2,4,5,7,11,14,15],10:[6,7,12,13,14],11:[6,8,9,13,15,16],12:[2,4,5,7,8,10,15,16],13:[1,3,5,7,8,10,11,14,16,17],14:[2,4,6,9,10,13,17], 15:[1,3,5,6,9,11,12,17],16:[2,3,4,11,12,13], 17:[1,4,13,14,15]}
        self.currentLockoutEscapes = []
        self.AlleyStateGraph = {0:{'A':[3, 4], 'B':[2]},
                        1:{'A':[3, 12, 14], 'B':[2, 16]},
                        2:{'A':[1, 16], 'B':[0]},
                        3:{'A':[1, 12, 14], 'B':[0, 4]},
                        4:{'A':[5, 6], 'B':[0, 3]},
                        5:{'A':[8, 11, 12], 'B':[4, 6]},
                        6:{'A':[7], 'B':[4, 5]},
                        7:{'A':[8, 9], 'B':[6]},
                        8:{'A':[7, 9], 'B':[5, 11, 12]},
                        9:{'A':[10], 'B':[7, 8]},
                        10:{'A':[9], 'B':[11, 13]},
                        11:{'A':[13, 10], 'B':[5, 8, 12]},
                        12:{'A':[5, 8, 11], 'B':[3, 14, 1]},
                        13:{'A':[10, 11], 'B':[14, 15]},
                        14:{'A':[13, 15], 'B':[1, 3, 12]},
                        15:{'A':[13, 14], 'B':[16]},
                        16:{'A':[15], 'B':[1, 2]}}
        self.init_savepath()
        
        
        
    def init_savepath(self):  
        dateTimestamp = datetime.datetime.now()
        year = str(dateTimestamp.year)[2:]
        month = self.padDate(dateTimestamp.month)
        day = self.padDate(dateTimestamp.day)
        hour = self.padDate(dateTimestamp.hour)
        minute = self.padDate(dateTimestamp.minute)
        pathstub = "C:\\Users\\whock\\Google Drive\\KnierimLab\\Ratterdam\\Data\\"
        self.fname = f"{pathstub}{self.rat}_{year}{month}{day}_{hour}-{minute}_Behavioral.json"
        
    def inputNewPosition(self, alley, side, rule):
        if rule == 1:
            newAlley = self.getNewAlley(alley, side)
            event = 'reward'
            self.currentLockoutEscapes = self.lockoutZone[alley]

        # elif rule == 'lockout':
        #     if alley in self.currentLockoutEscapes:
        #         newAlley = '1'
        #     else:
        #         newAlley = '0'
        else:
            newAlley = False
            event = 'pass'
        stamp = (datetime.datetime.now() - self.timeZero).total_seconds()
        for data, entryType in zip([alley, side, event, stamp], ['alleys', 'sides', 'events', 'ts']):
            self.dataRecord[entryType].append(data)


        return newAlley         


    def manual_sync(self):
        stamp = (datetime.datetime.now() - self.timeZero).total_seconds()
        for data, entryType in zip(['xxx', 'xxx', 'sync', stamp], ['alleys', 'sides', 'events', 'ts']):
            self.dataRecord[entryType].append(data)
        
    
        
    
    def getNewAlley(self, currentAlley, side):
        la = [i-1 for i in [1,2,3,5,6,7,8,9,10,11,13,14,15,17]]
        acceptables = set(la).difference(set([currentAlley] +
                                                      self.AlleyStateGraph[currentAlley][side] +
                                                  [i for i,x in enumerate(forg.alleyRules) if x == 1]

        ))
        return np.random.choice(list(acceptables), 1)[0]


    def padDate(self,x):
        if len(str(x)) == 1:
            return f"0{x}"
        else:
            return str(x)
        
    def save(self):
        with open(self.fname, 'w') as outfile:
            json.dump(self.dataRecord, outfile, indent=4)
        print("SAVED")
        
 
if __name__ == '__main__':
    rat = 'R781'
    profile = Profiler()
    tracker = PositionTracker(rat)
    # nActive, nChoices, nReplace, nAlleys
    forg = Foraging(2,1,10,17)    
    comm = CommChannel(rat,forg,'COM5')
    comm.communicate()

    
