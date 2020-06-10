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

import serial, random, datetime, numpy as np, time, datetime, json, pickle
from sys import argv
import ratterdamGUI as gui
import matplotlib.pyplot as plt

class CommChannel():

    def __init__(self, rat, TaskObj, commName):
        self.rat = rat
        self.ser = serial.Serial(commName, 115200, timeout = 1)
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
#            profile.stamp(0)
            alley = np.int8(self.TaskObj.currComm + 8) #i2c starts at 8

            incomingData = self.transceive(3,alley,5) # command 5 is a null throwaway command, same bytes/trans
#            print(incomingData)
            self.TaskObj.processUserInput_RF(incomingData)

            # if self.TaskObj.currComm == 10:
            #     print(f"user: {incomingData}")

            incomingData = self.transceive(0, alley, 9)
            
 #           print(incomingData)
            
            # if self.TaskObj.currComm == 10:
            #     print(f"track: {incomingData}")
            self.TaskObj.updateRules_DF(incomingData)
            #if self.TaskObj.currComm == 10:
            
            
            command = np.int8(self.TaskObj.nextState_DF())
            
            # if self.TaskObj.currComm == 10:
            #     print(f"command: {command}")

            throwaway = self.transceive(1, alley, command)
  #          print(throwaway)
            
            self.x +=1
            self.TaskObj.currCommIncrement()
 #           profile.look(0)
           

class Profiler():

    def __init__(self):
        self.stamps = {}

    def stamp(self, key):
        self.stamps[key] = datetime.datetime.now()

    def look(self,key):       
        print((datetime.datetime.now() - self.stamps[key]).total_seconds() * 1000)
        

    
class Foraging():

    def __init__(self, nActive, nChoices, nReplace, nAlleys = 17):
        ''' 
        nActive - number of active (rewarded) alleys at a given time. Either side accessible
        nChoices - number of choices animal makes before textures are swapped out
        nReplace - number of textures to swap out when a replacement occurs
        '''
        self.textureRecord = {i:[] for i in range(17)}
        self.stimulusExchangeRecord = []
        self.nActive = nActive
        self.numAlleys = nAlleys # change to smaller num for easier debugging. remember to change arduino code to match (for loop over the alleys ln 33 abouts)
        self.nChoices = nChoices
        self.choiceNum = 0
        self.alleySwapStatuses = [0]*17 #what is state of alley wrt alley swapping?
        self.nReplace = nReplace
        self.currComm =  0
        self.visits = [0]*17
        self.alleySwapCounts = [0]*17
        self.runStatus = 'P' # R - running, P - paused, T - terminate
        self.StimRProbCodes = {'A':7, 'B': 8, 'C': 9} # These codes correspond, in alley Ard code, to 0.75, 0.5, 0.2, respectively
        self.dataStates = [0,1]
        self.rewardNums = 0
        self.dataRecord = []
        self.textures = ['A', 'B', 'C']
        #self.texturePattern = np.random.choice(self.textures, self.numAlleys, replace=True)
        self.texturePattern = ['A' for i in range(self.numAlleys)] # for training, set all allyes to high R
        for i in range(17):
            self.updateStimulusRecord(i) 
        # rules lookup: send a unique code to alleys (via master) saying what their reward rule is
        # 0, 1, 3 - RF; 4, 7, 5, 8, 6, 9 - DF
        # Note two keys map to same val, the keys 3,4,5,6 are 'init' codes when reward state is first assigned
        # so I can tell, given a reward code, has it been sent yet or not? I see a 3 and I send a reward but know
        # that this is the first time (for this reward assignemnt) that i've done so. then set the rule to 1. Otherwise
        # I'm confused - i look at my rules and see, e.g. a 1 and look at alley and see a 0 and i dont know whether its 0 bc
        # i havent sent the rule yet or I already did and the rat got the reward (setting the val to 0). the init code disambiguates.
        self.alleyRulesLookup = {1:1, 0:0, 2:2, 3:1, 4:7, 7:7, 5:8, 8:8, 6:9, 9:9}
        self.initRules_DF = {'A':4, 'B':5, 'C': 6}
        # Initialize active alleys and texture pattern
        self.alleyRules = self.initializeRules_DF()
        self.nextTextures_DF = [self.selectNewStimuli(x) for i,x in enumerate(self.texturePattern)]
        self.printStatus_DF()


    def initializeRules(self):
        activeAlleys = np.random.choice(range(self.numAlleys), self.nActive, replace=False)
        alleyRules = np.array([0]*self.numAlleys)
        alleyRules[activeAlleys] = 1
        return alleyRules

    def initializeRules_DF(self):
        alleyRules = [self.getInitRule_DF(i) for i in range(17)]
        return alleyRules
    
    def getps(self, t):
        return [1./float(len(self.textures)-1) if i != self.textures.index(t) else 0 for i in range(len(self.textures))]

    def printTextures(self):
        for i,t in enumerate(self.texturePattern):
            print(f"Alley {i+1}: {t}")

    def selectNewStimuli(self, txt):
        cycle = {'A':'B', 'B':'C', 'C':'A'}
        anticycle = {'A':'C', 'B':'A', 'C':'B'}
        if random.random() < 0.75:
            nextStim = cycle[txt]
        else:
            nextStim = anticycle[txt]
        return nextStim


    def replaceTextures(self, i):
        '''Select which alley (i)  to replace textures and what textures to replace with'''
        self.texturePattern[i] = self.nextTextures_DF[i][0] # update current texture with the next one
        self.alleyRules[i] = self.getInitRule_DF(i)
        self.visits[i] = 0
        self.nextTextures_DF[i] = self.selectNewStimuli(self.texturePattern[i]) # select a new 'next one'
        self.updateStimulusRecord(i)
        self.printStatus_DF()

    def updateStimulusRecord(self, alley):
        self.textureRecord[alley].append(self.texturePattern[alley]) #swap has already been made, so this is 'new' txt
        self.stimulusExchangeRecord.append(((datetime.datetime.now()-timeZero).total_seconds(), alley))
        
    def nextState_RF(self):
        if self.runStatus == 'P' or self.runStatus == 'L':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]
            nextState = self.alleyRulesLookup[rule]
            if rule == 3:
                self.alleyRules[self.currComm] = 1 #update that youve now sent a reward command once
        return nextState


    def nextState_DF(self):
        if self.runStatus == 'P':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]
            nextState = self.alleyRulesLookup[rule] # this maps init codes, running codes, to running codes. so most of the time redundant, maybe simplify later
            if rule in [4, 5,  6]:
                rule = self.alleyRulesLookup[rule] # change init rule to 'normal' rule
        return nextState

    
    def  printStatus_RF(self):
        #print(self.alleyRules)
        t = self.texturePattern
        r = ['_' if i == 1 or i == 3 else '#' for i in self.alleyRules]
        print("-----------------------------------------------")
        print(f"     {t[0]}      {t[4]}      {t[6]}     \n    #{r[0]}#    #{r[4]}#    #{r[6]}#\n  {t[2]} {r[2]}##  {t[3]} {r[3]}##  {t[5]} {r[5]}#{r[7]} {t[7]}\n    #{r[1]}#    #{r[12]}#    #{r[8]}#\n     {t[1]}      {t[12]}     {t[8]}\n    ###    ###    ###\n {t[16]}  {r[16]}## {t[14]}  {r[14]}##  {t[11]} {r[11]}#{r[9]} {t[9]}\n    #{r[15]}#    #{r[13]}#    #{r[10]}#\n     {t[15]}     {t[13]}     {t[10]}")
        

    def printStatus_DF(self):
        print("------------------------------------------------------")
        for i,x in enumerate(self.texturePattern):
            print(f"Alley {i+1}: Texture {x}, Next txt: {self.nextTextures_DF[i][0]}, Count: {self.alleySwapCounts[i]}, Visits: {self.visits[i]}")

    
    def updateRules_RF(self, currentState):
        if self.runStatus == 'R':
            rule = self.alleyRules[self.currComm]

            # 0 - neither broken
            if '0' in currentState:
                pass

            # 1 - A is broken
            if '1' in currentState:
                maybeAlley = tracker.inputNewPosition_RF(self.currComm, 'A', rule)
                if maybeAlley:
                    self.alleyRules[self.currComm] = 0
                    self.alleyRules[maybeAlley] = 3
                    self.printStatus_RF()
            # 2 - B is broken
            if '2' in currentState:
                maybeAlley = tracker.inputNewPosition_RF(self.currComm, 'B', rule)
                if maybeAlley:
                    self.alleyRules[self.currComm] = 0
                    self.alleyRules[maybeAlley] = 3
                    self.printStatus_RF()

            # 3 - both are broken.
            if '3' in currentState:
                pass
        
    def updateRules_DF(self, currentState):
        maybeReleases = []
        if self.runStatus == 'R':
            # 0 - neither broken
            if '0' in currentState:
                maybeReleases = []
            if '1' in currentState:
                maybeReleases = tracker.inputNewPosition_DF(self.currComm, 'A',1)
            if '2' in currentState:
                maybeReleases  =  tracker.inputNewPosition_DF(self.currComm, 'B', 2)
            if '5' in currentState:
                maybeReleases = tracker.inputNewPosition_DF(self.currComm, 'A', 5)
                self.alleyRules[self.currComm] = 0
                self.hasVisited()
            if '6' in currentState:
                maybeReleases  = tracker.inputNewPosition_DF(self.currComm, 'B', 6)
                self.alleyRules[self.currComm] = 0
                self.hasVisited()
            if '7' in currentState:
                maybeReleases  = tracker.inputNewPosition_DF(self.currComm, 'A', 7)
                self.alleyRules[self.currComm] = 0
                self.hasVisited()
            if '8' in currentState:
                maybeReleases  = tracker.inputNewPosition_DF(self.currComm, 'B', 8)
                self.alleyRules[self.currComm] = 0
                self.hasVisited()

            # 
            if '3' in currentState:
                ''' if both broken, one must have broken first and should be what drives
                reward coding. as far as occupancy, if both broken record the invidiuals above'''
                maybeReleases = []


            if maybeReleases:
                # print(f"releases: {maybeReleases}")
                for i in maybeReleases:
                    self.alleyRules[i] = self.getInitRule_DF(i)


    def hasVisited(self):
        if self.visits[self.currComm] == 0:
            self.visits[self.currComm] = 1
            self.printStatus_DF()

                    
    def getInitRule_DF(self, i):
        """NB init doesnt mean init of session but init when rule
        is first sent, i.e. after txts swap"""
        return self.initRules_DF[self.texturePattern[i]]


    def processStimSwap(self, alley):
        if self.alleySwapStatuses[alley] == 0:
            self.alleySwapStatuses[alley] = 1
            self.alleyRules[alley] = 2
        elif self.alleySwapStatuses[alley] == 1:
            self.alleySwapStatuses[alley] = 0
            self.alleySwapCounts[alley] += 1
            self.replaceTextures(alley)
            
            
    def processUserInput(self, data):
        if '0' in data:
            self.runStatus = 'P'
            print("PAUSED")
        elif '1' in data:
            self.runStatus = 'R'
            print("RUNNING")
        elif '2' in data:
            self.runStatus = 'P'
        elif '4' in data:
            tracker.save()
        elif '3' in data:
            pass


    def padDate(self, x):
        if len(str(x)) == 1:
            return f"0{x}"
        else:
            return str(x)

        
    def save(self):
        dateTimestamp = datetime.datetime.now()
        year = str(dateTimestamp.year)[2:]
        month = self.padDate(dateTimestamp.month)
        day = self.padDate(dateTimestamp.day)
        hour = self.padDate(dateTimestamp.hour)
        minute = self.padDate(dateTimestamp.minute)
        data = {'stimuli':self.textureRecord, 'order':self.stimulusExchangeRecord}
        with open(f"R765_{year}{month}{day}_{hour}-{minute}_Stimuli.p", 'wb') as outfile:
                pickle.dump(data, outfile)
        print("SAVED STIMULI")

        
    def processUserInput_DF(self, data):
        if 'A' in data:
            alley = int(data[3:5])-1
            print(alley)
            self.processStimSwap(alley)
        elif '0' in data:
            self.runStatus = 'P'
            print("PAUSED")
        elif '1' in data:
            self.runStatus = 'R'
            print("RUNNING")
        elif '2' in data:
            self.runStatus = 'P'
        elif '4' in data:
            tracker.save()
            self.save()
        elif '3' in data:
            pass
    
    def currCommIncrement(self):
        if self.currComm == self.numAlleys - 1:
            self.currComm = 0
        else:
            self.currComm += 1
            

class PositionTracker():

    def __init__(self, task='foraging'):
        self.task = task
        self.dataRecord = {'alleys':[], 'sides':[], 'events':[], 'ts':[]}
        self.alleyVisitRec= {i:{'A':0,'B':0,'currentVisits':0, 'totalVisits':0} for i in range(17)}
        self.lockoutZone ={0: [1, 5, 6, 14, 16],
             1: [0, 4, 5, 8, 11, 13, 15],
             2: [3, 4, 12, 14, 15],
             3: [2, 5, 6, 8, 11, 13, 15, 16],
             4: [1, 2, 7, 8, 11, 12, 14],
             5: [0, 1, 3, 7, 9, 10, 13, 14],
             6: [0, 3, 8, 9, 11, 12],
             7: [4, 5, 10, 11, 12],
             8: [1, 3, 4, 6, 10, 13, 14],
             9: [5, 6, 11, 12, 13],
             10: [5, 7, 8, 12, 14, 15],
             11: [1, 3, 4, 6, 7, 9, 14, 15],
             12: [0, 2, 4, 6, 7, 9, 10, 13, 15, 16],
             13: [1, 3, 5, 8, 9, 12, 16],
             14: [0, 2, 4, 5, 8, 10, 11, 16],
             15: [1, 2, 3, 10, 11, 12],
             16: [0, 3, 12, 13, 14]}
        self.isCurrentlyLockedOut = [0]*17
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
        

    def inputNewPosition_RF(self, alley, side, rule):
        #implement better visit check alg by editing here. check
        # not only was this alley broken but was 1) its partner broken and
        # 2) 'invalidating' ring of alleys (if he looped around) not broken
        # though you need to split visit from reward. reward is still triggered
        #on 1st alley break. so can have R without 'visit'? need to rethink/work
        if rule == 1:
            newAlley = self.getNewAlley(alley, side)
            event = 'reward'
        else:
            newAlley = False
            event = 'pass'
        stamp = (datetime.datetime.now() - timeZero).total_seconds()
        for data, entryType in zip([alley, side, event, stamp], ['alleys', 'sides', 'events', 'ts']):
            self.dataRecord[entryType].append(data)


        return newAlley                

    def padDate(self, x):
        if len(str(x)) == 1:
            return f"0{x}"
        else:
            return str(x)


    def inputNewPosition_DF(self, alley, side, beamState):
        # if beamState != 3:
        #     for i,x in enumerate(self.isCurrentlyLockedOut):
        #         print(f"Alley {i+1} is {(lambda x: 'locked out' if x == 1 else 'active')(x)}")

        if beamState == 5 or beamState == 6:
            event = 'reward'
            self.isCurrentlyLockedOut[alley] = 1
        if beamState == 7 or beamState == 8:
            event = 'noReward'
            self.isCurrentlyLockedOut[alley] = 1
        else:
            event = 'pass'
        releases = []
        for i,a in enumerate(self.isCurrentlyLockedOut):
            if a == 1 and alley in self.lockoutZone[i]:
                self.isCurrentlyLockedOut[i] = 0
                releases.append(i)
        
        stamp = (datetime.datetime.now() - timeZero).total_seconds()
        for data, entryType in zip([alley, side, event, stamp], ['alleys', 'sides', 'events', 'ts']):
            self.dataRecord[entryType].append(data)

        return releases
    
    def getNewAlley(self, currentAlley, side):
        acceptables = set(range(0,17)).difference(set([currentAlley] +
                                                      self.AlleyStateGraph[currentAlley][side] +
                                                  [i for i,x in enumerate(forg.alleyRules) if x == 1]

        ))
        return np.random.choice(list(acceptables), 1)[0]

        
    def save(self):
        dateTimestamp = datetime.datetime.now()
        year = str(dateTimestamp.year)[2:]
        month = self.padDate(dateTimestamp.month)
        day = self.padDate(dateTimestamp.day)
        hour = self.padDate(dateTimestamp.hour)
        minute = self.padDate(dateTimestamp.minute)
        pathstub = "C:/Users/whock/Google Drive/KnierimLab/Ratterdam/Data/"
        rat = "R765"
        fname = f"{pathstub}{rat}/{rat}_{year}{month}{day}_{hour}-{minute}_Behavioral.json"
        with open(fname, 'w') as outfile:
            json.dump(self.dataRecord, outfile, indent=4)
        print("SAVED BEHAVIORAL")
        
 
if __name__ == '__main__':
    timeZero = datetime.datetime.now()
    profile = Profiler()
    tracker = PositionTracker()
    forg = Foraging(3,1,10,17)    
    comm = CommChannel('999',forg,'COM5')
    comm.communicate()

    
