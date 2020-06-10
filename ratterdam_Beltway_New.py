import serial, random, datetime, numpy as np, time, datetime, json, pickle,itertools
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
            alley = np.int8(self.TaskObj.currComm + 8) #i2c starts at 8

            incomingData = self.transceive(3,alley,5) # command 5 is a null throwaway command, same bytes/trans
            self.TaskObj.processUserInput_DF(incomingData)


            incomingData = self.transceive(0, alley, 9)
            self.TaskObj.updateRules_DF(incomingData)
            
            
            command = np.int8(self.TaskObj.nextState_DF())
            _ = self.transceive(1, alley, command)
            
            self.x +=1
            self.TaskObj.currCommIncrement()


class Alley():
    
    def __init__(self,alleyNum):
        self.alleyNum = alleyNum
        self.isLockedOut = False
        self.lapsCompleted = 0
        self.state = 0 # RF state is 0,1,3. 0-no R, 1-R and this has been sent, 3 - R but this info hasnt been sent
        self.rewardHistory = [] # list where each lap adds a 0 (was rewarded) or 1 (was rewarded)
        self.textures = ['A', 'B', 'C']
        
        self.initRules = {'A':3, 'B':3, 'C': 3}
        
    def createTrialHistory(self):
        """
        Initialize the session's stimulus list for each trial.
        Trials are grouped into 3s where each is a random w/o replacement draw
        of the three stims. E.g a list could be ABCCBABACACB...
        """
        thist = np.empty((0), dtype='<U1')
        
        for i in range(500): # make way more than you need
            cycle = np.random.choice(['A','B','C'],3,replace=False)
            thist = np.hstack((thist, cycle))
            
        self.trialHistory = thist
        self.updateStimulus()
        
        
    def updateStimulus(self):
        self.texture = self.trialHistory[self.lapsCompleted]
        self.lapsCompleted += 1
    
    
class Foraging_Beltway():
    
    def __init__(self, alleyList=[11,14,16,17,3,1,5,7,8,10]):
        self.numAlleys = len(alleyList)
        self.alleyList = alleyList
        self.lapsCompleted = 0
        self.rewardProb = 3 # expressed as actual number of rewarded alleys
        self.runStatus = 'P'
        self.commCycler = itertools.cycle(self.alleyList)
        self.currComm = next(self.commCycler)
        self.alleyRulesLookup = {1:1, 0:0, 2:2, 3:1, 4:7, 7:7, 5:8, 8:8, 6:9, 9:9}
        
        self.initRewardedAlleys()

        
    def createAlleys(self):
        """
        Create set of alley objects
        """
        
        self.alleys = {i: Alley(i) for i in self.alleyList}
        
        for alley in self.alleyList:
            self.alleys[alley].createTrialHistory()
            
            
    def initRewardedAlleys(self):
        """
        Initialize which alleys are rewarded for lap 0
        """
        rewardedAlleys = np.random.choice(self.alleyList, self.rewwardProb, replace=False)
        for alley in rewardedAlleys:
            self.alleys[alley].state = 3
            self.alleys[alley].rewardHistory.append(1)
            
        nonrewardedAlleys = [alley for alley in self.alleyList if alley not in rewardedAlleys]
        for alley in nonrewardedAlleys:
            self.alleys[alley].state = 0
            self.alleys[alley].rewardHistory.append(0)
            
    def selectRewardedAlleys(self):
        """
        On a lap a number of alleys given by self.rewardProb will be selected
        to be rewarded. Select them here and set their init reward state (3).
        Balance the rewarding s.t. alleys are not rewarded on subsequent laps
        """
        
        rewardHist = [(alley.rewardHistory[-1], alley.alleyNum) for alley in self.alleys]
        rewardedAlleys = np.random.choice([i[1] for i in rewardHist if i[0] == 0],self.rewardProb, replace=False)
        for alley in rewardedAlleys:
            self.alleys[alley].state = 3
            self.alleys[alley].rewardHistory.append(1)
            
        nonrewardedAlleys = [alley for alley in self.alleyList if alley not in rewardedAlleys]
        for alley in nonrewardedAlleys:
            self.alleys[alley].state = 0
            self.alleys[alley].rewardHistory.append(0)
        
            
            
    def processUserInput(self, data):
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
            #tracker.save() TODO - save fx 
            self.save()
        elif '3' in data:
            pass
        
        
    def updateRules(self, currentState):
        if self.runStatus == 'R':
            
            # 0 - neither broken
            if '0' in currentState:
                pass

            # 1 - A is broken, 2-B is broken
            if '1' in currentState or '2' in currentState:
                self.alleys[self.currComm].state = 0

            # 3 - both are broken.
            if '3' in currentState:
                pass
            
                        
    def nextState(self):
        if self.runStatus == 'P' or self.runStatus == 'L':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleys[self.currComm].state
            nextState = self.alleyRulesLookup[rule]
            if rule == 3:
                self.alleys[self.currComm].state = 1 #update that youve now sent a reward command once
        return nextState
    
    def instructAlleyStimUpdate(self):
        for alley in self.alleyList:
            self.alleys[alley].updateStimulus
    
    def currCommIncrement(self):
        nextAlley =  next(self.commCycler)
        if nextAlley == self.alleyList[0]:
            self.lapsCompleted += 1
            self.selectRewardedAlleys()
            self.instructAlleyStimUpdate()
            
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        