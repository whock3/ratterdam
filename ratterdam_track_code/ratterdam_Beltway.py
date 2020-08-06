"""
18-03-07
William Hockeimer
Knierim Lab
Ratterdam Project - "Beltway" subtask.

In this task only outer perimeter of ratterdam is available to the rat, consisting
of 10 alleys (one being start box). Each alley has a active lickport with some prob of reward and a texture
that is swapped between laps. After a lap around the perimeter a confinment start box
is placed around rat while textures are swapped.

Comm logic is that each alley is taken in turn for a read-read-write process.
First read is for user input, second for alley status, and after registering any 
changes to status the new state is send via write.

Changes to basic ratterdam design:
    - Only poll alleys in perimeter. Maybe consider in order of direction of travel
    - lockout relief is lap based. ie. until new lap starts.


old info:
Communication script - python to arduino mega in Ratterdam V2
Control logic: Each trial is a string, 17 characters long (for 17 alleys)
Each character byte is an instruction for the alley number corresponding to that idx
i.e 5th entry for alley 5. Code is: A (A+,B-), B(B+,A-), E(A-,B-), N(A0,B0), O (A+,B+)
This is sent to master arduino which in turn sends each char to the child ards.
String passed back to this python script is a 17 char long script with statuses at time
of trial outcome. Statuses - N - waiting, A - A is broken, A - B is broken
"""

import serial, datetime, numpy as np, time, json, pickle, itertools
from collections import OrderedDict
from sys import argv
import csv
import matplotlib.pyplot as plt

class CommChannel():

    def __init__(self, rat, TaskObj, commName):
        self.rat = rat
        self.ser = serial.Serial(commName, 115200, timeout = 1)
        self.TaskObj = TaskObj

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
            self.TaskObj.processUserInput(incomingData)


            incomingData = self.transceive(0, alley, 9)
            self.TaskObj.updateRules(incomingData)
            
            
            command = np.int8(self.TaskObj.nextState())
            _ = self.transceive(1, alley, command)
            
            self.TaskObj.currCommIncrement()
           

class RatterdamBeltway():
    
    def __init__(self, rat, direction):
        self.timeZero = datetime.datetime.now()
        #self.savepath = "C:\\Users\\whockei1\\Google Drive\\KnierimLab\\Ratterdam\\Beltway\\BehavioralData\\"
        self.savepath = 'C:\\Users\\admin\\Desktop\\ratterdam_code\\'
        self.rat = rat
        self.direction = direction
        self.currentLap = 0
        self.lapTS = [] # ts in python time for beginning of a lap (actually ITI preceeding lap)
        
        self.beltwayAlleyList = [i-1 for i in [16, 17, 3, 1, 5, 7, 8, 10, 11]]
        
        #initialize alleys
        self.alleyManager = AlleyManager()
        self.alleys = self.alleyManager.generateAlleyData()
        
        
        self.currCommCycler = itertools.cycle(self.beltwayAlleyList)
        self.currComm = next(self.currCommCycler)
        
        self.alleyRulesLookup = {1:1, 0:0, 2:2, 3:1, 4:7, 7:7, 5:8, 8:8, 6:9, 9:9}
        self.rewardProbs = {1:0.15, 2:0.7, 3:0.15} # num rewards: prob of that num/lap
        self.rewardWindow = 4 # number of laps to look back when enforcing even reward sampling
        self.rewardPaucityThresh = 0.5
        self.runStatus = 'P'
        self.initTimestamp()
        
        
        
    def processUserInput(self, userInput):
        userInput = str(userInput)
        if 'P' in userInput:
            self.runStatus = 'P'
            print("PAUSED")
        elif 'R' in userInput:
            self.runStatus = 'R'
            print("RUNNING")
        elif 'X' in userInput:
            self.runStatus = 'P'
            self.nextLap()
            self.save()
        elif 'S' in userInput:
            self.save()
            print("SAVED")
        elif 'T' in userInput:
            self.printStatus()
        elif 'Q' in userInput:
            self.writeDataForGUI(ratquit=True)
            self.save()
            
    def updateRules(self, currentState):
        
        
        if self.runStatus == 'R':
        
            
            rule = self.alleys[self.currComm].rule
            
            if '0' in currentState:
                pass
            
            
            elif ('1' in currentState or '2' in currentState) and rule == 1:
                #note 1,,3 refers to which IR is broken (3==both, not possible but kept in case of bugs). But for beltway there's only one way he can travel at a time so no need to keep track.
                self.alleys[self.currComm].rule = 0
            
            
    def nextState(self):
        if self.runStatus == 'P':
            nextState = 0
        elif self.runStatus == 'R':
            rule = self.alleys[self.currComm].rule
            nextState = self.alleyRulesLookup[rule]
            
            #init rule is 1, once youve sent change to 1. to i can tell between a rule i just set and havent sent vs have sent
            if rule == 3:
                self.alleys[self.currComm].rule = 1 
        return nextState
            
    
    def nextLap(self):
        """
        Select which alleys to reward and update texture pattern
        Print status. 
        """
        self.currentLap += 1
        
        self.lapTS.append((datetime.datetime.now() - self.timeZero).total_seconds())
        
        self.updateRewards()
        
        self.printStatus()
        
        self.writeDataForGUI(ratquit='False')
        
    
    def writeDataForGUI(self,ratquit = False):
        """ Write a simple .txt with txts and r locs
        File consists of 1 line. 10 csvs. 9 pairs of one letter
        txt code and 0,1 for r or not. 10th entry is 'X' for exit task
        or 'Z' for pleaZe continue.
        """
        line = ''
        for alley in self.beltwayAlleyList:
            nextTxt = self.alleys[alley].textureHistory[self.currentLap-1]
            reward = self.alleys[alley].rule
            line += f"{nextTxt}{reward},"
        if not ratquit:
            line += 'Z'
        else:
            line += 'X'
        
        with open(self.savepath+'beltway_gui_data.txt', "w", newline='') as file:
            file.write(line)
            

        
    def textureRewardCorrAtAlley(self, alley):
        failToggle = False
        pairs = zip(self.alleys[alley].rewardHistory, self.alleys[alley].textureHistory[:self.currentLap]) # dont need currentlap-1 bc slices are right exclusive,
        txtCounts = {a:0 for a in ['A','B','C']}
        for item in pairs:
            txtCounts[item[1]] += item[0]
            
        for txt in ['A','B','C']:
            if sum(self.alleys[alley].rewardHistory) > 0:
                if txtCounts[txt]/sum(self.alleys[alley].rewardHistory) > 0.5:
                    failToggle = True
                
        return failToggle
        
        

    def getLessFrequentlyRewardedAlleys(self, txtByRBalance=True):
        lessRewardedAlleys = []
        for alley in self.beltwayAlleyList:
            rewardHist = self.alleys[alley].rewardHistory[-1*self.rewardWindow:]
            
            if txtByRBalance:
                if (sum(rewardHist)/len(rewardHist) < self.rewardPaucityThresh) and not self.textureRewardCorrAtAlley(alley):
                    lessRewardedAlleys.append(alley)
                    
            else:
                if sum(rewardHist)/len(rewardHist) < self.rewardPaucityThresh:
                    lessRewardedAlleys.append(alley)
                
                
        return lessRewardedAlleys
            
    
    def updateRewards(self):
        
        numRewards = np.random.choice(list(self.rewardProbs.keys()), size=1, p=list(self.rewardProbs.values()))
                
            
        if self.currentLap < self.rewardWindow:
            rewardedAlleys = np.random.choice(self.beltwayAlleyList, numRewards)
        else:
            
            # Get alleys that have been less often rewarded in lookback window
            # if the balanceTxtR toggle is true, also check that the txt on the less
            # often chosen alleys hasnt been associated with r too much
            # if this isnt possible to achieve just get a less commonly rewarded alley
            
            possibleAlleys = self.getLessFrequentlyRewardedAlleys(txtByRBalance=False)
            if len(possibleAlleys) > 0:
                rewardedAlleys = np.random.choice(possibleAlleys, numRewards)
            else:
                possibleAlleys = self.getLessFrequentlyRewardedAlleys(txtByRBalance=False)
                rewardedAlleys = np.random.choice(possibleAlleys, numRewards)
                
            
            
        for alley in self.beltwayAlleyList:
            if alley in rewardedAlleys:
                self.alleys[alley].rule = 3
                self.alleys[alley].rewardHistory.append(1)
            else:
                self.alleys[alley].rule = 0
                self.alleys[alley].rewardHistory.append(0)
                
#        for alley in rewardedAlleys:
#            self.alleys[alley].rule = 3 # note the init rule is 3, not 1
#            self.alleys[alley].rewardHistory.append(1)
#            
#        #following set operation gets non rewarded alleys by taking set differences
#        for alley in set(self.beltwayAlleyList).difference(set(rewardedAlleys)):
#            self.alleys[alley].rule = 0
#            self.alleys[alley].rewardHistory.append(0)
                

            
    
    def printStatus(self):
        
        print(f"============== Lap {self.currentLap} ===============")
        for i,alley in enumerate(self.beltwayAlleyList):
            if i%3 == 0:
                print("------------------------------------------------------")
            nextTxt = self.alleys[alley].textureHistory[self.currentLap-1]
            reward = self.alleys[alley].rule
            print(f"Alley {alley+1}, Texture {nextTxt}, Reward? {reward}")
        print("====================================================")
            
        
        
    def currCommIncrement(self):
        
        self.currComm = next(self.currCommCycler)
        
    def padDate(self, x):
        if len(str(x)) == 1:
            return f"0{x}"
        else:
            return str(x)
        
    def initTimestamp(self):
        """init timestamp so you can reopen a file on autosave and not make
        nee file per lap"""
        dateTimestamp = datetime.datetime.now()
        year = str(dateTimestamp.year)[2:]
        month = self.padDate(dateTimestamp.month)
        day = self.padDate(dateTimestamp.day)
        hour = self.padDate(dateTimestamp.hour)
        minute = self.padDate(dateTimestamp.minute)
        
        self.fname = f"{self.rat}_{year}{month}{day}_{hour}-{minute}_BehavioralData.csv"
        
        
    def save(self):
        with open(self.savepath+self.fname, "w", newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["direction", self.direction])
            writer.writerow(["rat", self.rat])
            writer.writerow(["time zero", self.timeZero])
            writer.writerow(["lap starts"] + self.lapTS)
            for alley in self.beltwayAlleyList:
                rewards = list(self.alleys[alley].rewardHistory)
                textures = list(self.alleys[alley].textureHistory[:self.currentLap])
                writer.writerow([f"alley {alley} rewards"] + rewards)
                writer.writerow([f"alley {alley} textures"] + textures)
                



class AlleyManager():

    def __init__(self):
        self.alleyIDs = [i-1 for i in [16, 17, 3, 1, 5, 7, 8, 10, 11]]
        self.alleys = OrderedDict({a:[] for a in self.alleyIDs})
        for i in range(40): # arg to range * 3 = num trials. so 40 trial blocks = 120 trials
            self.makeTrialBlock()
        
        
    def isValidTrial(self, x,s=4):
        valid=True
        for i in range(x.shape[0]):
            chunk = x[0+i:s+i]
            for t in ['A','B','C']:
                count = np.where(chunk==t)[0].shape[0]
                if count == s:
                    valid = False
        return valid
    
    def addTrials(self, candidate):
        for i,a in enumerate(self.alleyIDs):
            self.alleys[a].append(candidate[i])
    
    def makeTrialBlock(self):
        isValid = False
        while not isValid:
            candidate = np.random.choice(['A','B','C'], 9)
            isValid = self.isValidTrial(candidate)
            
        self.addTrials(candidate)
        
        alleyTxtsRec = {a:['A','B','C'] for a in self.alleyIDs}
        
        for a in self.alleyIDs:
            alleyTxtsRec[a].remove(self.alleys[a][-1])
        
        for i in range(2):
            c = 0
            isValid = False
            if not isValid and c < 100:
                candidate = []
                for a in self.alleyIDs:
                    candidate.append(np.random.choice(alleyTxtsRec[a]))
                isValid = self.isValidTrial(np.asarray(candidate))
                c+=1
            else:
                candidate = []
                for a in self.alleyIDs:
                    candidate.append(np.random.choice(alleyTxtsRec[a]))
            
            self.addTrials(candidate)
            
            for a in self.alleyIDs:
                alleyTxtsRec[a].remove(self.alleys[a][-1])
                
        
    def generateAlleyData(self):
        alleys = {a:Alley(a) for a in self.alleyIDs}
        for alley in self.alleyIDs:
            alleys[alley].textureHistory = self.alleys[alley]
            
        return alleys
            
        
             
                
class Alley():
    
    def __init__(self, alleyNum):
        self.alleyNum = alleyNum
        self.isRewarded = False
        self.rule = 0 #rewarded or not this lap 
        self.rewardHistory = [] #will be of len(actual laps ran)
        self.textureHistory = []
        

    
if __name__ == '__main__':
    beltway = RatterdamBeltway('999', 'CW')  
    comm = CommChannel('999', beltway, 'COM9')
    comm.communicate()