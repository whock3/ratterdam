"""
Ratterdam Realtime Behavior Visualization
WH Late March 2018

Creates dynamic plot of rat trajectory thru maze
Annotated with info about session
Code logic: update-clear-plot. Least backend-y

"""
import numpy as np, matplotlib.pyplot as plt, random, os, matplotlib.patches as patches
import time

class RattGUI():

    def __init__(self):
        self.alleyCoords = {1:[1.5, 4.5], 2:[1.5, 2.5], 3:[0.5, 3.5], 4:[2.5, 3.5], 5:[3.5, 4.5],
                            6:[4.5, 3.5,], 7:[5.5, 4.5], 8:[6.5, 3.5], 9:[5.5, 2.5], 10:[6.5, 1.5],
                            11:[1.5, 4.5], 12:[4.5, 1.5], 13:[3.5, 2.5], 14:[3.5, 0.5], 15:[2.5, 1.5],
                            16:[1.5, 0.5], 17:[0.5, 1.5]}
        self.blockCoords = {1:[1,3], 2:[3,3], 3:[5,3], 4:[5,1], 5:[3,1], 6:[1,1]} # origin is BL for mpl patch
        self.blockRects = {}
        self.boundsCoords = [0,0,7,5] # origin x, origin y, width, height
        self.alleyStateColors = ['k']*17
        self.annLocs = {1:[0.15, 0.39, 0],
                        2:[0.15, 0.31, 0],
                        3:[0.11, 0.35, 90],
                        4:[0.31, 0.35, 90],
                        5:[0.35, 0.39, 0],
                        6:[0.39, 0.35, 90],
                        7:[0.55, 0.39, 0],
                        8:[0.59, 0.35, 90],
                        9:[0.55, 0.31, 0],
                        10:[0.59, 0.15, 90],
                        11:[0.55, 0.11, 0],
                        12:[0.51, 0.15, 90],
                        13:[0.35, 0.19, 0],
                        14:[0.35, 0.11, 0], 
                        15:[0.19, 0.15, 90],
                        16:[0.15, 0.11, 0],
                        17:[0.11, 0.15, 90]}
                        
        self.alleyHistory = {i:{'rewards':0, 'passes':0} for i in range(17)}
        self.texturePattern = []
        self.fig = plt.figure(figsize=(28, 22))
        self.axes = self.fig.add_subplot(111,aspect='equal')
        self.axes.set_ylim([0,0.41])
        self.axes.set_xlim([0,0.63])
        self.x = []
        self.y = []
        self.offset =0 
        self.markers = [] # to denote event type at traj data pts 
        
    def drawRatterdam(self):
        for blockNum, origin in self.blockCoords.items():
            r = patches.Rectangle((origin[0]/10, origin[1]/10), 1/10, 1/10, edgecolor='b')
            self.axes.add_patch(r)
        plt.axis('off')
        plt.suptitle("Ratterdam GUI", fontsize=30)
        self.axes.set_ylim([0,0.41])
        self.axes.set_xlim([0,0.63])
        plt.pause(0.01)
        plt.draw()

    def readData(self, trackEvent, activeAlleys, texturePattern):
        self.annotate_activeAlleys(activeAlleys)
        self.texturePattern = texturePattern
        alley = int(trackEvent.split("-")[0])
        event = trackEvent.split("-")[1]
        x,y = self.alleyCoords[alley+1] #using subset alleys, offset where they starti
        self.x.append(x/10)
        self.y.append(y/10)
        if event == '0':
            self.alleyHistory[alley]['rewards'] += 1
        elif event == '7':
            self.alleyHistory[alley]['passes'] += 1

        self.updateGUI()

    def annotate_activeAlleys(self,activeAlleys):
        '''Given an array of alley states, display on plot
        which are active'''
        for i,a in enumerate(activeAlleys):
            if a in ['1', '3']:
                self.alleyStateColors[i+1] = 'r'

    def plotAlleyStats(self):
        for i, k in self.annLocs.items():
            self.axes.text(k[0], k[1], f"{self.texturePattern[i-1]}, {self.alleyHistory[i-1]['passes']}/{self.alleyHistory[i-1]['rewards']}",
                           color=self.alleyStateColors[i-1] ,size=24, rotation=k[2], ha='center', va='center')

    def plotTraj(self):
        for i in range(1):
            self.axes.plot(self.x,self.y)
            self.axes.set_ylim([0,0.5])
            self.axes.set_xlim([0,0.7])

    def updateGUI(self):
        self.axes.cla()
        self.drawRatterdam()
        self.plotTraj()
        self.plotAlleyStats()
        self.axes.set_ylim([0.05,0.42])
        self.axes.set_xlim([0.05,0.63])
        plt.pause(0.01)
        plt.draw()
        
    def mock(self):
        '''test some mock data'''
        a = [4,5,7,8,9,12,14,15,2,3,1,4,13,9,8,7]
        ac = [self.alleyCoords[i] for i in a]
        for i in ac:
            self.x.append(i[0]/10)
            self.y.append(i[1]/10)
            active = [7,8,14,3]
            self.plotTraj()
            self.annotate_activeAlleys(active)
            plt.pause(0.01)
            time.sleep(0.1)

    


