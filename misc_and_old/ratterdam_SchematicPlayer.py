
"""
Created on Tue Apr 23 12:35:26 2019

@author: whockei1
"""

import numpy as np, socket, os, sys, cv2, re
from scipy.stats import sem
import matplotlib.colors as colors
from importlib import reload
import matplotlib.pyplot as plt
from PIL import Image

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.widgets import Button
from matplotlib import path



import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_CoreDataStructures as core
import ratterdam_Directionality as Dir
import ratterdam_visBasic as Vis
import ratterdam_Defaults as Def

import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef


ratPlot = Vis.BasicRateMaps()

cnames = ['black', 'blue', 'green', 'red', 'firebrick', 'cornflowerblue', 'orchid', 'darkcyan', 'midnightblue', 'saddlebrown', 'darkviolet', 'seagreen', 'indianred', 'goldenrod', 'orange', 'olive']
binWidth = wmDef.binWidth
cmap = util.makeCustomColormap()


class Unit():
    """
    Wrapper class because rep field ID algorithm looks
    for instance.spikes and instance.position
    """
    
    def __init__(self, s, p):
        self.spikes = s
        self.position = p
        self.fields = []
        self.colors = []
        self.smoothing = 1
        self.repUnit = RateMapClass.RateMap(self)
        
    def findFields(self):
        for i,pf in enumerate(self.repUnit.PF[:]):
            contour = path.Path(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)))
            PinC = self.position[contour.contains_points(self.position[:,1:])]
            posVisits = util.getVisits(PinC[:,0])
            field_FR = []
            field_TS = [] # take middle ts 
            for visit in posVisits:
                spk = self.spikes[np.logical_and(self.spikes[:,0] > visit[0], self.spikes[:,0] < visit[-1])]
                field_FR.append((spk.shape[0]/len(visit))*30)
                field_TS.append(((visit[0]-self.position[0,0])/1e6)/60)
            field_FR = util.weird_smooth(np.asarray(field_FR), self.smoothing)
            self.fields.append(np.column_stack((field_TS, field_FR)))



def plot_all_alleyBounds(ax):
    for alley in range(1,18):
        ratPlot.drawAlleyBounds(ax, alley)
        

class RatterdamAnimation():
    
    def __init__(self, unit, df):
        self.fig, self.ax = plt.subplots()
        self.increment_time = 5e6 # size of window
        self.offset_time = 1e6 # size of slide
        self.current_epoch_time = self.start_time # this refers to timestamp of the beginning of the currently plotted period of activity
        self.animation_finished = False
        self.paused = False
        self.play_spikes = True
        self.df = df
        # load clust
        self.position, self.spikes = unit.position, unit.spikes  
        
        
    def defineSessionBounds(self):
        with open(self.df+"sessionEpochInfo.txt","r") as f:
            lines = f.readlines()
            self.start_time, self.end_time = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
        
        
    def forward(self, event):
        self.current_epoch_time += 5e6
        
    def prev(self, event):
        self.current_epoch_time -= 5e6
        
    def myprint(self, event):
        print(self.current_epoch_time)
        
    def pause(self, event):
        if self.paused == True:
            self.paused = False
        else:
            self.paused = True
        
            
    def visualizeSession_scatter(self):
                
        while not self.animation_finished:    
            if not self.paused:

                self.ax.set_xlim([0,640])
                self.ax.set_ylim([0,480])
        
                plot_all_alleyBounds(self.ax)
                
                if self.play_spikes:
                
                    occs = util.getTsinInverval(self.position, self.current_epoch_time, self.current_epoch_time + self.increment_time)
                    occs = util.getPosFromTs(occs[:,0], self.position)
                    
                    spikes = util.getTsinInverval(self.spikes, self.current_epoch_time, self.current_epoch_time + self.increment_time)
                    spikes = util.getPosFromTs(spikes[:,0], self.position)
                    
                    self.ax.plot(occs[:,0], occs[:,1], color='k', alpha=0.7)
                    self.ax.scatter(spikes[:,0], spikes[:,1],s=56,c='r',alpha=0.7)
                self.ax.set_title(self.current_epoch_time,fontsize=24)
                
                if not self.paused:
                    self.current_epoch_time = self.current_epoch_time + self.offset_time
                
                if self.current_epoch_time >self.endtime: 
                    self.animation_finished = True
                
                plt.pause(0.1)
                self.ax.clear()
                
    
    def run(self):
        
        
        self.axpause = self.fig.add_axes([0.48, 0.05, 0.1, 0.075])
        bpause = Button(self.axpause, "Pause/Play")
        bpause.on_clicked(self.pause)
        
        self.axstore = self.fig.add_axes([0.59, 0.05, 0.1, 0.075])
        bstore = Button(self.axstore, "Print")
        bstore.on_clicked(self.myprint)
        
        self.axnext = self.fig.add_axes([0.81, 0.05, 0.1, 0.075])
        bnext = Button(self.axnext, 'Next')
        bnext.on_clicked(self.forward)
        
        self.axprev = self.fig.add_axes([0.7, 0.05, 0.1, 0.075])
        bprev = Button(self.axprev, 'Prev')
        bprev.on_clicked(self.prev)

        self.visualizeSession_scatter()
        
rat = "R808"
day = "D6"
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clust = "TT15\\cl-maze1.1"
p, s = util.loadRepeatingUnit(df, clust)
unit = Unit(s,p)
unit.findFields()
rattAm = RatterdamAnimation(unit, df)
rattAm.run()  
