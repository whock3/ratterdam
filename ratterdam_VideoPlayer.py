# -*- coding: utf-8 -*-
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

if socket.gethostname() == 'Tolman':
    codeDirBase = 'C:\\Users\\whockei1\\Google Drive'
elif socket.gethostname() == 'DESKTOP-BECTOJ9':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'

sys.path.insert(0, codeDirBase + '\\KnierimLab\\Ratterdam\\Code')
import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_CoreDataStructures as core
import ratterdam_Directionality as Dir
import ratterdam_visBasic as Vis
import ratterdam_Defaults as Def


ratPlot = Vis.BasicRateMaps()

def plot_all_alleyBounds(ax):
    for alley in range(1,18):
        ratPlot.drawAlleyBounds(ax, alley)
        

class RatterdamAnimation():
    
    def __init__(self, vfile, clustname):
        self.fig, self.ax = plt.subplots()
        self.df = vfile
        self.start_time, _ = util.loadRatterdamEpochTimes(vfile)
        self.increment_time = 5e6 # size of window
        self.offset_time = 1e6 # size of slide
        self.current_epoch_time = self.start_time # this refers to timestamp of the beginning of the currently plotted period of activity
        self.animation_finished = False
        self.paused = False
        self.play_spikes = True
        self.cap = cv2.VideoCapture(vfile)
        # load clust
        self.position, self.spikes = util.loadRepeatingUnit(self.df, clustname)
        
    
        
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
            
    def parseSMIFile(self):
        """
        SMI file is generic (i.e. used beyond NL) type to store subtitles
        Neuralynx stores the NL ts of a video frame as subtitles
        So extract them into a np array
        """
        with open(self.df+"VT1.smi","r") as file:
            lines = file.readlines()
        tsLines = [line for line in lines if 'SYNC' in line] # SYNC tag being the unique one for ts data AFAIK
        pattern = re.compile("<.*?>")
        NLts = [int(re.sub(pattern,'',entry).strip()) for entry in tsLines]
        self.NLts = np.array(NLts)
        
            
    def getCurrentVideoFrame(self):
        start, end = self.current_epoch_time, self.current_epoch_time + self.increment_time
        idx, = np.where((self.NLts > start)&(self.NLts < end))
        print(self.NLts[:])
        self.cap.set(1, idx[-1])
        ret, frame = self.cap.read()
        frame = np.flip(frame, axis=1)
        pic = Image.fromarray(frame)
        pic.thumbnail((640,480))
        frame = np.asarray(pic)
        return frame
        
            
    def visualizeSession_scatter(self):
        
        self.parseSMIFile()
        
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
                
                frame = self.getCurrentVideoFrame()
                self.ax.imshow(frame[:,:,2],origin='upper')
                if not self.paused:
                    self.current_epoch_time = self.current_epoch_time + self.offset_time
                
                if self.current_epoch_time > 1907952676606: #position[-1,0]:
                    self.animation_finished = True
                
                plt.pause(0.1)
                self.ax.clear()
                    
    def visualizeSession_rm(self):
        while not self.animation_finished:  
            if not self.paused:
                self.ax.set_xlim([0,640])
                self.ax.set_ylim([0,480])

                plot_all_alleyBounds(self.ax)

                occs = util.getTsinInverval(self.position, self.current_epoch_time, self.current_epoch_time + self.increment_time)
                occs = util.getPosFromTs(occs[:,0], self.position)

                spikes = util.getTsinInverval(self.spikes, self.current_epoch_time, self.current_epoch_time + self.increment_time)
                spikes = util.getPosFromTs(spikes[:,0], self.position)

                self.ax.plot(occs[:,0], occs[:,1], color='k', alpha=0.7)
                self.ax.scatter(spikes[:,0], spikes[:,1],s=56,c='r',alpha=0.7)
                self.ax.set_title(self.current_epoch_time,fontsize=24)

                self.current_epoch_time = self.current_epoch_time + self.offset_time

                if self.current_epoch_time > self.position[-1,0]:
                    self.animation_finished = True

                plt.pause(0.05)
                self.ax.clear()
                
    
    def run(self):
        
        
        self.axpause = self.fig.add_axes([0.48, 0.05, 0.1, 0.075])
        bpause = Button(self.axpause, "Pause")
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

vfile =  'E:\\Ratterdam\\R859\\R859_RatterdamOpen_D3\\'
clustname = 'TT1_0001\\cl-maze1.1'
rattAm = RatterdamAnimation(vfile, clustname)
rattAm.run()  
