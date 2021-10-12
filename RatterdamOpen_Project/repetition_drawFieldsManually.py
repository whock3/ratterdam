# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 15:00:59 2021

@author: whockei1

Repetition
Pre-processing and manual curation
Re-drawing place fields that detection algorithm drew poorly
Draw ratemap with original field perims from alg, connect event
handler to return clicks and click around field. Return a list that can
be saved to a file in the recording day directory.
"""

#%% 
from matplotlib import pyplot as plt
import numpy as np
import ratterdam_RepetitionCoreFx as RepCore
import utility_fx as util
from matplotlib import path
import williamDefaults as wmDef
import os, json


rat = 'R886'
day = 'D1'
df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustname = "TT10\\cl-maze1.3"

cmap = util.makeCustomColormap()

unit = RepCore.loadRepeatingUnit(rat, day, clustname)

class FieldDrawer():
    
    def __init__(self, df,unit):
        self.df = df
        self.unit = unit
        self.coords = []
        self.fields = {}

    def onclick(self, event):
        b = f"{event.button}"
        if b == '1':
            self.coords.append([event.xdata, event.ydata])
            plt.scatter(event.xdata, event.ydata, c='r',marker='*',s=20,zorder=99)
            self.fig.canvas.draw()
            
    
    def addField(self, fnum):
        self.fields[fnum] = self.coords + [self.coords[0]] # close the border
        self.coords = []
        
    
    def clearField(self):
        self.coords = []
        
        
    def viewPosInField(self,border):
        contour = path.Path(border)
        pinc = self.unit.position[contour.contains_points(self.unit.position[:,1:])]
        plt.scatter(pinc[:,1],pinc[:,2])
        
    def saveFields(self):
        
        if 'manuallyRedrawnFields' not in os.listdir(df):
            os.mkdir(self.df+"manuallyRedrawnFields\\")
        
        u = self.unit.name.split('\\')[0]+self.unit.name.split('\\')[1]
    
        with open(self.df+f"manuallyRedrawnFields\\{u}_manuallyRedrawnFields.json", "w") as f:
            json.dump(self.fields, f)
    
    
    def drawFields(self):
    
        self.fig, self.ax = plt.subplots(figsize=(8,6))
        self.ax.imshow(self.unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                      cmap=cmap, vmax=np.nanpercentile(self.unit.repUnit.rateMap2D, 98),
               extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
        self.ax.set_title(f"{self.unit.name}, cutoff = {round(np.nanpercentile(self.unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
        self.ax.axis('equal')
        self.ax.set_ylim([0,480])
        self.ax.set_xlim([0, 640])
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
    
        for i, field in enumerate(self.unit.smoothedFields):       
                self.ax.plot(self.unit.perimeters[i][:,0], self.unit.perimeters[i][:,1],color=self.unit.colors[i], label=f"Field {i}")
        self.ax.legend()

        self.fig.canvas.mpl_connect("button_press_event",self.onclick)


    
# def onpress(self, event):
#     """
#     Key press controls for simple field drawing 
#     c - clear current field list
#     x - save into existing list of lists 
#     z - print current list of field coords and draw 
#     """
#     k = f"{event.key}"
#     if k == 'c':
#         self.field = []
#     if k == 'x':
#         allfields.append(coords)
#     if k == 'z':
#         c = np.asarray(coords)
#         plt.plot(c[:,0], c[:,1], c='k',linewidth=2)
    
    
