# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 12:49:40 2021

@author: whockei1

Script to compute population ratemap correlations over time


Method 1: is similar to Leutgeb 2007 fig 3a. ("pattern sep in DG and CA3")
basically create rms in time window. stack them across population. correlate
pixel stack same location across time windows
"""

import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import newAlleyBounds as nab
import repeatingPC as repPC
import math
import bisect
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib import path
import copy
import traceback
import numpy.ma as ma

rat, day = 'R859', 'D2'
ratborders = nab.loadAlleyBounds(rat, day)
datapath = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
population, turns = RepCore.loadRecordingSessionData(rat, 
                                                     day,
                                                     activityThreshType='spikes',
                                                     activityThresh=50)


#%% Method 1: bin-wise population vector correlation over time like leutgeb 2007
unit = population[list(population.keys())[0]] # get a random unit for its positiuon data
w = 10*60*1e6
s = 5*60*1e6
#make a bunch of windows then trim to length of session 
wb = [[unit.position[0,0]+(i*s),unit.position[0,0]+w+(i*s)] for i in range(100)]
windows = [i for i in wb if i[1] < unit.position[-1,0]]

#windows = np.linspace(u.position[0,0], u.position[-1,0],num=13)

allRMStacks = []

for winA, winB in windows:
    
    rms = []
    for unitname, unit in population.items():
        
        behavior = unit.position[(unit.position[:,0]>winA)&(unit.position[:,0]<=winB)]
        spikes = unit.spikes[(unit.spikes[:,0]>winA)&(unit.spikes[:,0]<=winB)]
        rm = util.makeRM(spikes, behavior, bins=[30, 50])
        rms.append(rm)
        
    rms = np.asarray(rms)
    
    allRMStacks.append(rms)
    
avgcorrs = []
for a,b in zip(allRMStacks[:-1], allRMStacks[1:]):
    corrs = []
    for i in range(a.shape[1]):
        for j in range(a.shape[2]):
            corrs.append(np.corrcoef(a[:,i,j], b[:,i,j])[0,1])
    corrs = np.asarray(corrs)
    avgcorrs.append(np.mean(corrs[~np.isnan(corrs)]))
        

#%% correlate individual cells' ratemaps over time
cellcorrs = []
for i in range(len(population.keys())):
    
    cellcorr = []
    
    for a,b in zip(allRMStacks[:-1], allRMStacks[1:]):

        maska = ma.masked_invalid(a[i,:,:].flatten())
        maskb = ma.masked_invalid(b[i,:,:].flatten())
        m = (~maska.mask & ~maskb.mask)
        cellcorr.append(ma.corrcoef(a[i,:,:].flatten()[m], b[i,:,:].flatten()[m]).data[0,1])
        
    cellcorrs.append(cellcorr)
        
    
#%% Plot

fig, ax = plt.subplots(8,6)
for i,c in enumerate(cellcorrs):
    fig.axes[i].plot(util.weird_smooth(np.asarray(c),1),marker='.')
    fig.axes[i].set_title(population[list(population.keys())[i]].name)
    fig.axes[i].set_ylim([0,1])

plt.suptitle(f"{rat} {day} Windowed Within Cell Correlations")
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()