# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:27:26 2021

@author: whockei1

Script to generate ratemaps from periods of stable winding number activity

Winding number (https://en.wikipedia.org/wiki/Winding_number)
is a measure of loops around a fixed observer point. We (IA/YB) calculate
WNs around each block of the maze and calculate their cumulative number
cw/ccw to give a sense of (one kind of) behavioral stereotypy

Q here is what, if any, relation is there between WN and repeating activity - 
i.e. can WN give insight into effect of stereotypy on repetition
"""

import numpy as np
import utility_fx as util
import os
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
import pandas as pd
from matplotlib import path

#%% Load data 
rat, day = 'R781', 'D4'
population, turns = RepCore.loadRecordingSessionData(rat, day)
cmap = util.makeCustomColormap()

# R781 D3, applies more or less to each block
# windingBounds = [0, 3055050000000.0, 3.0560*1e12, 3056500000000.0, 3057590000000.0, 3058120562236]

# R808 D6
# windingBounds = [0, 2.125*1e9, 3.46*1e9, 4.46*1e9]

#R808 D7
#windingBounds = [0, 2.61*1e9, 3.24*1e9, 5230154400]

df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
winding = np.load(df+f"{rat}{day}_winding.npz")
wdf = pd.DataFrame(winding['winding_data'])
ts = wdf['Timestamp']

#%% Run winding RMs

clustname = "TT5\\cl-maze1.2"
unit = population[clustname]
windingBounds = np.linspace(unit.position[0,0],unit.position[-1,0],num=12)


allbins = []
allrms = []
for b in range(len(windingBounds)-1):
    
    start, stop = windingBounds[b], windingBounds[b+1]
    spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)]
    pos = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=stop)]
    rm = util.makeRM(spikes,pos,bins=[40,60])
    allrms.append(rm)
    allbins.extend(rm.flatten())
    
mymax = np.nanpercentile(allbins, 98)
 
ncols = 3
fig, _ax = plt.subplots(int(np.ceil(len(windingBounds)/ncols)), ncols)

for b in range(len(windingBounds)-1):
    interval = [((windingBounds[b]-windingBounds[0])/1e6)/60, ((windingBounds[b+1]-windingBounds[0])/1e6)/60]
    rm = allrms[b]
    fig.axes[b].imshow(rm, aspect='auto', interpolation='None',origin='lower', cmap=cmap, vmax = mymax)
    fig.axes[b].set_title(f"{interval[0]} - {interval[1]}")
plt.suptitle(f"{rat}{day} {clustname}")    

plt.figure()
for i in range(1,7):
    plt.plot(wdf["Timestamp"], wdf[f"Block_{i}"],label=i)
plt.legend()
    
for b in windingBounds[1:]:
    plt.vlines(b,-50,20,color='k')
    
plt.title(f"{rat}{day} Winding Numbers")