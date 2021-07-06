# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 13:36:04 2021

@author: whockei1

Aggregating data based on turns to visualize intuitively any relation between
firing and direction/turn/trajectory.

"""

#%% Imports
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
import math
import bisect
import pandas as pd

#%% Setup
rattrans = alleyTrans.R781
pos, turns = alleyTrans.alleyTransitions(unit.position, rattrans, graph=False)
turns = pd.DataFrame(turns)
turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']

turns = pd.DataFrame(data=turns)
turns.dropna(inplace=True) 


#%% Turn based rate maps - setup
# Group data by allo dir pre (4 plots) and allo dir post (4 plots) +/- 1.5 turns
#(the half turn is bc the turns end at the intersection so punch out a bit more
# in time to get the full +/1 1 turn)
#Reminder about code order: N,E,S,W,  F,R,B,L


turnRMS = {'Pre':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
                  }, 
           
           'Post':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
                  }
           }

for field in unit.fields:
    for visit in field:
        turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
        try:
            # first/last turn has no previous/next turn, so just ignore it
            if turnIdx > 2 and turnIdx < turns.shape[0]-2:
                spikesPre = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
                spikesPost = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

                turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'],spikesPre))
                turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'],spikesPost))
    
                occPre = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
                occPost = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

                turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'],occPre))
                turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'],occPost))
        except:
            print("Likely invalid code found")

#%% Plot turn rms
for d in ['Pre', 'Post']:
    fig, ax = plt.subplots(2,2)
    for i,code in enumerate([('1','N'),('2','E'),('3','S'),('4','W')]):
        s,o = turnRMS[d][code[0]]['Spikes'], turnRMS[d][code[0]]['Pos']
        rm = util.makeRM(s,o)
        fig.axes[i].imshow(rm,aspect='auto',interpolation='None',cmap=cmap,vmax=7,origin='lower')
        fig.axes[i].set_title(code[1])
    plt.suptitle(f"{d}-turn Bearing")
    
    




