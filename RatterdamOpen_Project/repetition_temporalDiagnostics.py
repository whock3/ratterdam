# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 16:56:04 2021

@author: whockei1

Script to generate diagnostic plots of raw neural and behavioral data
over time for each unit.

"""

#%% Imports
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
import math
import bisect
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import repeatingPC as repPC
import copy
import more_itertools

#%% Load data and defaults

rat, day = 'R859', 'D1'
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\timeDiagnostics\\'
if not os.path.exists(savepath):
    os.makedirs(savepath)
    
population, turns = RepCore.loadRecordingSessionData(rat, day)
rewards = RepCore.readinRewards(rat, day)
ratborders = nab.loadAlleyBounds(rat, day)
allocodedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
egocodedict = {'1':'S','2':'R','3':'B','4':'L','0':'X'}
cmap = util.makeCustomColormap()

# Remove turnarounds/pivots
ballisticTurnIdx = []
for i in range(1,turns.shape[0]-1):
   row = turns.iloc[i]
   inter = row['Inter']
   if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter and turns.iloc[i-1].Inter != inter:
       ballisticTurnIdx.append(i)

refturns = copy.deepcopy(turns) # keep a copy without filtering.
turns = turns.iloc[np.asarray(ballisticTurnIdx)]

#%% 
for unitname, unit in population.items():
    
    print(f"Starting {unit.name}")
    u = unit.name.split('\\')[0]+unit.name.split('\\')[1]

    repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
    #making these plots for all fields so combine all the overlaps
    alloverlaps = []
    for i in overlaps:
        for j in i:
            alloverlaps.append(str(j))
            
    # this is redundant in loop as all units share position data. not a huge deal. 
    start, stop = unit.position[0,0], unit.position[-1,0]
    # 6 windows is 10mins per window assuming an hour session. pass 7 bc we want 6 intervals so give 7 endpoints
    windowBounds = np.linspace(start, stop, 7)
    
    fig, ax = plt.subplots(6,4,figsize=(13,15))
    
    for i,(wA, wB) in enumerate(zip(windowBounds[:-1], windowBounds[1:])):
        # get windowed data 
        wpos = unit.position[(unit.position[:,0]>wA)&(unit.position[:,0]<=wB)]
        wspk = unit.spikes[(unit.spikes[:,0]>wA)&(unit.spikes[:,0]<=wB)]
        wrm = util.makeRM(wspk, wpos,bins=[50,70])
        wturns = turns[(turns['Ts exit'].astype(float)>wA)&(turns['Ts exit'].astype(float)<=wB)]
        wr = rewards[(rewards>wA)&(rewards<=wB)]
        
        
        # matrix to count different types of turns within time window. E.g matrix['N']['W'] is the number of north to west turns
        alloTurnMatrix = {i:{j:0 for j in ['N','E','S','W','X']} for i in ['N','E','S','W','X']}
        egoTurnMatrix = {i:0 for i in ['S','R','B','L','X']}
        
        for _, turn in wturns.iterrows():
            # only want turns that involve the cell's fields 
            if turn['Alley-'] in alloverlaps or turn['Alley+'] in alloverlaps or turn['Inter'] in alloverlaps:
                alloTurnMatrix[allocodedict[turn['Allo-']]][allocodedict[turn['Allo+']]] += 1
                egoTurnMatrix[egocodedict[turn['Ego']]] += 1
         
        
        ax[i,0].imshow(wrm, origin='lower',aspect='auto',interpolation='None',cmap=cmap)
        ax[i,0].set_title(f"Ratemap {round(((wA-windowBounds[0])/1e6)/60,2)}-{round(((wB-windowBounds[0])/1e6)/60,2)}mins")
         
        alloTurnMatrix = pd.DataFrame(alloTurnMatrix)
        ax[i,1].imshow(alloTurnMatrix)
        ax[i,1].set_xticks(range(5))
        ax[i,1].set_yticks(range(5))
        ax[i,1].set_xticklabels(['N','E','S','W','X'])
        ax[i,1].set_yticklabels(['N','E','S','W','X'])
        ax[i,1].set_xlabel("Pre-turn")
        ax[i,1].set_ylabel("Post-turn")
        for (q, j), z in np.ndenumerate(alloTurnMatrix):
            ax[i,1].text(j, q, z, ha='center', va='center')
        ax[i,1].set_title("Allocentric Turn Counts")
            
        ax[i,2].bar(range(len(egoTurnMatrix.keys())),egoTurnMatrix.values())
        ax[i,2].set_xticks(range(5))
        ax[i,2].set_xticklabels(egoTurnMatrix.keys())
        ax[i,2].set_title("Egocentric Turn Counts")    
        
        util.drawTrack(ax=ax[i,3],rat=rat,day=day)
        for r in wr:
            rts = util.takeClosest(unit.position[:,0],r)
            rp=unit.position[np.where(unit.position[:,0]==rts)[0][0]]
            ax[i,3].scatter(rp[1],rp[2],c='r',s=25,zorder=99)
        ax[i,3].set_title("Rewards in Window")
    
    
    plt.suptitle(f"{rat}{day} {unit.name}")
    plt.tight_layout()
    plt.savefig(savepath+u+"_timeDiagnostics.pdf",dpi=300)
    plt.close()
    
    print(f"Finished {unit.name}")