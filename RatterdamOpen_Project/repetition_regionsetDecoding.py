# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 13:03:34 2021

@author: whockei1

Ratterdam repetition - topological region decoding
Break track into sets of regions. All regions within a set share the same
connectivity. Goal is to decode trajectory from each region, based on what
types of trajectories are available at each region 
"""

import sklearn as skl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score


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

# region sets for 7-19-21 decoding
# region_sets = {'RS1':[1,3,5,14,11,8],
#                'RS2':[4,13],
#                'RS3':[2,0,16,15,10,9,7,6],
#                'RS4':['F','G'],
#                'RS5':['A','D','I','L'],
#                'RS6':['B','C','J','K']
#                }

#region sets for 7-20-21 decoding
# region_sets = {'RS1':[3,5,14,11],
#                 'RS2':[1,8],
#                 'RS3':[2,0,16,15,10,9,7,6,4,13],
#                 'RS4':[12]
#                 }

#region sets for 7-21-21 decoding
region_sets = {'RS1':[12],  #decode traj and dir
                'RS2':[3,5],  # decode traj and dir
                'RS3':[14,11], # decode traj and dir
                'RS4':[0,4,6,15,13,10], #decode dir (E-W)
                'RS5':[2,16,7,9],  # decode dir (N-S)
                'RS6':[0,4,6,15,13,10,1,12,8],
                'RS7':[2,16,3,14,5,11,7,9]
                }

alldata = []

timestamp = util.genTimestamp()
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}


for rat,day in zip(['R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886'],['D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2']):
    print(rat,day)

    ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859, 'R765':nab.R765, 'R886':nab.R886}[rat]
    savepath = "E:\\Ratterdam\\repetition_decoding\\8-23-21_directionDecoding\\"
    datapath = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    
    # Filter to remove turn-around trajectories. These cause same label
    # to map onto different behaviors
    ballisticTurnIdx = []
    for i in range(1,turns.shape[0]-1):
        row = turns.iloc[i]
        inter = row['Inter']
        
        #logic checking turn arounds: code 3 for ego means within the turn he turned around. e.g. 13-j-13. ignore it.
        # for turn arounds that span two turns (most of them), think of it like the turn around consists of a turn in
        # and then a turn back out. each 'leg' of the turnaround has an entry associated w it in the turns db
        # so as we iter over the db we check 1 ahead and 1 behind to make sure were not part of a turn around currently.
        # caveat is you lose a 'straight-through' leg if its part of a turn around (.e.g the leg out he traverses through
        # an alley) and this could theoretically be used in directionality decoding
        if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter and turns.iloc[i-1].Inter != inter:
            ballisticTurnIdx.append(i)
            
    turns = turns.iloc[ballisticTurnIdx]
    
    #%% Create data arrays
    #Here's the logic. If you want the acivity from when the animal was on a given
    # alley on a certain pass, you take the ts entry of turn n to the ts exit of
    # turn n+1 and that corresponds to time spent on alley+. 
    currentDir, previousDir, nextDir, turnsIn, turnsOut = [], [], [], [], []
    currentAlley = []
    traj = []
    X = np.empty((0, len(population)))
    
    for t in range(1,turns.shape[0]-1):
        
        turn_nm1 = turns.iloc[t-1]
        turn = turns.iloc[t]
        turn_np1 = turns.iloc[t+1]
        
        # ballisticTurns removes turn-around turns. But since trajectory
        # labels depend on the n-1, n+1 turns, the turns adjacent to the removed
        # turns must also be ignored, at an unfortunate loss of data
        
        # if turn['Alley+'] == turn_np1['Alley-']:
        #     pass
                
        # cD= turn['Ego'] # currentDir value this turn
        # pD = turn_nm1['Ego'] # previous ''
        # nD = turn_np1['Ego'] # next ''
        cD = turn['Allo+']
        pD = turn['Allo-']
        nD = turn_np1['Allo+']
        currentAlley.append(turn['Alley+']) # use this later to get visits to regions in a certain set 
        
        start, stop = float(turn['Ts entry']), float(turn_np1['Ts exit'])
        duration = (stop-start)/1e6
        
        popvector = []
        for unitname, unit in population.items():
            #basic stuff, get the firing rate (# spikes / time on alley ) for each unit and append
                spike_count = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]
                rate = spike_count / duration
                popvector.append(rate)
    
        popvector = np.asarray(popvector)
        X = np.vstack((X, popvector))
        currentDir.append(cD)
        previousDir.append(pD)
        nextDir.append(nD)
        traj.append(f"{pD}{cD}{nD}")
        turnsIn.append(f"{pD}{cD}")
        turnsOut.append(f"{cD}{nD}")
        
    currentDir = np.asarray(currentDir)
    nextDir = np.asarray(nextDir)
    previousDir = np.asarray(previousDir)
    traj = np.asarray(traj)
    currentAlley = np.asarray(currentAlley)
    turnsIn = np.asarray(turnsIn)
    turnsOut = np.asarray(turnsOut)
    
        
    
    #%% Run random forest
    
    for regionsetlabel, regionset in region_sets.items():
        if regionsetlabel in ['RS6', 'RS7']:
            targets, targetlabels = [currentDir], ['CurrentDirection']
       
            
            for target, targetlabel in zip(targets, targetlabels):
                
                print(targetlabel, regionsetlabel)
                subset =[int(i) in regionset for i in currentAlley]
                Xsubset, targetsubset = X[subset,:], target[subset]
                        
                if Xsubset.shape[0] > 5:
                
                    realoobs, shuffoobs = [], []
                    nreps = 1
                    nshuffs = 1000
                    ntrees = 1000
                    
                    # real
                    for i in range(nreps):
                        clf = RandomForestClassifier(n_estimators=ntrees, oob_score=True)
                        clf.fit(Xsubset,targetsubset)
                        realoobs.append(clf.oob_score_)
                    
                    # shuffle
                    for i in range(nshuffs):
                        Yshuff = np.random.permutation(targetsubset)
                        clf = RandomForestClassifier(n_estimators=ntrees, oob_score=True)
                        clf.fit(Xsubset,Yshuff)
                        shuffoobs.append(clf.oob_score_)
                        
                    alldata.append((rat, day, targetlabel, regionsetlabel, realoobs, shuffoobs))
                    plt.figure(figsize=(8,6))
                    shuff95 = round(np.percentile(shuffoobs,95),2)
                    realmean = round(np.mean(realoobs),2)
                    plt.title(f"{rat}{day} {targetlabel} {regionsetlabel} - {nshuffs}s,{nreps}r,{ntrees}t, 95th:{shuff95},realmean:{realmean}")
                    plt.xlabel("OOB Score", fontsize=16)
                    plt.ylabel("Frequency", fontsize=16)
                    plt.hist(shuffoobs, bins=25, color='k')
                    plt.vlines(realmean,0, 100, color='r')
                    plt.savefig(savepath+f"{timestamp}_{rat}{day}_{targetlabel}_{regionsetlabel}_RFDecoding.png")
                    plt.close()
                    
                else:
                    print(f"{rat}{day} {targetlabel} had insufficient sampling for decoder")
            
    