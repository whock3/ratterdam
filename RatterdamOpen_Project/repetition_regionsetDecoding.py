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
                'RS5':[2,16,7,9]  # decode dir (N-S)
                }

alldata = []

timestamp = util.genTimestamp()

for rat,day in zip(['R859','R859','R781','R781','R808','R808'],['D1','D2','D3','D4','D6','D7']):
    print(rat,day)

    ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859}[rat]
    savepath = "E:\\Ratterdam\\repetition_decoding\\"
    datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    clustList, clustQuals = util.getClustList(datapath)
    population = {}
    qualThresh = 3
    
    for i,clust in enumerate(clustList):
        
        if clustQuals[i] >= qualThresh:
       
            print(clust)
            unit = RepCore.loadRepeatingUnit(datapath, clust, smoothing=1)                                   
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
            
            
    # Session endpoints data
    with open(datapath+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
        start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
        
    nepoch=3
    intervals = np.linspace(start,end,nepoch+1)
    
    #%% Create turn df 
    pos, turns = alleyTrans.alleyTransitions(unit.position, ratborders, graph=False)
    turns = pd.DataFrame(turns)
    turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']
    
    turns = pd.DataFrame(data=turns)
    turns.dropna(inplace=True) 
    #%% Create data arrays
    #Here's the logic. If you want the acivity from when the animal was on a given
    # alley on a certain pass, you take the ts entry of turn n to the ts exit of
    # turn n+1 and that corresponds to time spent on alley+. 
    currentDir, previousDir, nextDir = [], [], []
    currentAlley = []
    traj = []
    X = np.empty((0, len(population)))
    
    for t in range(1,turns.shape[0]-1):
        
        turn_nm1 = turns.iloc[t-1]
        turn = turns.iloc[t]
        turn_np1 = turns.iloc[t+1]
                
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
        
    currentDir = np.asarray(currentDir)
    nextDir = np.asarray(nextDir)
    previousDir = np.asarray(previousDir)
    traj = np.asarray(traj)
    currentAlley = np.asarray(currentAlley)
    
        
    
    #%% Run random forest
    
    for regionsetlabel, regionset in region_sets.items():
        if regionsetlabel in ['RS1','RS2','RS3']:
            targets, targetlabels = [traj, currentDir], ['Trajectories', 'CurrentDir']
        else:
            targets, targetlabels = [currentDir], ['CurrentDir']
            
        for target, targetlabel in zip(targets, targetlabels):
            print(targetlabel, regionsetlabel)
            subset =[int(i) in regionset for i in currentAlley]
            Xsubset, targetsubset = X[subset,:], target[subset]
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
                
            alldata.append((rat, day, regionsetlabel, realoobs, shuffoobs))
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
            
        