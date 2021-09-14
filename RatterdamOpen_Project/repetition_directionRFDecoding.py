# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:49:30 2021

@author: whockei1
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

#%% Load data
for rat,day in zip(['R859','R859','R781','R781','R808','R808'],['D1','D2','D3','D4','D6','D7']):
    print(rat,day)
    
    ratborders = nab.loadAlleyBounds(rat, day)
    savepath = "E:\\Ratterdam\\repetition_decoding\\"
    datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    clustList, clustQuals = util.getClustList(datapath)
    population = {}
    qualThresh = 3
    
    for i,clust in enumerate(clustList):
        
        if clustQuals[i] >= qualThresh:
       
            print(clust)
            unit = RepCore.loadRepeatingUnit(datapath, clust, smoothing=1)
            
            # check to see if the overlapping regions for each field (i.e. each
            # entry in overlaps) contains >1 alley and if so delete the spikes
            #from that field.
            # repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
            # for p,ov in enumerate(overlaps):
            #     nalleys = sum([True if type(i)==int else False for i in ov]) # alleys are numeric, intersections alphabetical
            #     if nalleys >1:
            #         f = path.Path(unit.perimeters[p])
            #         o = f.contains_points(unit.spikes[:,1:])
            #         unit.spikes = unit.spikes[~o,:]
                    
                    
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
    alleyType = []
    currentDir, previousDir, nextDir = [], [], []
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
        if cD in ['0', '4', '6', '1', '12', '8', '15' '13', '10']:
            alleyType.append('H')
        else:
            alleyType.append('V')
        
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
    alleyType = np.asarray(alleyType)
        
    
    #%% Run random forest
    X = X[alleyType=='V',:]
    currentDir = currentDir[alleyType=='V']
    
    for target, label in zip([currentDir],["CurrentDirection"]):
        print(target)
        realoobs, shuffoobs = [], []
        nreps = 10
        nshuffs = 500
        ntrees = 1000
        
        # real
        for i in range(nreps):
            clf = RandomForestClassifier(n_estimators=ntrees, oob_score=True)
            clf.fit(X,target)
            realoobs.append(clf.oob_score_)
        
        # shuffle
        for i in range(nshuffs):
            Yshuff = np.random.permutation(target)
            clf = RandomForestClassifier(n_estimators=ntrees, oob_score=True)
            clf.fit(X,Yshuff)
            shuffoobs.append(clf.oob_score_)
            
        plt.figure(figsize=(8,6))
        shuff95 = round(np.percentile(shuffoobs,95),2)
        realmean = round(np.mean(realoobs),2)
        plt.title(f"{rat}{day} {label} - {nshuffs}s,{nreps}r,{ntrees}t, 95th:{shuff95},realmean:{realmean}")
        plt.xlabel("OOB Score", fontsize=16)
        plt.ylabel("Frequency", fontsize=16)
        plt.hist(shuffoobs, bins=25, color='k')
        plt.vlines(realmean,0, 100, color='r')
        plt.savefig(savepath+f"21-07-16_{rat}{day}_{label}_RFDecoding.png")
        plt.close()
        

    
    
    
    
    
    