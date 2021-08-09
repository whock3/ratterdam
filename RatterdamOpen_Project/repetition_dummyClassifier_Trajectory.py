# -*- coding: utf-8 -*-
"""
Created on Tuesday July 27th ay 1:14pm 

@author: whockei1

Repetition Project
Creating and testing naive classifiers to establish performance benchmarks
for real classifier (e.g. RF) performance. Especially an issue with imbalanced
data like we have here. 

Classifier will be written by hand as sci-kit learn doesn't deal with a latent
group structure (i.e. track regions) afaik. 

Concern is there is data leak

This file is just for decoding trajectory 
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
from collections import Counter

#%% Load data, define defaults/data structs
#Load one unit, bc we are just looking at class label biases and any interaction
# between that and location

savepath = "E:\\Ratterdam\\repetition_decoding\\8-9-21_NaiveClassifierTrajectory_newrats\\"
alldata = []
timestamp = util.genTimestamp()
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
region_sets = {'RS1':[12],  #decode traj and dir
                'RS2':[3,5],  # decode traj and dir
                'RS3':[14,11] # decode traj and dir         
                }
    
for rat, day in zip(['R886', 'R886'],['D1', 'D2']):
    print(rat, day, "Beginning Naive Classifier")
    population, turns = RepCore.loadRecordingSessionData(rat, day)

            
    #%% Group data by trajectory 
    
    #Here's the logic. If you want the acivity from when the animal was on a given
    # alley on a certain pass, you take the ts entry of turn n to the ts exit of
    # turn n+1 and that corresponds to time spent on alley+. 
    currentDir, previousDir, nextDir = [], [], []
    currentAlley = []
    traj = []
    
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
        
        currentDir.append(cD)
        previousDir.append(pD)
        nextDir.append(nD)
        traj.append(f"{pD}{cD}{nD}")
        
    currentDir = np.asarray(currentDir)
    nextDir = np.asarray(nextDir)
    previousDir = np.asarray(previousDir)
    traj = np.asarray(traj)
    currentAlley = np.asarray(currentAlley)
    
    
    #%% Calculate trajectory biases 
    nruns = 1000
    
    for label, rs in region_sets.items():
        print(rat, day, label)
        naivePerfs = []
        
        for n in range(nruns):
            singleRunOutcomes = []
            
            for region in rs:
    
                # filter and count
                alleytrajs = traj[currentAlley==str(region)]
                c = Counter(alleytrajs)
                
                #normalize 
                for k,v in c.items():
                    c[k] = v/alleytrajs.shape[0]
                
            
                for t,real in enumerate(alleytrajs):
                    o=np.random.choice(list(c.keys()),size=1,p=list(c.values()))[0]
                    if o == real:
                        singleRunOutcomes.append(True)
                    else:
                        singleRunOutcomes.append(False)
                        
            naivePerfs.append(sum(singleRunOutcomes)/len(singleRunOutcomes))
            
        pct95 = np.percentile(naivePerfs, 95)
        alldata.append((rat, day, rs, naivePerfs))
        plt.figure(figsize=(15,12))
        plt.hist(naivePerfs,bins=25)
        plt.vlines(pct95, 0, 200)
        plt.title(f"{rat}{day} {label} Naive Stratified Traj. Decoder {nruns}x 95%ile {round(pct95,2)}", fontsize=16)
        plt.xlabel("Decoder performance, single run", fontsize=16)
        plt.ylabel("Frequency", fontsize=16)
        plt.savefig(savepath+f"{timestamp}_{rat}{day}_{label}_NaiveStratifiedTrajectory.png",dpi=300)
        
            
                
        
        
        











