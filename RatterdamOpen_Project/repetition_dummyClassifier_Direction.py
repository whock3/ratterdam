# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 16:30:55 2021

@author: whockei1

Repetition Project
Creating and testing naive classifiers to establish performance benchmarks
for real classifier (e.g. RF) performance. Especially an issue with imbalanced
data like we have here. 

Classifier will be written by hand as sci-kit learn doesn't deal with a latent
group structure (i.e. track regions) afaik. 

Concern is there is data leak so decode based on other, uninteresting 
data and see how classifier performs

This file is just for decoding directionality 
"""


import sklearn as skl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score

import numpy as np
import utility_fx as util
import os, json 
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

#%% Load data
#Load one unit, bc we are just looking at class label biases and any interaction
# between that and location

timestamp = util.genTimestamp()

#region sets for 9-8-21: vertical vs horizontal alleys

region_sets = {'RS6':[0,4,6,15,13,10,1,12,8],
                'RS7':[2,16,3,14,5,11,7,9]
                }

codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
savepath = "E:\\Ratterdam\\repetition_decoding\\2022-03-23_decoding\\"

naive_perfs_datasets = {}

# rat_list = ['R765',
#             'R765',
#             'R781', 
#             'R781', 
#             'R808', 
#             'R808', 
#             'R859', 
#             'R859', 
#             'R886', 
#             'R886']

# day_list = ['RFD5',
#             'DFD4',
#             'D3', 
#             'D4',
#             'D6',
#             'D7',
#             'D1',
#             'D2',
#             'D1',
#             'D2']

rat_list = ['R886']
day_list = ['D1']

for rat,day in zip(rat_list,day_list):
    
    rewards = RepCore.readinRewards(rat, day)

    naive_perfs_datasets[f"{rat}{day}"] = {'RS6':None, 'RS7':None}
    ratborders = nab.loadAlleyBounds(rat, day)
    datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'    
    
    turns, unit = RepCore.loadTurns(rat, day)
    
    ballisticTurnIdx = []
    for t, turn in turns.iterrows():
        
        if t < turns.shape[0]-1:
        
        
            inter = turn['Inter']
            ts_start, ts_end = float(turn['Ts entry']), float(turns.iloc[t+1]['Ts exit'])
            
            isReward = np.where(np.asarray([(ts_start < i < ts_end) for i in rewards])==True)[0]
            
            # edit 10/2 removing check that last turn's inter wasnt the same,
            # i.e if alley- had a turnaround. since we are looking at things
            # in terms of alley+, only remove a turn if thats where a turnaround was
            if turn['Ego'] != '3' and turns.iloc[t+1].Inter != inter and isReward.shape[0]==0:
                ballisticTurnIdx.append(t)
    
    refturns = copy.deepcopy(turns) # keep a copy without filtering.
    turns = turns.iloc[np.asarray(ballisticTurnIdx)]
    
    #%% Set up classifier 
    
    dirbias = {str(i):{'N':0,'S':0,'E':0,'W':0} for i in range(17)}
    
    for alley in dirbias.keys():
        for d in ['1', '2', '3', '4']:
            count = sum(turns[turns['Alley+']==alley]['Allo+']==d) # expression in sum() gives boolean series if direction was d or not. sum gives count thereof
            dirbias[alley][codedict[d]] = count
        
    for alley in dirbias.keys():
        totalvisits = sum(dirbias[alley].values()) # total num passes in alley
        for d in ['N', 'S', 'E', 'W']:
            dirbias[alley][d] = dirbias[alley][d]/totalvisits # num visits in a dir / total visits, gets you prob
    
            
    #%% Run Naive Stratified Decoder
            
    nreps = 1000
    for rslabel, rs in region_sets.items():
        print(rslabel)
        naive_perf = []
        for n in range(nreps):
            outcomes = []
            for alley in rs:
                alley = str(alley)
                turnsubset = turns[turns['Alley+']==str(alley)]
                for _,turn in turnsubset.iterrows():
                    c=np.random.choice(['N','S','E','W'],size=1,p=[dirbias[alley]['N'],dirbias[alley]['S'],dirbias[alley]['E'],dirbias[alley]['W']])[0]
                    real = codedict[turn['Allo+']]
                    if real == c:
                        outcomes.append(1)
                    else:
                        outcomes.append(0)
            naive_perf.append(sum(outcomes)/len(outcomes))
        print(np.percentile(naive_perf,95))
        # plt.figure(figsize=(12,12))
        # plt.hist(naive_perf,color='k',bins=25)
        # plt.vlines(np.percentile(naive_perf,95), 0, 200)
        # plt.title(f"{rat}{day} {rslabel} Stratified Classifier, {nreps}x, 95th%ile = {round(np.percentile(naive_perf,95),2)}%", fontsize=18)
        # plt.ylabel("Frequency", fontsize=16)
        # plt.xlabel("Performance Guessing Label", fontsize=16)
        #plt.savefig(savepath+f"{timestamp}_{rat}{day}_{rslabel}_NaiveStratifiedDecoding.png", dpi=300)
        #plt.close()
        
        naive_perfs_datasets[f"{rat}{day}"][rslabel] = naive_perf 
        
    
       
with open(savepath+"naiveClassifierData.json", "w") as f:
    json.dump(naive_perfs_datasets, f)
        
        
        
        
        

        