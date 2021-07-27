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

Concern is there is data leak

This file is just for decoding directionality 
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
#Load one unit, bc we are just looking at class label biases and any interaction
# between that and location

rat = 'R859'
day = 'D2'

ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859}[rat]
savepath = "E:\\Ratterdam\\repetition_decoding\\naive_stratified_decoding\\"
datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(datapath)
population = {}
qualThresh = 3
    

unit = RepCore.loadRepeatingUnit(datapath, clustList[0], smoothing=1)       
pos, turns = alleyTrans.alleyTransitions(unit.position, ratborders, graph=False)
turns = pd.DataFrame(turns)
turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']

turns = pd.DataFrame(data=turns)
turns.dropna(inplace=True) 

#%% Set up classifier 
timestamp = util.genTimestamp(
    )
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
region_sets = {'RS1':[12],  #decode traj and dir
                'RS2':[3,5],  # decode traj and dir
                'RS3':[14,11], # decode traj and dir
                'RS4':[0,4,6,15,13,10], #decode dir (E-W)
                'RS5':[2,16,7,9]  # decode dir (N-S)
                }


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
    plt.figure(figsize=(12,12))
    plt.hist(naive_perf,color='k',bins=25)
    plt.vlines(np.percentile(naive_perf,95), 0, 200)
    plt.title(f"Naive Classifier (Stratified) for {rslabel}, {nreps}x trial outcome simulation", fontsize=18)
    plt.ylabel("Frequency", fontsize=16)
    plt.xlabel("Performance Guessing Label", fontsize=16)
    plt.savefig(savepath+f"{timestamp}_{rat}{day}_{rslabel}_NaiveStratifiedDecoding.png", dpi=300)
    plt.close()
    
    
#%% Testing 
nreps = 1000
rslabel = 'RS1'
rs = ['12']
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
plt.figure(figsize=(12,12))
plt.hist(naive_perf,color='k',bins=25)
plt.vlines(np.percentile(naive_perf,95), 0, 200)
plt.title(f"Naive Classifier (Stratified) for {rslabel}, {nreps}x trial outcome simulation", fontsize=18)
plt.ylabel("Frequency", fontsize=16)
plt.xlabel("Performance Guessing Label", fontsize=16)

    
        
        
        
        
        
        
        
        
        
        
        
        