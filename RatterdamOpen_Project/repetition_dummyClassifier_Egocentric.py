# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 16:31:07 2021

@author: whockei1

Naive classifier for egocentric turn decoding

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


codedict = {'1':'S','2':'R','3':'B','4':'L','0':'X'}

region_sets = {
                'RS6':[0,4,6,15,13,10,1,12,8],
                'RS7':[2,16,3,14,5,11,7,9],
                'RS8':[1,3,14,12,5,11,8]
                }
timestamp = util.genTimestamp()
savepath = 'E:\\Ratterdam\\repetition_decoding\\21-09-14_decoding\\'

for rat, day in zip(['R765', 'R781', 'R781', 'R808','R808', 'R859', 'R859', 'R886', 'R886'], ['RFD5', 'D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2']):
   
    turns, unit = RepCore.loadTurns(rat, day)
    
    #%% Calculate biases 
    
    dirbias = {str(i):{'S':0, 'R':0, 'B':0, 'L':0} for i in range(17)}
    
    for alley in dirbias.keys():
        for d in ['1', '2', '3', '4']:
            count = sum(turns[turns['Alley+']==alley]['Ego']==d) # expression in sum() gives boolean series if direction was d or not. sum gives count thereof
            dirbias[alley][codedict[d]] = count
        
    
    for alley in dirbias.keys():
        totalvisits = sum(dirbias[alley].values()) # total num passes in alley
        for d in ['S', 'R', 'B', 'L']:
            dirbias[alley][d] = dirbias[alley][d]/totalvisits # num visits in a dir / total visits, gets you prob
    
    #%% Run classifier 
    
    nreps = 1000
    for rslabel, rs in [['RS8', region_sets['RS8']]]:
        print(rslabel)
        naive_perf = []
        for n in range(nreps):
            outcomes = []
            for alley in rs:
                alley = str(alley)
                turnsubset = turns[turns['Alley+']==str(alley)]
                for _,turn in turnsubset.iterrows():
                    c=np.random.choice(['S','R','B','L'],size=1,p=[dirbias[alley]['S'],dirbias[alley]['R'],dirbias[alley]['B'],dirbias[alley]['L']])[0]
                    real = codedict[turn['Ego']]
                    if real == c:
                        outcomes.append(1)
                    else:
                        outcomes.append(0)
            naive_perf.append(sum(outcomes)/len(outcomes))
        print(np.percentile(naive_perf,95))
        plt.figure(figsize=(12,12))
        plt.hist(naive_perf,color='k',bins=25)
        plt.vlines(np.percentile(naive_perf,95), 0, 200)
        plt.title(f"{rat}{day} {rslabel} Egocentric Stratified Classifier, {nreps}x, 95th%ile = {round(np.percentile(naive_perf,95),2)}%", fontsize=18)
        plt.ylabel("Frequency", fontsize=16)
        plt.xlabel("Performance Guessing Label", fontsize=16)
        plt.savefig(savepath+f"{timestamp}_{rat}{day}_{rslabel}_NaiveStratifiedEgocentricDecoding.png", dpi=300)
        plt.close()