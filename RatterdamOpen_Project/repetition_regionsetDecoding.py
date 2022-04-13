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
import copy
import traceback
import pickle 
import ratterdam_DataFiltering as Filt 

region_sets = {
                'RS6':[0,4,6,15,13,10,1,12,8],
                'RS7':[2,16,3,14,5,11,7,9]            
                }

alldata = {}

timestamp = util.genTimestamp()
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  

rat_list = ['R765',
            'R765',
            'R781', 
            'R781', 
            'R808', 
            'R808', 
            'R859', 
            'R859', 
            'R886', 
            'R886']

day_list = ['RFD5',
            'DFD4',
            'D3', 
            'D4',
            'D6',
            'D7',
            'D1',
            'D2',
            'D1',
            'D2']

# rat_list = ['R886']
# day_list = ['D1']

for rat,day in zip(rat_list,day_list):
    
    rewards = RepCore.readinRewards(rat, day)

    
    if rat not in alldata:
        alldata[rat] = {}
    alldata[rat][day] = {}
   
    print(f"Beginning decoding {rat} {day}...")
    if True:
        ratborders = nab.loadAlleyBounds(rat, day)
        savepath = "E:\\Ratterdam\\repetition_decoding\\2022-04-11_decoding\\"
        datapath = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
        
        population = superpop[rat][day]['units']
        turns = superpop[rat][day]['turns']
        refturns = superpop[rat][day]['refturns']
        
        #%% Create data arrays
        #Here's the logic. If you want the acivity from when the animal was on a given
        # alley on a certain pass, you take the ts entry of turn n to the ts exit of
        # turn n+1 and that corresponds to time spent on alley+. 
        currentDir, previousDir, nextDir, turnsIn, turnsOut = [], [], [], [], []
        currentAlley = []
        egoTurn = []
        traj = []
        X = np.empty((0, len(population)))
        
        for t, turn in turns.iterrows():
            if t< turns.shape[0]-1:
            
                turn_nm1 = refturns.iloc[t-1]         
                #turn_np1 = refturns.iloc[t+1]
                
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
                #nD = turn_np1['Allo+']
                
                # time interval below, based on alley+ of each turn, is how it's previously been done
                # for egocentric decoding, however, it becomes retrospective coding
                # start, stop = float(turn['Ts entry']), float(turn_np1['Ts exit'])
                
                start, stop =  float(turn['Ts entry']), float(refturns.iloc[t+1]['Ts exit'])
                
                #added 2-25-22 to check if there's a reward and not include the trial if so 
                isReward = np.where(np.asarray([(start < i < stop) for i in rewards])==True)[0]
                
                alley = turn['Alley+']
                border = ratborders.alleyInterBounds[alley]
                
                if isReward.shape[0]>0:
                    pass
                else:
                    
                    # cD= turn['Ego'] # currentDir value this turn
                    # pD = turn_nm1['Ego'] # previous ''
                    # nD = turn_np1['Ego'] # next ''
                    cD = turn['Allo+']
                    pD = turn['Allo-']
                    #nD = turn_np1['Allo+']
                    currentAlley.append(turn['Alley+']) # use this later to get visits to regions in a certain set 
                
    
                    duration = (stop-start)/1e6
                    
                    popvector = []
                    for unitname, unit in population.items():
                        
                        # 4/11/22: changing the feature variable to be normalized FR
                        #instead of FR. Normalized by session behavior through alley (and occ corrected too)
                        
                        #basic stuff, get the firing rate (# spikes / time on alley ) for each unit and append
                        # spike_count = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]
                        # rate = spike_count / duration
                        
                        normedRate = Filt.normalizeTrajectory(unit, start, stop, None, alley, border)
                        
                        popvector.append(normedRate)
                
                    popvector = np.asarray(popvector)
                    X = np.vstack((X, popvector))
                    currentDir.append(cD)
                    previousDir.append(pD)
                    #nextDir.append(nD)
                    #traj.append(f"{pD}{cD}{nD}")
                    turnsIn.append(f"{pD}{cD}")
                    
                    # this if statement is here because if the next turn, whose allo+
                    # would be the next turn direction, is not in ballistic turn idx
                    # then the turn featured a turn around and we don't want to conflate
                    # the next heading from a ballistic turn with the next heading where the animal
                    # turned around
                    #if t+1 in turns.index:            
                    #    turnsOut.append(f"{cD}{nD}")
                    egoTurn.append(turn['Ego'])
                    
                
        currentDir = np.asarray(currentDir)
        nextDir = np.asarray(nextDir)
        previousDir = np.asarray(previousDir)
        traj = np.asarray(traj)
        currentAlley = np.asarray(currentAlley)
        turnsIn = np.asarray(turnsIn)
        turnsOut = np.asarray(turnsOut)
        egoTurn = np.asarray(egoTurn)
        X = np.nan_to_num(X) # added 4/11/22 bc normalizing can give nans/infs 
        
        
        #%% Run random forest
        
        for regionsetlabel, regionset in [['RS6',region_sets['RS6']],['RS7',region_sets['RS7']]]:


            targets, targetlabels = [currentDir],['CurrentDirection']
            
            alldata[rat][day][regionsetlabel] = {}
           
                
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
                        
                    alldata[rat][day][regionsetlabel]['targets'] = [targetlabel,regionset]
                    alldata[rat][day][regionsetlabel]['oobs'] = {'Real':realoobs, 'Shuffle':shuffoobs}
                    
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
            
        print(f"Finished decoding {rat} {day}")
    # except Exception:
    #    print(traceback.format_exc())