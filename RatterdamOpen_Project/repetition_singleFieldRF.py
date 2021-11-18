# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:25:59 2021

@author: whockei1

Single unit decoding via random forest
Decoding direction, turn in/out, and trajectory at each field. 
Cannot pool fields for reasons that field properties like gain, sampling, etc
can bias things. Also turns and trajs are direction-gated so you get any encoding
above and beyond the directional signal (if any)
"""

import sklearn as skl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score

import numpy as np
import sys
import utility_fx as util
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import newAlleyBounds as nab
import repeatingPC as repPC
import math
import bisect
import pandas as pnd
from matplotlib import path
import copy
import pickle 
import repetition_decoderFunctions as Decoder
from sklearn import svm, preprocessing, metrics
from collections import Counter

 
with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superpopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)
    
rat, day = 'R859','D2'
turns = superpop[rat][day]['turns']
refturns = superpop[rat][day]['refturns']
pop = superpop[rat][day]['units']
ratborders = nab.loadAlleyBounds(rat, day)

#%% Parameters

params = {'ntrees':1000,'nreps':1,'nshuff':1000}
nsamples = 3
nbehaviors = 2

verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]
trackperimeter = [str(i) for i in [0,4,6,7,9,10,13,15,16,2]]
trackinterior = [str(i) for i in [1,3,14,12,5,11, 8]]

#sys.stdout = open("E:\\Ratterdam\\repetition_decoding\\21-10-14_decoding\\21-10-14_decoding.txt","w")

#%% Run 
print(rat, day)
for unitname, unit in pop.items():
    print(unitname)
    verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
    horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]
    
    for i,(perim, foverlap) in enumerate(zip(unit.perimeters, unit.overlaps)):
        try:
            print(f"Field {i}")
            for alleyfield in foverlap:
                
                #only looking at alleys right now
                if type(alleyfield)==int:
                    print(f"Alley {alleyfield}")
        
                    border = ratborders.alleyInterBounds[str(alleyfield)]
                    
                    contour = path.Path(perim)
                    field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                    field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                    
                    #break alley into number of segments, used to 2d hist each turn's activity
                    numSections = 4
                    if str(alleyfield) in verticals:
                        sections = np.linspace(border[1][0],border[1][1],num=numSections+1)
                        walls = border[0]
                    elif str(alleyfield) in horizontals:
                        sections = np.linspace(border[0][0],border[0][1],num=numSections+1)
                        walls = border[1]
                        
                        
                    currentDirection, nextDirection, previousDirection = [], [], []
                    turnIn, turnOut, trajectory = [], [], []
                    
                    responses = np.empty((0,numSections))
                    
                    # Get all behavioral events
                    for t, turn in turns.iterrows():
                        if int(turn['Alley+']) == alleyfield:
                            currd = Def.allocodedict[turn['Allo+']]
                            nd = Def.allocodedict[refturns.iloc[t+1]['Allo+']]
                            pd = Def.allocodedict[turn['Allo-']]
                            
                            currentDirection.append(currd)
                            nextDirection.append(nd)
                            previousDirection.append(pd)
                            
                            turnIn.append(f"{pd}{currd}")
                            turnOut.append(f"{currd}{nd}")
                            trajectory.append(f"{pd}{currd}{nd}")
                            
                            winStart, winEnd = float(turn['Ts entry']), float(refturns.iloc[t+1]['Ts exit'])
                            
                            winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                            winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                        
                            if str(alleyfield) in horizontals:
                                sh = np.histogram2d(winSpikes[:,2], winSpikes[:,1], bins=[walls,sections])
                                ph = np.histogram2d(winPos[:,2], winPos[:,1], bins=[walls,sections])
                            elif str(alleyfield) in verticals:
                                sh = np.histogram2d(winSpikes[:,2], winSpikes[:,1], bins=[sections, walls])
                                ph = np.histogram2d(winPos[:,2], winPos[:,1], bins=[sections, walls])
                                
                            
                            rate = sh[0]/ph[0]
                            responses = np.vstack((responses, rate.reshape(1,rate.shape[0])))
                                                    
                            
                    currentDirection = np.asarray(currentDirection)
                    nextDirection = np.asarray(nextDirection)
                    previousDirection = np.asarray(previousDirection)
                    
                    turnIn = np.asarray(turnIn)
                    turnOut = np.asarray(turnOut)
                    trajectory = np.asarray(trajectory)
                    
                    
                    # Now filter 
                    
                    # 'X' is a lost turn bc algorithm couldnt ID the turn properly.
                    goodidx = np.where(currentDirection!='X')[0]
                    currentDirection = currentDirection[goodidx]
                    nextDirection = nextDirection[goodidx]
                    previousDirection = previousDirection[goodidx]
                    turnIn = turnIn[goodidx]
                    turnOut = turnOut[goodidx]
                    trajectory = trajectory[goodidx]
                    responses = responses[goodidx]
                    
                    responses[np.where(~np.isfinite(responses))] = 0
                    responses = preprocessing.StandardScaler().fit_transform(responses)
                    
                    
                    
                    #%% Current Direction 
                    
                    
                    # Iterate over behaviors and include each with minimum # samples
                    # If the behavior type has sufficient num behaviors, decode it.
                    c = Counter(currentDirection)
                    validbehaviors = {}
                    for k,v in c.items():
                        if v >= nsamples:
                            validbehaviors[k] = v
                    if len(validbehaviors) >= nbehaviors:
                    
                        realCurrent = Decoder.runRandomForest(responses, 
                                                              currentDirection, 
                                                              **params)
                        
                        shuffCurrent = []
                        for i in range(params['nshuff']):
                            shuffCurrent.append(Decoder.runRandomForest(responses, 
                                                                        np.random.permutation(currentDirection), 
                                                                        **params))
                        
                        print(f"Real current direction decoding: {realCurrent}")
                        print(f"95th percentile shuffle current direction decoding: {np.nanpercentile(shuffCurrent,95)}")
                    #%% Turns in
                    
                    orientations = np.unique(currentDirection)
                    if len(orientations)>2:
                        print(f"{rat}{day} {unitname} {alleyfield} ERROR - too many current directions")
                        
                    for cdir in orientations:
                        
                        #filter by direction - we are trying to decode a signal beyond any directional signal
                        filtTurnIn = turnIn[currentDirection==cdir]
                        filtResponses = responses[currentDirection==cdir,:]
                        
                        #filter behaviors by threshold and only analyze if there's enough sampling
                        # within and between behaviors
                        c = Counter(filtTurnIn)
                        validbehaviors = {}
                        for k,v in c.items():
                            if v >= nsamples:
                                validbehaviors[k] = v
                        if len(validbehaviors) >= nbehaviors:
                            # keep valid behaviors only 
                            sampledResponses, sampledTurnIn = [], []
                            for r,s in zip(filtResponses, filtTurnIn):
                                if s in validbehaviors.keys():
                                    sampledResponses.append(r)
                                    sampledTurnIn.append(s)
                            
                            real = Decoder.runRandomForest(sampledResponses, 
                                                           sampledTurnIn, 
                                                           **params)
                            
                            shuff = []
                            for i in range(params['nshuff']):
                                shuff.append(Decoder.runRandomForest(sampledResponses, 
                                                                            np.random.permutation(sampledTurnIn), 
                                                                            **params))
                            
                            print(f"Real {cdir} turn in decoding: {real}")
                            print(f"95th percentile shuffle turn in decoding: {np.nanpercentile(shuff,95)}")
                            
                            
                    #%% Turns out
                    
                    orientations = np.unique(currentDirection)
                    if len(orientations)>2:
                        print("ERROR - too many current directions")
                        
                    for cdir in orientations:
                        print(cdir)
                        
                        #filter by direction - we are trying to decode a signal beyond any directional signal
                        filtTurnOut = turnOut[currentDirection==cdir]
                        filtResponses = responses[currentDirection==cdir,:]
                        
                        #filter behaviors by threshold and only analyze if there's enough sampling
                        # within and between behaviors
                        c = Counter(filtTurnOut)
                        validbehaviors = {}
                        for k,v in c.items():
                            if v >= nsamples:
                                validbehaviors[k] = v
                        if len(validbehaviors) >= nbehaviors:
                            # keep valid behaviors only 
                            sampledResponses, sampledTurnOut = [],[]
                            for r,s in zip(filtResponses, filtTurnOut):
                                if s in validbehaviors.keys():
                                    sampledResponses.append(r)
                                    sampledTurnOut.append(s)
                            
                            real = Decoder.runRandomForest(sampledResponses, 
                                                           sampledTurnOut, 
                                                           **params)
                            
                            shuff = []
                            for i in range(params['nshuff']):
                                shuff.append(Decoder.runRandomForest(sampledResponses, 
                                                                     np.random.permutation(sampledTurnOut), 
                                                                     **params))
                                
                            print(f"Real {cdir} turn out decoding: {real}")
                            print(f"95th percentile shuffle turn out decoding: {np.nanpercentile(shuff,95)}")
                            
                    
                    
                    
        except:
            pass