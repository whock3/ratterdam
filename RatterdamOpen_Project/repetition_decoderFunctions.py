# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 13:14:37 2021

@author: whockei1

Functions to run decoders, both RF and naive classifier
These will be run in a separate batch script
"""

import numpy as np, pandas as pd, copy
from sklearn.ensemble import RandomForestClassifier
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import utility_fx as util
import json, os 
from matplotlib import pyplot as plt 

def loadNeuralData(rats,days):
    """
    Iterate over list of rats and days
    For each, load the population/turns using RepCore.loadRecordingSesionData()
    with default params. Assign pop and turns to nested dict
    population[rat][day]['units'] and ...['turns'].
    Return population
    """
    population = {}
    for rat, day in zip(rats, days):
        if rat not in population.keys():
            population[rat] = {}
        if day not in population[rat].keys():
            population[rat][day] = {}
            
        pop,turns = RepCore.loadRecordingSessionData(rat,day)
        population[rat][day]['units'] = pop
        turns, refturns = filterVisitDf(turns)
        population[rat][day]['turns'] = turns
        population[rat][day]['refturns'] = refturns
        
    return population
        

def runRandomForest(X,Y,**kwargs):
    
    oobs = []
    if 'ntrees' in kwargs.keys():
        ntrees = kwargs['ntrees']
    else:
        ntrees = 1000
    if 'nreps' in kwargs.keys():
        nreps = kwargs['nreps']
    else:
        nreps = 1
        
    
    for i in range(nreps):
        clf = RandomForestClassifier(n_estimators=ntrees, oob_score=True)
        clf.fit(X,Y)
        oobs.append(clf.oob_score_)
        
    return oobs

def filterVisitDf(turns):
    """
    Take turn df and filter turns which involve 
    a turn around behavior
    Return: filtered turn df
    """
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
            
    refturns = copy.deepcopy(turns)
    turns = turns.iloc[ballisticTurnIdx]
    
    return turns, refturns


def createLabelArrayAlley(turns, refturns, targetlabel):
    """
    Create array of target labels of behavioral
    variables to decode. Labels based on alley+
    being current location.
    
    inputs:
        turns = filtered turn df (remove turnarounds)
        refturns = all turns df
    Options:
        current direction   = "CurrDir"
        previous direction  = "PrevDir"
        next direction      = "NextDir"
        turn in             = "TurnIn"
        turn out            = "TurnOut"
        current ego turn    = "CurrEgo"
        next ego turn       = "NextEgo"
        prev ego turn       = "PrevEgo"
        trajectory          = "Traj"
        
        Result will be (n,4) df with columns:
        col 0 - target label
        col 1 - current alley (for subsetting regions)
        col 2 - ts start of this behavioral event
        col 3 - ts end of '' ''
    """
    ballisticTurnIdx = [i for i,_ in refturns.iterrows()]
    labels, currAlleys, starts, stops = [], [], [], []
    for t, turn in turns.iterrows():
        
        turn_nm1 = refturns.iloc[t-1]         
        turn_np1 = refturns.iloc[t+1]
        
        # This gets us activty at Alley+ for each turn
        start, stop = turn['Ts entry'], turn_np1['Ts exit']
        
        currAlley = turn['Alley+']
        starts.append(start)
        stops.append(stop)
        currAlleys.append(currAlley)
        
        if targetlabel == 'CurrDir':
            labels.append(turn['Allo+'])
        elif targetlabel == 'PrevDir':
            labels.append(turn['Allo-'])
        elif targetlabel == 'NextDir':
            #check if next turn was ballistic bc
            if turn_np1.name in ballisticTurnIdx:
                labels.append(turn_np1['Allo+'])
        elif targetlabel == 'TurnIn':
            labels.append(f"{turn['Allo-']}{turn['Allo+']}")
        elif targetlabel == 'TurnOut':
            if turn_np1.name in ballisticTurnIdx:
                labels.append(f"{turn['Allo+']}{turn_np1['Allo+']}")
        elif targetlabel == 'CurrEgo':
            labels.append(turn['Ego'])
        elif targetlabel == 'NextEgo':
            if turn_np1.name in ballisticTurnIdx:
                labels.append(turn_np1['Ego'])
        elif targetlabel == 'PrevEgo':
            if turn_nm1.name in ballisticTurnIdx:
                labels.append(turn_nm1['Ego'])
        elif targetlabel == 'Traj':
            if turn_np1.name in ballisticTurnIdx:
                labels.append(f"{turn['Allo-']}{turn['Allo+']}{turn_np1['Allo+']}")
        
    labeldf = pd.DataFrame(data = list(zip(labels,
                                           currAlleys,
                                           starts,
                                           stops
                                           )),
                           columns = ["Label",
                                      "CurrAlley",
                                      "StartTs",
                                      "StopTs"
                                      ]
                            
                            )
    
    return labeldf
 

def createNeuralResponseArray(pop, labels, group='All'):
    """
    Input: pop - single recording day. dict where keys are neuron names in
                TTx\\cl-maze1.t format. keys are Unit() class instances
        labels - df (n,4) cols = label,curralley,ts_start,ts_end
        group - string "All", "Repeating", "Non-Repeating" indicating what kind of cells to include
    Create (n,t) array n= n neurons, t = t behavioral events. Each value
    is average response of neuron n to occurance of behavioral event t (avg meaning
    average over time of single behavioral event, not averaging across events)
    
    NB This does not filter by where a neuron's field is, it gets all neurons whose
    quality and overall activity have passed thresholds for the whole recording day
    
    return X array (n,t)
    """
    X = []
    for start, stop in zip(labels['StartTs'].astype(float), labels['StopTs'].astype(float)):
        pv  = []
        for _,unit in pop.items():
            if group == 'All':
                pv.append(unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]/((stop-start)/1e6))
            elif group == 'Repeating' and unit.repeating == True:
                pv.append(unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]/((stop-start)/1e6))
            elif group == 'Non-Repeating' and unit.repeating == False:
                pv.append(unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]/((stop-start)/1e6))

        X.append(pv)
    return np.asarray(X)
        

def calcDirBias(turns, regions):
    # for vert alleys E,W will be 0 and vice-versa
    dirbias = {i:{'N':0,'S':0,'E':0,'W':0} for i in regions}  
    for alley in dirbias.keys():
        for d in ['1', '2', '3', '4']:
            count = sum(turns[turns['Alley+']==alley]['Allo+']==d) # expression in sum() gives boolean series if direction was d or not. sum gives count thereof
            dirbias[alley][Def.allocodedict[d]] = count
        
    for alley in dirbias.keys():
        totalvisits = sum(dirbias[alley].values()) # total num passes in alley
        for d in ['N', 'S', 'E', 'W']:
            dirbias[alley][d] = dirbias[alley][d]/totalvisits # num visits in a dir / total visits, gets you prob
    
    return dirbias 

def naiveDirectionClassifier(turns, regions, **kwargs):
    dirbias = calcDirBias(turns, regions)
    if 'nreps' in kwargs.keys():
        nreps = kwargs['nreps']
    else:
        nreps = 1000
    perfs = []
    for n in range(nreps):
        outcomes = []
        for alley in dirbias.keys():
            alley = str(alley)
            turnsubset = turns[turns['Alley+']==str(alley)]
            for _,turn in turnsubset.iterrows():
                c=np.random.choice(['N','S','E','W'],size=1,p=[dirbias[alley]['N'],dirbias[alley]['S'],dirbias[alley]['E'],dirbias[alley]['W']])[0]
                real = Def.allocodedict[turn['Allo+']]
                if real == c:
                    outcomes.append(1)
                else:
                    outcomes.append(0)
        perfs.append(sum(outcomes)/len(outcomes))
        
    return perfs
        
    
    
    
        
def saveResults(**kwargs):
    """
    Accept a variable kwargs dict of information to save
    'data' must be one key with a 1d array or list of values
    'savepath' is path to dir to save things.
    'realPerf' is real score, usually 95th percentile mean
    'shuffPerf' is shuffled score, usually 95th percentile of shuffle
    'target' is target label
    """

    if 'data' not in kwargs.keys():
        print("Error - no field 'data' in input to saveResults")
    if 'savepath' not in kwargs.keys():
        print("Error - no save path given to saveResults")
        
    if os.path.isfile(kwargs['savepath']+"data.json"):
        with open(kwargs['savepath'],'r') as f:
            olddata = json.load(f)
            
        newdata = {}
        for k,v in olddata:
            newdata[k] = v
        for k,v in kwargs['data']:
            #json cant do np arrays
            if type(v)==np.ndarray:
                v = list(v)
            newdata[k] = v
            
        with open(kwargs['savepath']+"data.json", 'w') as f:
            json.dump(newdata, f)
    
    else:
        newdata = {}
        for k,v in kwargs['data']:
            #json cant do np arrays
            if type(v)==np.ndarray:
                v = list(v)
            newdata[k] = v
        with open(kwargs['savepath']+"data.json", 'w') as f:
            json.dump(newdata, f)
        
    title = f"{kwargs['rat']}{kwargs['day']} {kwargs['decoderType']} {kwargs['target']} {kwargs['regionlabel']}"

    fig, ax = plt.subplots(figsize=(10,10))
    ax.hist(kwargs['data'],bins=25)
    ax.set_xlabel("Performance")
    ax.set_ylabel("Frequency")
    if 'realPerf' in kwargs.keys():
        ax.vlines(0,100,kwargs['realPerf'],color='r')
        title = title + f"Real={round(kwargs[realPerf],2)} "
    if 'shuffPerf' in kwargs.keys():
        ax.vlines(0,100,kwargs['shuffPerf'],color='k')
        title = title + f"Shuff={round(kwargs[shuffPerf],2)} "
    ax.set_title(title,fontsize=10)
    plt.savefig(path=kwargs['savePath']+f"{kwargs['rat']}{kwargs['day']}_{kwargs['decoderType']}{kwargs['target']}_{kwargs['regionlabel']}.png",dpi=300)
    
        
        
  