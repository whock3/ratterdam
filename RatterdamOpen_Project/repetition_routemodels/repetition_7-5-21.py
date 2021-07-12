# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:21:08 2021

@author: whockei1

Goal: Create dataframe with turns involving a given field divided according
to the direction of traversal through that field.

"""

#%% Imports
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
import repeatingPC as repPC
import math
import bisect
import pandas as pd

#%% Load data
rat = "R859"
day = "D2"
ratBorders = alleyTrans.R859
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(datapath)
population = {}
qualThresh = 3

for i,clust in enumerate(clustList):
    
    if clustQuals[i] >= qualThresh:
   
        try:
            print(clust)
            unit = RepCore.loadRepeatingUnit(datapath, clust, smoothing=1)
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
        except:
            pass
        
        
# Session endpoints data
with open(datapath+"sessionEpochInfo.txt","r") as f:
    lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    
nepoch=3
intervals = np.linspace(start,end,nepoch+1)
        
#%% Create turn df

rattrans = alleyTrans.R859
pos, turns = alleyTrans.alleyTransitions(unit.position, rattrans, graph=False)
turns = pd.DataFrame(turns)
turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']

turns = pd.DataFrame(data=turns)
turns.dropna(inplace=True) 

#%% Functions to support above

def addTrajVariables(turnIdx, turns, turnarm):
    """
    Input turnIdx - index of the turn in the turns df
          turns - df of the turn data for each turn
          turnarm - str ("pre", "post") indicating whether the field
                      is on the pre-alley or the post-alley of the turn
                      (for an intersection field the value should be "pre")
          
    Fx pulls out the turn-related regressors for the model.
    Returns: list of variables
    """
    turn = turns.iloc[turnIdx]
    m1turn = turns.iloc[turnIdx-1] # n-1 turn
    p1turn = turns.iloc[turnIdx+1] # n+1 turn
    
    if turnarm == 'pre':
        bearing = turn['Allo-']
        m1bearing = m1turn['Allo-'] # n-1 bearing, n.b. this is different from m1turn. 
        p1bearing = turn['Allo+']
    elif turnarm == 'post':
        bearing = turn['Allo+']
        m1bearing = turn['Allo-']
        p1bearing = p1turn['Allo+']
        
    return [bearing, m1bearing, p1bearing]

#%%
dirdata = []
for unitname, unit in population.items():
    print(unitname)
    repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,alleyTrans.R781)
    
    for fi, field in enumerate(unit.fields):
        print(fi)
        foverlap = overlaps[fi]
        if len(foverlap) > 0:
            
            alleyoverlap = [i for i in foverlap if type(i)==int] #alleys are numbered, intersections are lettered by convention
            interoverlap = [i for i in foverlap if type(i)==str] # '' 
            
            if len(alleyoverlap) > 1:
                print("CASE NOT BEING ANALYZED RIGHT NOW - field in multiple alleys. Consider later.")
            else:
                # fix this - this is nonexhaustive and redundant at the same time
                if len(alleyoverlap) == 0 and len(interoverlap) ==1:
                    overlaptype = 'intersection'
                elif len(alleyoverlap)==1 and len(interoverlap) ==1:
                    overlaptype = 'alleyint'
                elif len(alleyoverlap) == 1 and len(interoverlap) == 0:
                    overlaptype = 'alley'
                    
                for visit in field:

                    turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
                    turndata = [np.nan]
                    if turnIdx > 0 and turnIdx < turns.shape[0]-1:
                        turn = turns.iloc[turnIdx]
                        
                        if overlaptype == 'alley':
                            if str(alleyoverlap[0]) == turn['Alley-']:
                                turndata = addTrajVariables(turnIdx, turns, "pre")
                            elif str(alleyoverlap[0]) == turn['Alley+']:
                                turndata = addTrajVariables(turnIdx, turns, "post")
                            else:
                                print("Turn ignored")
                        
                        elif overlaptype == 'intersection':
                            if turn['Ego'] == 1: # means forward through intersection which is all we're going to look at for now 7/5/21
                               turndata = addTrajVariables(turnIdx, turns, "pre") 
                        
                    
                        elif overlaptype == 'alleyint':
                            #print(turn['Alley-'], turn['Inter'], turn['Alley+'], turn['Ego'])
                            if str(alleyoverlap[0]) == turn['Alley-']:
                                turndata = addTrajVariables(turnIdx, turns, "pre")
                            elif str(alleyoverlap[0]) == turn['Alley+']:
                                turndata = addTrajVariables(turnIdx, turns, "post")
                            elif str(interoverlap[0]) == turn['Inter'] and turn['Ego'] == 1:
                                turndata = addTrajVariables(turnIdx, turns, "pre")
                            else:
                                print("Turn ignored") 
            
                    epoch = bisect.bisect_left(intervals, visit[0]) # returns what period of the session the visit was in. bisect returns insertion index to maintain order.
                    dirdata.append((unit.name, fi, visit[1], *turndata, epoch))
                    

df = pd.DataFrame(data=dirdata, columns=['unit','field','rate','dirC','dirM1', 'dirP1', 'epoch'])
df.dropna(inplace=True)
stamp = util.genTimestamp()
filename = f"{stamp}_{rat}{day}_dirLMER_{Def.velocity_filter_thresh}vfilt.csv"
df.to_csv(f"E:\\Ratterdam\\R_data_repetition\\{filename}", header=True, index=False)



