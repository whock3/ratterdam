# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 14:13:38 2021

@author: whockei1

Script to set up csv for LMER in R 
Repetition Project
Specifically looking at turns in ego/allo coordinates and effect on firing
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
import newAlleyBounds as nab
import math
import bisect
import pandas as pd

#%% Support functions




#%% Load Data


rat = "R859"
day = "D2"
ratBorders = nab.loadAlleyBounds(rat, day)
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(df)
population = {}
qualThresh = 3

for i,clust in enumerate(clustList):
    
    if clustQuals[i] >= qualThresh:
   
        try:
            print(clust)
            unit = RepCore.loadRepeatingUnit(df, clust, smoothing=1)
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
        except:
            pass
        
# position data is same for all units, so we just juse the last one processed here
pos, turns = alleyTrans.alleyTransitions(unit.position, ratBorders, graph=False)

    
    
#%%



# Create pandas data frame and save it 
rates, cells, fields, alloPrevPre, egoPrev, alloPrevPost, alloCurrPre, egoCurr, alloCurrPost, alloNextPre, egoNext, alloNextPost = [],[],[],[],[],[],[],[],[],[],[],[]

for unitname, unit in population.items():
    
    for i, field in enumerate(unit.fields):
        
        for v,visit in enumerate(field):
            
            # turns[:,3] is ts of entrance to intersection part of turn. good enough for now
            # to match up field visit and turn ts. code just finds arg that mins the diff btwen the two
            turnIdx = np.argmin(np.abs(turns[:,3].astype(np.double)-visit[0]))
            
            # first/last turn has no previous/next turn, so just ignore it
            if turnIdx > 0 and turnIdx < turns.shape[0]-1:
                
                rates.append(visit[1])
                cells.append(unitname)
                fields.append(i)
                
                turn = turns[turnIdx-1:turnIdx+2]
                alloPrevPre.append(str(turn[0,0]))
                egoPrev.append(str(turn[0,1]))
                alloPrevPost.append(str(turn[0,2]))
                
                alloCurrPre.append(str(turn[1,0]))
                egoCurr.append(str(turn[1,1]))
                alloCurrPost.append(str(turn[1,2]))
                
                alloNextPre.append(str(turn[2,0]))
                egoNext.append(str(turn[2,1]))
                alloNextPost.append(str(turn[2,2]))
            
                
        

data = {'rate':rates, 
        'cell':cells, 
        'field':fields, 
        'alloPrevPre':alloPrevPre,
        'egoPre':egoPrev,
        'alloPrevPost':alloPrevPost,
        'alloCurrPre':alloCurrPre,
        'egoCurr':egoCurr,
        'alloCurrPost':alloCurrPost,
        'alloNextPre':alloNextPre,
        'egoNext':egoNext,
        'alloNextPost':alloNextPost
        }

alldata = pd.DataFrame(data=data)
alldata.dropna(inplace=True) 

alldata.drop(alldata[alldata['alloPrevPre'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['alloPrevPost'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['egoPre'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['alloCurrPre'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['alloCurrPost'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['egoCurr'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['alloNextPre'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['alloNextPost'] == '0'].index, inplace = True)
alldata.drop(alldata[alldata['egoNext'] == '0'].index, inplace = True)


#save data
stamp = util.genTimestamp()
filename = f"{stamp}_{rat}{day}_TurnLMER_{Def.velocity_filter_thresh}vfilt_.csv"
               
alldata.to_csv(f"E:\\Ratterdam\\R_data_repetition\\{filename}", header=True, index=False)
