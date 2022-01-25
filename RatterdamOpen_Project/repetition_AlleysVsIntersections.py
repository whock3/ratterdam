# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 16:35:21 2022

@author: whockei1

Looking at directionality/path coding in alleys vs intersections
"""

import pandas as pd, numpy as np
from matplotlib import pyplot as plt 
import ratterdam_Defaults as Def 
import utility_fx as util
from collections import Counter


# Saving figures to: E:\Ratterdam\temp\PathTuningByLocation 


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)


interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
interdf = pd.read_csv(interdatapath)

#%% Directionality 

alleyDiffs = []
interDiffs = []


# alleys
for fid, field in alleydf.groupby("FieldID"):
    #field can spill over multiple alleys within perimeter or interior alleysets
    oriens = np.unique(field.Orientation)
    for o in oriens:
        ofield = field[field.Orientation==o]
        dirs = np.unique(ofield.CurrDir)
        if len(dirs) > 2:
            print(f"ERROR - too many directions for {o} {fid}")
        diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
        alleyDiffs.append(diff)
        
# intersections 
for fid, field in interdf.groupby("FieldID"):
    dirs = np.unique(field.CurrEgo)
    meanDirDiffs = []
    for direction in dirs:
        meanDirDiffs.append(field[field.CurrEgo==direction].Rate.mean())
    interDiffs.append(max(meanDirDiffs)-min(meanDirDiffs))
        
 
alleyDiffs = np.asarray(alleyDiffs)
interDiffs = np.asarray(interDiffs)

plt.hist(alleyDiffs,bins=50,density=True,label='Alleys')
plt.hist(interDiffs,bins=50,alpha=0.5,density=True,label='Intersections')
plt.title(" Mean Directional Firing Rate Difference, Per Field", fontsize=22)
plt.ylabel('Density', fontsize=16)
plt.xlabel('Difference in FR, Direction A - B', fontsize=16)
plt.legend()




#%%


for iName, inter in interdf.groupby("Inters"):
    meanDiffs = []
    for dName, egoD in inter.groupby("CurrEgo"):
        meanDiffs.append(egoD.Rate.mean())
    print(iName)
    print(max(meanDiffs)-min(meanDiffs))
    
    
for aName, alley in alleydf.groupby("Alleys"):
    meanDiffs = []
    for dName, d in alley.groupby("CurrDir"):
        meanDiffs.append(d.Rate.mean())
    print(aName)
    print(max(meanDiffs)-min(meanDiffs))