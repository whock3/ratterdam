# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 15:00:00 2022

@author: whockei1

Compare interior versus perimeter alleys or intersections. Comparisons relate
to directional or route related responses. Rationale for comparing alleys
and intersections is that they have different choice degrees of freedom, in theory
"""

import pandas as pd, numpy as np
from matplotlib import pyplot as plt 
import ratterdam_Defaults as Def 
import utility_fx as util
from collections import Counter


# Saving figures to: E:\Ratterdam\temp\PathTuningByLocation 


datapath = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

perimeterDf = df[df.Location=='P']
interiorDf = df[df.Location=='I']

#%% Directionality 

perimeterDiffs = []
interiorDiffs = []


for alleyTypeDf, alleyTypeDiffs  in zip([perimeterDf, interiorDf],[perimeterDiffs, interiorDiffs]):
    for fid, field in alleyTypeDf.groupby("FieldID"):
        #field can spill over multiple alleys within perimeter or interior alleysets
        oriens = np.unique(field.Orientation)
        for o in oriens:
            ofield = field[field.Orientation==o]
            dirs = np.unique(ofield.CurrDir)
            if len(dirs) > 2:
                print(f"ERROR - too many directions for {o} {fid}")
            diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
            alleyTypeDiffs.append(diff)
        
 
perimeterDiffs = np.asarray(perimeterDiffs)
interiorDiffs = np.asarray(interiorDiffs)

plt.hist(interiorDiffs,bins=50,density=True,label='Interior Alleys')
plt.hist(perimeterDiffs,bins=50,alpha=0.5,density=True,label='Perimeter Alleys')
plt.title("Signed Mean Directional Firing Rate Difference, Per Field", fontsize=22)
plt.ylabel('Density', fontsize=16)
plt.xlabel('Difference in FR, Direction A - B', fontsize=16)
plt.legend()



#%% SEM of samples by direction

perimeterDiffs = []
interiorDiffs = []


for alleyTypeDf, alleyTypeDiffs  in zip([perimeterDf, interiorDf],[perimeterDiffs, interiorDiffs]):
    for fid, field in alleyTypeDf.groupby("FieldID"):
        #field can spill over multiple alleys within perimeter or interior alleysets
        alleyTypeDiffs.append(field.Rate.std())
        
 
perimeterDiffs = np.asarray(perimeterDiffs)
interiorDiffs = np.asarray(interiorDiffs)

plt.hist(interiorDiffs,bins=25,density=True,label='Interior Alleys')
plt.hist(perimeterDiffs,bins=25,alpha=0.5,density=True,label='Perimeter Alleys')
plt.title("Average Std, Per Field. Orientation-Filtered if necessary", fontsize=22)
plt.ylabel('Density', fontsize=16)
plt.xlabel('Average Std', fontsize=16)
plt.legend()
