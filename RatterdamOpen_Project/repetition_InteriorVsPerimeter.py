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


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)


interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
interdf = pd.read_csv(interdatapath)

alleyperimeterDf = alleydf[alleydf.Location=='P']
alleyinteriorDf = alleydf[alleydf.Location=='I']

#repeating is true or false. If false it means single fields only.
# not looking at multiple field nonrepeating now (1-27-22) bc it's messy
repeating = False

#%% Directionality 

perimeterAlleyDiffs = []
interiorAlleyDiffs = []


for alleyTypeDf, alleyTypeDiffs  in zip([alleyperimeterDf, alleyinteriorDf],[perimeterAlleyDiffs, interiorAlleyDiffs]):
    for fid, field in alleyTypeDf.groupby("FieldID"):
        if repeating == True:
            if np.unique(field.Repeating)[0] == repeating:
                #field can spill over multiple alleys within perimeter or interior alleysets
                oriens = np.unique(field.Orientation)
                for o in oriens:
                    ofield = field[field.Orientation==o]
                    dirs = np.unique(ofield.CurrDir)
                    if len(dirs) > 2:
                        print(f"ERROR - too many directions for {o} {fid}")
                    diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                    alleyTypeDiffs.append(diff)
        elif repeating == False:
            if np.unique(field.NumFields)[0] == 1:
                #field can spill over multiple alleys within perimeter or interior alleysets
                oriens = np.unique(field.Orientation)
                for o in oriens:
                    ofield = field[field.Orientation==o]
                    dirs = np.unique(ofield.CurrDir)
                    if len(dirs) > 2:
                        print(f"ERROR - too many directions for {o} {fid}")
                    diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                    alleyTypeDiffs.append(diff)
            
        
 
perimeterAlleyDiffs = np.asarray(perimeterAlleyDiffs)
interiorAlleyDiffs = np.asarray(interiorAlleyDiffs)

plt.figure()
plt.hist(interiorAlleyDiffs,bins=25,density=True,label='Interior Alleys')
plt.hist(perimeterAlleyDiffs,bins=25,alpha=0.5,density=True,label='Perimeter Alleys')
plt.title(f"Unigned Mean Directional Firing Rate Difference, Per Field, repeating = {repeating}", fontsize=22)
plt.ylabel('Density', fontsize=16)
plt.xlabel('Difference in FR, Direction A - B', fontsize=16)
plt.legend()



# #%% SEM of samples by direction

# perimeterDiffs = []
# interiorDiffs = []


# for alleyTypeDf, alleyTypeDiffs  in zip([perimeterDf, interiorDf],[perimeterDiffs, interiorDiffs]):
#     for fid, field in alleyTypeDf.groupby("FieldID"):
#         #field can spill over multiple alleys within perimeter or interior alleysets
#         alleyTypeDiffs.append(field.Rate.std())
        
 
# perimeterDiffs = np.asarray(perimeterDiffs)
# interiorDiffs = np.asarray(interiorDiffs)

# plt.hist(interiorDiffs,bins=25,density=True,label='Interior Alleys')
# plt.hist(perimeterDiffs,bins=25,alpha=0.5,density=True,label='Perimeter Alleys')
# plt.title("Average Std, Per Field. Orientation-Filtered if necessary", fontsize=22)
# plt.ylabel('Density', fontsize=16)
# plt.xlabel('Average Std', fontsize=16)
# plt.legend()



#%% Interior vs perimeter - Intersections

interiorInters = ['F', 'G'] # every other intersection is perimeter

interperimeterDiffs = []
interinteriorDiffs = []

for fid, field in interdf.groupby("FieldID"):
    if repeating == True:
        if np.unique(field.Repeating)[0] == repeating:
        # most fields will overlap only 1 intersection (if any), but some are long enough to overlap 2
            inters = np.unique(field.Inters)
            for inter in inters:
                ifield = field[field.Inters==inter]
                dirs = np.unique(ifield.CurrEgo)
                meanDirDiffs = []
                for direction in dirs:
                    meanDirDiffs.append(ifield[ifield.CurrEgo==direction].Rate.mean())
                if inter in interiorInters:
                    interinteriorDiffs.append(max(meanDirDiffs)-min(meanDirDiffs))
                else:
                    interperimeterDiffs.append(max(meanDirDiffs)-min(meanDirDiffs))
    elif repeating == False:
        if np.unique(field.NumFields)[0] == 1:
        # most fields will overlap only 1 intersection (if any), but some are long enough to overlap 2
            inters = np.unique(field.Inters)
            for inter in inters:
                ifield = field[field.Inters==inter]
                dirs = np.unique(ifield.CurrEgo)
                meanDirDiffs = []
                for direction in dirs:
                    meanDirDiffs.append(ifield[ifield.CurrEgo==direction].Rate.mean())
                if inter in interiorInters:
                    interinteriorDiffs.append(max(meanDirDiffs)-min(meanDirDiffs))
                else:
                    interperimeterDiffs.append(max(meanDirDiffs)-min(meanDirDiffs))        
                
plt.figure()
plt.hist(interinteriorDiffs,bins=25,density=True,label='Interior Intersections')
plt.hist(interperimeterDiffs,bins=25,alpha=0.5,density=True,label='Perimeter Intersections')
plt.title(f"Average Mean Diff, Per Field, repeating = {repeating}", fontsize=22)
plt.ylabel('Density', fontsize=16)
plt.xlabel('Average Diff', fontsize=16)
plt.legend()



plt.figure()
plt.violinplot([perimeterAlleyDiffs, interiorAlleyDiffs, interperimeterDiffs, interinteriorDiffs])
plt.xticks(ticks=[1,2,3,4],labels=["Perimeter Alleys", "Interior Alleys", "Perimeter Inters.", "Interior Inters"], fontsize=16)
plt.title(f"Unsigned Mean Difference in Directional Firing At Different Maze Locations repeating = {repeating}",fontsize=25)
