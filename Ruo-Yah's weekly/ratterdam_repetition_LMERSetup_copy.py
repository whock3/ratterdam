# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 20:05:57 2021
@author: whockei1
"""

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
import math
import bisect
import pandas as pd
from alleyTransitions2 import alleyTransitions, closestTurnToVisit
from newAlleyBounds import R781, R808, R859, R886
from directionality import rateDir2D


# Define functions and parameters

tslice = 1e6*0.25 # time window before and after visit start used to make segment
nepochs = 3 # number of epochs in time to break the session into

def getPositionSlice(pos,begin,end):
    posslice =  pos[(pos[:,0]>=begin)&(pos[:,0]<end)]
    posslice = posslice[(posslice[:,1]!=0.0)&(posslice[:,2]!=0.0)]
    return posslice

def angle_trunc(a):
    while a < 0.0:
        a += math.pi * 2
    return a

def getAngleBetweenPoints(x_orig, y_orig, x_landmark, y_landmark):
    deltaY = y_landmark - y_orig
    deltaX = x_landmark - x_orig
    ang_rad = angle_trunc(math.atan2(deltaY, deltaX))
    return ang_rad*180/math.pi

def calcFieldEpochs(unit, field):
    maxtime = max([max(f[:,0]) for f in unit.fields]) # in nlts
    mintime = min([min(f[:,0]) for f in unit.fields]) # in nlts
    intervals  = np.linspace(mintime,maxtime,nepochs+1)    
    visitEpochs = [bisect.bisect_left(intervals,t) for t in field[:,0]]
    
    return visitEpochs

def calcFieldDirections(pos, field):
    dirs = []
    begins, ends = field[:,0] - tslice, field[:,0] + tslice
    for b,e in zip(begins, ends):
        
        ps = getPositionSlice(pos, b, e)
        ang = getAngleBetweenPoints(ps[0,1], ps[0,2], ps[-1,1], ps[-1,2])
        dirs.append(ang)
    
    return dirs
    


# Load data into dictionary

rat = "R886"
day = "D3"
#savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(df)
dfRep = "E:\\Ratterdam\\R859\\ratterdam_tabulations\\"
clustListRep = util.readRepeatingCells(f"{rat}{day}_tabulations.csv", dfRep)
population = {}
qualThresh = 0

clustListGoodQual = []
for i,clust in enumerate(clustList):
    if clustQuals[i] >= qualThresh:
        clustListGoodQual.append(clust)
clustListRep = list(set(clustListGoodQual)&set(clustListRep))
for i,clust in enumerate(clustListRep):
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
    
#%% 
# Create pandas data frame and save it 
rates, cells, fields, epochs, directions = [], [], [], [], []
prevTurnAllo, prevTurnEgo, nextTurnAllo, nextTurnEgo = [], [], [], []

if rat == "R781":
    ratAlleyBounds = R781
elif rat == "R808":
    ratAlleyBounds = R808
elif rat == "R859":
    ratAlleyBounds = R859
elif rat == "R886":
    ratAlleyBounds = R886

#for unitname, unit in population.items():
    #_, turns = alleyTransitions(unit.position, ratAlleyBounds)
    #break

count = [0,0,0]
totalCells = 0
oneField = 0
for unitname, unit in population.items():
    try:
        ANOVA = rateDir2D(unit)
        for i in range(3):
            if ANOVA["PR(>F)"][i] < 0.05:
                count[i] += 1
        totalCells += 1
    except:
        if len(unit.perimeters) <= 1:
            oneField += 1
        else:
            print(f"failed ANOVA: {unitname}")
print(f"{count}, rep: {len(clustListRep)}, totalwANOVA: {totalCells}, one field: {oneField}")
    #for i, field in enumerate(unit.fields):
        
        #rates.extend(field[:,1])
        #cells.extend([unitname]*field.shape[0])
        #fields.extend([i]*field.shape[0])
        #directions.extend(calcFieldDirections(unit.position, field))
        #epochs.extend(calcFieldEpochs(unit, field))
        
        #pta, pte, nta, nte = closestTurnToVisit(turns, field)
        #prevTurnAllo.extend(pta)
        #prevTurnEgo.extend(pte)
        #nextTurnAllo.extend(nta)
        #nextTurnEgo.extend(nte)
"""
data = {'rate':rates, 
        'cell':cells, 
        'field':fields,
        'direction':directions,
        'epoch':epochs,
        'prevTurnAllo':prevTurnAllo,
        'prevTurnEgo':prevTurnEgo,
        'nextTurnAllo':nextTurnAllo,
        'nextTurnEgo':nextTurnEgo
        }

alldata = pd.DataFrame(data=data)
alldata.dropna(inplace=True) 

#save data
stamp = util.genTimestamp()
filename = f"{stamp}_{rat}{day}_{Def.velocity_filter_thresh}vfilt_.csv"
               
alldata.to_csv(f"C:\\Users\\Ruo-Yah Lai\\Desktop\\My folder\\College\\Junior\\K lab research\\Graphs\\LMER\\{filename}", header=True, index=False)
"""