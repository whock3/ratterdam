# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 16:09:37 2021

@author: whockei1

Script to iterate over days and rats to create a csv for R. For intersection data 

Difference between this and file with repetition_createSuperpopAlleyBehavioralModel
is here we get activity in intersection, in other file activiy in alley.
Related, we get current ego turn and in other we get current allo. also get
next/prev allo and ego vars in both. 
(doesnt make sense to get current allo dir
in inter bc youre often turing,which is 2 dirs together)

For each datset, get every pass through alleys overlapping each field (defined
as alley+ for convention). Get neural activity as intersection of activity on
that intersection and in field bound. Then label with behavioral variables. This is one sample. 
"""

#%% Imports
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
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import repeatingPC as repPC
import copy
import more_itertools
from matplotlib import path

def drawRegion(ax, bounds,color):
    """
    Bounds are the corners of a region on the track
    Format is [[xmin, xmax], [ymin, ymax]]
    Ax is an axis to which to add the region
    Color can be anything in practice it will be the rate 
    """
    x0,y0 = bounds[0][0], bounds[1][0]
    w,h = bounds[0][1] - bounds[0][0], bounds[1][1] - bounds[1][0]
    ax.add_patch(Rectangle((x0,y0),w,h,color=color))
    ax.autoscale_view() # for some reason the axes dont update automatically, so run this


def drawTrack(ax=None):
    if ax == None:
        ax = plt.gca()
    for i in range(17):
        drawRegion(ax, ratborders.alleyInterBounds[str(i)],'lightgrey')
              
#want each cell and field across all the days/rats we're analyzing to be coded with unique
# label. So lmer RE component doesnt think, e.g. field 4 from two different repeating
# cells is the same thing. 
cellcounter, fieldcounter = 0, 0
previousDirection , nextDirection, currentEgo, turnNums, rates = [], [], [], [], []
cellnames, fieldnums, cellIDs, fieldIDs = [], [], [], []
repeating = []
startTimes = []
rats, days, inters = [], [], []
nfields = [] # number of fields a cell has. want to look at things by fieldedness and # fields 
traversal = [] # list of bools corresponding to whether rat went thru alley (True) or turned around (False)
allocodedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
egocodedict = {'1':'S','2':'R','3':'B','4':'L','0':'X'}

#define alleys in groups to label visits accordingly later
verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]
trackperimeter = [str(i) for i in [0,4,6,7,9,10,13,15,16,2]]
trackinterior = [str(i) for i in [1,3,14,12,5,11, 8]]



for rat, day in zip(['R765','R781','R781','R808','R808','R859','R859','R886','R886'], ['RFD5','D3','D4','D6','D7','D1','D2','D1','D2']):
    df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    
    superpopInterDf = []
    
    ratborders = nab.loadAlleyBounds(rat, day)
    rewards = RepCore.readinRewards(rat, day)

    codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
    
    # Remove turnarounds/pivots
    ballisticTurnIdx = []
    for i in range(1,turns.shape[0]-1):
       row = turns.iloc[i]
       inter = row['Inter']
       # edit 10/2 removing check that last turn's inter wasnt the same,
       # i.e if alley- had a turnaround. since we are looking at things
       # in terms of alley+, only remove a turn if thats where a turnaround was
       if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter:
           ballisticTurnIdx.append(i)
    
    refturns = copy.deepcopy(turns) # keep a copy without filtering.
    turns = turns.iloc[np.asarray(ballisticTurnIdx)]
    
    
    for clustName, unit in population.items():
        repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
        
        for fnum, (field, perim, foverlap) in enumerate(zip(unit.fields, unit.perimeters, overlaps)):
        
            foverlap = [str(i) for i in foverlap]
            # Field size data, used to filter visits according to threshold pct'ages of field size
            maxx,minx, maxy, miny = max(perim[:,0]), min(perim[:,0]), max(perim[:,1]), min(perim[:,1])
            dista = np.sqrt(((minx-minx)**2)+((maxy-miny)**2))
            distb = np.sqrt(((minx-maxx)**2)+((miny-miny)**2))
            
            contour = path.Path(perim)
            
            # For alleys only - loop over using alley+ as convention for 'current' alley
            # check if the visit to the field is valid
            
            # R model will just use the IDs which will be unique for each cell and field.
            # but I still want the real name and field num there for inspection purposes
            for tnum, turn in refturns.iterrows():
                if tnum < refturns.shape[0]-1:
                
                    # Because turn definition is staggered as it moves across track,
                    # pick either alley+ or alley- by convention and make the ts, labels
                    # match that choice 
                    
                    if turn['Inter'] in foverlap: 
                        
                        ts_start, ts_end = float(turn['Ts entry']), float(refturns.iloc[tnum+1]['Ts exit'])
                        behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
                        behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
                        
                        filtOutcome = RepCore.filterVisit(dista,distb,behav,perim,
                                                          length_thresh=0.2,
                                                          dist_thresh=0.1,
                                                          dist_point_thresh=2,
                                                          inside_point_thresh=2)
                        
                        if filtOutcome == True:
                            rats.append(rat)
                            days.append(day)
                            turnNums.append(tnum)
                            inters.append(turn['Inter'])
                            
                            # 1-20-22: copied code from script to create superpop for alley data
                            # That script notes rewards so they can be removed later. They may confound response
                            # But you don't get rewarded at Intersections, so don't note them. may add back if I change mind.
                            
                            # isReward = np.where(np.asarray([(ts_start < i < ts_end) for i in rewards])==True)[0]
                            # if isReward.shape[0]>0:
                            #     visitRewards.append(True)
                            # else:
                            #     visitRewards.append(False)
                            
                            previousDirection.append(allocodedict[turn['Allo-']])
                            nextDirection.append(allocodedict[turn['Allo+']])
                            currentEgo.append(egocodedict[turn['Ego']])
                            repeating.append(repeat) # label if cell is repeating so we can group by repeating status later in R
                            startTimes.append(ts_start)
                            
                            if tnum in ballisticTurnIdx:
                                traversal.append(True)
                            else:
                                traversal.append(False)
                            
                            
                            # get spikes on alley (using the ts as endpoints) and filter by field perim to get spikes in field
                            rates.append(contour.contains_points(unit.spikes[(unit.spikes[:,0]>ts_start)&(unit.spikes[:,0]<= ts_end),1:]).shape[0]/((ts_end-ts_start)/1e6))
                            cellnames.append(f"{clustName}")
                            fieldnums.append(fnum)
                            cellIDs.append(cellcounter)
                            fieldIDs.append(fieldcounter)
                            nfields.append(len(unit.fields))
                        
                        
                        
            fieldcounter += 1
            
        cellcounter += 1
                

    
df = pd.DataFrame(data=list(zip(rats,
                                days,
                                cellnames,
                                cellIDs,
                                fieldnums,
                                fieldIDs,
                                nfields,
                                inters,
                                previousDirection,
                                currentEgo,
                                nextDirection,
                                traversal,
                                repeating,
                                startTimes,
                                rates
                                    )), 
columns=["Rat",
         "Day",
         "CellName",
         "CellID",
         "FieldNum",
         "FieldID",
         "NumFields",
         "Inters",
         "PrevDir",
         "CurrEgo",
         "NextDir",
         "Traversal",
         "Repeating",
         "StartTimes",
         "Rate"])

df.dropna(inplace=True)
# drop rows with missing data, coded as a '0' in the turn df (and then an X in the lookup table)
for var in ['PrevDir','NextDir', 'CurrEgo']:
    df.drop(df[df[var]=='X'].index, inplace=True)
    

stamp = util.genTimestamp()
filename = f"{stamp}_superPopInterBehaviorResponse_{Def.velocity_filter_thresh}vfilt.csv"
df.to_csv(f"E:\\Ratterdam\\R_data_repetition\\{filename}", header=True, index=False)
