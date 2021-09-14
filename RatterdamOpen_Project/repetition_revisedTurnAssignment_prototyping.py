# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 11:58:33 2021

@author: whockei1

Revised way of finding turns
Assign any turn that starts or ends on region field overlaps 
So most (all?) fields will have multiple turns associated w each visit
And spikes/occ samples from those segments of time will be used separately 
to define activity for each turn associated with visit to that field 
(Current way of doing it gets closest turn to 1st pt of field visit)

Will also exclude turnarounds/pivoting optionally, depending on analysis 
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
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import repeatingPC as repPC
import copy
import more_itertools


#%% Helper fx 

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
              
                    
        
        
#%% Load turns 
rat, day = 'R859', 'D1'
df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
turns,_ = RepCore.loadTurns(rat, day)
ratborders = nab.loadAlleyBounds(rat, day)
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}

# Remove turnarounds/pivots
ballisticTurnIdx = []
for i in range(1,turns.shape[0]-1):
   row = turns.iloc[i]
   inter = row['Inter']
   if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter and turns.iloc[i-1].Inter != inter:
       ballisticTurnIdx.append(i)

refturns = copy.deepcopy(turns) # keep a copy without filtering.
turns = turns.iloc[np.asarray(ballisticTurnIdx)]


#%% Load unit and associated data
clustName = 'TT7\\cl-maze1.3'
unit = RepCore.loadRepeatingUnit(df, clustName)
repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
#%% Load Field data
fnum = 1
field = unit.fields[fnum]
perim = unit.perimeters[fnum]
foverlap = overlaps[fnum]
foverlap = [str(i) for i in foverlap]
# Field size data, used to filter visits according to threshold pct'ages of field size
maxx,minx, maxy, miny = max(perim[:,0]), min(perim[:,0]), max(perim[:,1]), min(perim[:,1])
dista = np.sqrt(((minx-minx)**2)+((maxy-miny)**2))
distb = np.sqrt(((minx-maxx)**2)+((miny-miny)**2))

#%% Iterate over turns and assign to fields, get rates and labels

trajectories, directionsIn, directionsOut, turnNums = [], [], [], []

for tnum, turn in turns.iterrows():
    if turn['Alley-'] in foverlap or turn['Alley+'] in foverlap or turn['Inter'] in foverlap:
        
        # Get behavior and spikes using refturns which has all the turns bc we need the flanking turns
        # not because we are looking at nonlocal data (like for trajectory decoding) but because the 
        # turn df only has ts associated with the middle of the turn basically 
        ts_start, ts_end = float(refturns.iloc[tnum-1]['Ts entry']), float(refturns.iloc[tnum+1]['Ts exit'])
        behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
        behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
          
        filtOutcome = RepCore.filterVisit(dista,distb,behav,perim,length_thresh=0.3,dist_thresh=0.1,dist_point_thresh=3,inside_point_thresh=3)
        if filtOutcome == True:
            turnNums.append(tnum)
            trajectories.append(behav)
            directionsIn.append(turn['Allo-'])
            directionsOut.append(turn['Allo+'])
        
trajectories = np.asarray(trajectories)
directionsIn = np.asarray(directionsIn)
directionsOut = np.asarray(directionsOut)
turnNums = np.asarray(turnNums)

#%% Plot all data 

ncols = 12
nturns = trajectories.shape[0] # all the above arrays should have same shape
fig, ax = plt.subplots(int(np.ceil(nturns/ncols)),ncols,figsize=(15,8))

for i, (traj, dIn, dOut, tnum) in enumerate(zip(trajectories, directionsIn, directionsOut, turnNums)):
    cax = fig.axes[i]
    drawTrack(ax=cax)
    
    cax.plot(perim[:,0], perim[:,1],c='r')
    cax.plot(traj[:,1], traj[:,2],c='k',zorder=99)
    cax.scatter(traj[0,1], traj[0,2], c='g', marker='o',s=20,zorder=99)
    cax.scatter(traj[-1,1], traj[-1,2], c='r', marker='o',s=20,zorder=99)
    cax.set_title(f"{i},{tnum}, {codedict[dIn]}->{codedict[dOut]}",fontsize=10)
    
for ii in fig.axes:
    ii.axis("off")


#%% Plot individual turns, mostly for diagnostic 

plt.figure()
drawTrack()
tnum = 417
turn = refturns.iloc[tnum]
ts_start, ts_end = float(refturns.iloc[tnum-1]['Ts entry']), float(refturns.iloc[tnum+1]['Ts exit'])
behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
plt.plot(perim[:,0], perim[:,1],c='r',zorder=99)
plt.plot(behav[:,1], behav[:,2], c='k')
plt.scatter(behav[0,1],behav[0,2],c='g',marker='o',s=20,zorder=99)
plt.scatter(behav[-1,1],behav[-1,2],c='r',marker='o',s=20,zorder=99)
plt.title(f"{unit.name} Field {fnum} Turn {turn.name}")

    
#%% Goal here is to take all turns involving a field and group them by
# the axis of travel.

# group turns into sets of consecutive turns, each group being the
# turns that involve the field for a given visit to that gield
visitTurns = [list(group) for group in more_itertools.consecutive_groups(turnNums)]

#%% Create directional segments for each place field traversal
# get firing rate during that period and add to dataframe 

directionSegmentDf = []

for group in visitTurns:

    visitData = []
            
    # get start, end timestamps for each alley along the traversal through the field
    # The alleys must overlap w the field, i.e we are not getting the alleys pre/post traversal
    for idx in group:
        turn = refturns.iloc[idx]
        if turn['Alley-'] in foverlap:
            if len(visitData)>=1 and refturns.iloc[idx-1]['Ts entry'] not in visitData[-1]:
                visitData.append((refturns.iloc[idx-1]['Ts entry'], turn['Ts exit'], codedict[turn['Allo-']]))
            elif visitData == []:
                visitData.append((refturns.iloc[idx-1]['Ts entry'], turn['Ts exit'], codedict[turn['Allo-']]))
        if turn['Alley+'] in foverlap:
            if len(visitData)>=1 and turn['Ts entry'] not in visitData[-1]:
                visitData.append((turn['Ts entry'], refturns.iloc[idx+1]['Ts exit'], codedict[turn['Allo+']]))
            elif visitData == []:
                visitData.append((turn['Ts entry'], refturns.iloc[idx+1]['Ts exit'], codedict[turn['Allo+']]))
                
        # check passes through inter when field doesnt overlap involved alleys. Wont have double counting here. 
        if turn['Alley-'] not in foverlap and turn['Alley+'] not in foverlap and turn['Inter'] in foverlap:
            # intersection traversals with a 90deg turn are ambiguous as to their direction without more involved analysis. 
            if turn['Allo-'] == turn['Allo+']:
                visitData.append((refturns.iloc[idx-1]['Ts entry'], refturns.iloc[idx+1]['Ts exit'], codedict[turn['Allo-']]))
                
    
    
    
    # Now group data by common direction, i.e. if there are two consecutive alleys
    # traversed going north, then get the starting ts of the first and the ending
    # ts of the second and call that the "north leg" of the traversal. NB this should
    # work for arbitrary # of legs in the same direction and repeated legs in the same direction
    # e.g. N,E,N for at least one alley length each
    segmentedHeadings = []
    heading = None
    for i,seg in enumerate(visitData):
        
        if heading == None and i < len(visitData)-1:
            segmentedHeadings.append([])
            heading = seg[2]
            startSegTs = seg[0]
            segmentedHeadings[-1].append(heading)
            segmentedHeadings[-1].append(startSegTs)
        
        
        if heading == None and i == len(visitData)-1:
            segmentedHeadings.append([])
            heading = seg[2]
            startSegTs = seg[0]
            endSegTs = seg[1]
            segmentedHeadings[-1].append(heading)
            segmentedHeadings[-1].append(startSegTs)
            segmentedHeadings[-1].append(endSegTs)
            
        elif heading == seg[2] and i == len(visitData)-1:
            endSegTs = visitData[i][1]
            segmentedHeadings[-1].append(endSegTs)
            
        elif heading != seg[2] and i == len(visitData)-1:
            endSegTs = visitData[i-1][1]
            segmentedHeadings[-1].append(endSegTs)
            heading = seg[2]
            startSegTs = seg[0]
            endSegTs = seg[1]
            segmentedHeadings.append([])
            segmentedHeadings[-1].append(heading)
            segmentedHeadings[-1].append(startSegTs)
            segmentedHeadings[-1].append(endSegTs)
            
        elif heading != seg[2] and i < len(visitData)-1:
            endSegTs = visitData[i-1][1]
            segmentedHeadings[-1].extend(endSegTs)
            heading = seg[2]
            startSegTs = seg[0]
            segmentedHeadings.append([])
            segmentedHeadings[-1].append(heading)
            segmentedHeadings[-1].append(startSegTs)
         
            
    #For each directional segment of a place field traversal, compute firing rate
    # note that the timestamps are from alley entrances/exits. And although the alleys in
    # question, by definition, overlap with the place field it doesnt mean the field
    # extends all the way to the alley edge. meaning there is some length of trajectory
    # and associated spiking data thats on the overlapping alley but not in the field
    # as defined by the session-averaged border. I dont think this is a huge deal
    # bc there shouldnt be much spiking and the session avg field bound doesnt necessarily
    # map onto where the field stops on a given traversal - WH 9/2/21 
    
    for segment in segmentedHeadings:
    
        direction = segment[0]
        ts_start, ts_end = segment[1], segment[2]
        spikes = unit.spikes[(unit.spikes[:,0]>ts_start)&(unit.spikes[:,0]<=ts_end)]
        rate = spikes.shape[0]/((ts_end-ts_start)/1e6)
        
        directionSegmentDf.append([unit.name, fnum, direction, ts_start])
        
    
df = pd.DataFrame(data=directionSegmentDf, columns=["Unit","Field","Direction","StartTime"])
df.dropna(inplace=True)

stamp = util.genTimestamp()
filename = f"{stamp}_{rat}{day}_segmentedDirections_{Def.velocity_filter_thresh}vfilt.csv"
df.to_csv(f"E:\\Ratterdam\\R_data_repetition\\{filename}", header=True, index=False)
