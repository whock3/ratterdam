# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 23:33:13 2021

@author: whockei1

Scrap file
"""
#%% Functions to load
def checkCheating(unit, turn, turn_np1):
    """
    deprecated. solving problem a different way
    
    Check to see if any of the unit's fields overlap with either of the non-local
    regions involved in this turn (i.e. the pre and post alleys were decoding dir
     from) and only if not do we let this unit contribute data to this turn
    """
    threshold = 0.10
    anyExcessiveOverlap = False
    for perim in unit.perimeters:
        for a in [turn['Alley-'], turn_np1['Alley+']]:
            alley = ratborders.alleyInterBounds[str(a)]
            xmin,xmax = alley[0][0],alley[0][1]
            ymin,ymax = alley[1][0],alley[1][1]
            alleyPerim = np.array([[xmin, ymin], [xmax, ymin],
                                           [xmax, ymax], [xmin, ymax]])
            alleySize = (xmax-xmin) * (ymax-ymin)
            overlap1 = repPC.overlap(perim, alleyPerim)
            if overlap1 > alleySize*threshold:
                anyExcessiveOverlap = True
            
    return anyExcessiveOverlap



#%% Turn based rate maps - setup
# Group data by allo dir pre (4 plots) and allo dir post (4 plots) +/- 1.5 turns
#(the half turn is bc the turns end at the intersection so punch out a bit more
# in time to get the full +/1 1 turn)
#Reminder about code order: N,E,S,W,  F,R,B,L


# turnRMS = {'Pre':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
#                   }, 
           
#            'Post':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
#                   }
#            }

# for field in unit.fields:
#     for visit in field:
#         turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
#         try:
#             # first/last turn has no previous/next turn, so just ignore it
#             if turnIdx > 2 and turnIdx < turns.shape[0]-2:
#                 spikesPre = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
#                 spikesPost = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

#                 turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'],spikesPre))
#                 turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'],spikesPost))
    
#                 occPre = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
#                 occPost = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

#                 turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'],occPre))
#                 turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'],occPost))
#         except:
#             print("Likely invalid code found")

# for d in ['Pre', 'Post']:
#     fig, ax = plt.subplots(2,2)
#     for i,code in enumerate([('1','N'),('2','E'),('3','S'),('4','W')]):
#         s,o = turnRMS[d][code[0]]['Spikes'], turnRMS[d][code[0]]['Pos']
#         rm = util.makeRM(s,o)
#         fig.axes[i].imshow(rm,aspect='auto',interpolation='None',cmap=cmap,vmax=7,origin='lower')
#         fig.axes[i].set_title(code[1])
#     plt.suptitle(f"{d}-turn Bearing")
    


#%% Working on turn visualization, getting rid of turn arounds, how to bridge gaps, etc
# 8-25-21

nturns = 50
nrows = 5
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}

fig, ax = plt.subplots(nrows,int(nturns/nrows),figsize=(12,8))

start = 100
for t in range(start, start+nturns):
    turn = turns.iloc[t]
    i = t-start+1
    ts_start, ts_end = refturns.iloc[turn.name-1]['Ts entry'], refturns.iloc[turn.name+1]['Ts exit']
    behav = unit.position[(unit.position[:,0]>float(ts_start))&(unit.position[:,0]<=float(ts_end))]
    behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
    
    for r in range(17):
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[str(r)],'lightgrey')
    for rr in ascii_uppercase[:12]:
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[rr],'lightgrey')
        
    for region in ['Alley-','Inter','Alley+']:
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[str(turn[region])],'lightblue')
        
    fig.axes[i-1].plot(behav[:,1], behav[:,2], color='k', zorder=99)
    fig.axes[i-1].scatter(behav[0,1],behav[0,2],color='green',marker='o',s=50,zorder=99)
    fig.axes[i-1].scatter(behav[-1,1],behav[-1,2],color='red',marker='o',s=50,zorder=99)
    
    label = '->'.join([codedict[i] for i in [turn['Allo-'],turn['Allo+']]])
    fig.axes[i-1].set_title(f"Turn {turn.name}, {label}")
    
for ii in range(len(fig.axes)):
    fig.axes[ii].axis("off")
    #fig.axes[ii].set_aspect('equal', adjustable='box')
plt.suptitle(f"{rat} {day} Turns {turns.iloc[start].name} - {turns.iloc[t].name} (Turn arounds removed)")

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
          

codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
                    
        
    
#%% Imports and load data for
import numpy as np
import utility_fx as util
import os
from matplotlib import pyplot as plt
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import newAlleyBounds as nab
import pandas as pd
from matplotlib.patches import Rectangle
import repeatingPC as repPC
import copy
from string import ascii_uppercase
import datetime

rat, day = 'R781', 'D3'
ratborders = nab.loadAlleyBounds(rat, day)
turns, unit = RepCore.loadTurns(rat, day)

#%% Filter turns 
# Filter to remove turn-around trajectories. These cause same label
# to map onto different behaviors
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


ballisticTurnRows = [turns.iloc[i].name for i in range(turns.shape[0])]

#%%
sessionStart = float(refturns.iloc[0]['Ts exit']) # technically missing time of first leg of first turn, but we dont analyze this turn anyway
nturns = 50
ncols = 10
fig, ax = plt.subplots(int(np.ceil(nturns/ncols)),ncols, figsize=(12,8))
start = 1
for t in range(start, start+nturns):
    turn = refturns.iloc[t]
    i = t-start+1
    
    turnTime = ((float(turn['Ts exit']) - sessionStart) / 1e6)
    ts_start, ts_end = refturns.iloc[t-1]['Ts entry'], refturns.iloc[t+1]['Ts exit']
    behav = unit.position[(unit.position[:,0]>float(ts_start))&(unit.position[:,0]<=float(ts_end))]
    behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
    spikes = unit.spikes[(unit.spikes[:,0]>float(ts_start))&(unit.spikes[:,0]<=float(ts_end))]

    
    if turn.name in ballisticTurnRows:
        bgcolor = 'lightcoral'
        turncolor = 'red'
    else:
        bgcolor = 'cornflowerblue'
        turncolor = 'blue'
    
    for r in range(17):
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[str(r)], bgcolor)
    for rr in ascii_uppercase[:12]:
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[rr], bgcolor)
        
    for region in ['Alley-','Inter','Alley+']:
        drawRegion(fig.axes[i-1],ratborders.alleyInterBounds[str(turn[region])], turncolor)
        
    fig.axes[i-1].plot(behav[:,1], behav[:,2], color='k', zorder=99)
    fig.axes[i-1].scatter(behav[0,1],behav[0,2],color='green',edgecolor='k',marker='o',s=75,zorder=99)
    fig.axes[i-1].scatter(behav[-1,1],behav[-1,2],color='red',edgecolor='k',marker='o',s=75,zorder=99)
    fig.axes[i-1].scatter(spikes[:,1], spikes[:,2], color='y', marker='o', s=75, zorder=99)
    
    
    label = '->'.join([codedict[i] for i in [turn['Allo-'],turn['Allo+']]])
    fig.axes[i-1].set_title(f"{turn.name}, {str(datetime.timedelta(seconds=round(turnTime)))}, {label}")
    
for ii in range(len(fig.axes)):
    fig.axes[ii].axis("off")
    #fig.axes[ii].set_aspect('equal', adjustable='box')
plt.suptitle(f"{rat} {day} {unit.name} Turns {refturns.iloc[start].name} - {refturns.iloc[t].name} (Reds=Ballistic, Blues=Pivots)")



#%% Visualize turns nearest to field visit and spikes overlaid (assumes data is loaded)
unit = population['TT9\\cl-maze1.3']
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
field = unit.fields[0]
perim = unit.perimeters[0]
turnIdx = []
allTurnIdx =[]
for visit in field:
    turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
    allTurnIdx.append(turnIdx)

allTurnIdx = np.asarray(allTurnIdx)
ncols=10
fig, ax = plt.subplots(int(np.ceil(allTurnIdx.shape[0]/ncols)),ncols,figsize=(12,8))
for i,tidx in enumerate(allTurnIdx):
    turn = turns.iloc[tidx]
    ts_start, ts_end = float(turns.iloc[tidx-1]['Ts entry']), float(turns.iloc[tidx+1]['Ts exit'])
    behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
    behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
    spikes = unit.spikes[(unit.spikes[:,0]>ts_start)&(unit.spikes[:,0]<=ts_end)]
    rate = spikes.shape[0]/((ts_end-ts_start)/1e6)
    drawTrack(ax=fig.axes[i])
    fig.axes[i].plot(behav[:,1], behav[:,2],color='k',zorder=99)
    fig.axes[i].plot(perim[:,0], perim[:,1],color='red',zorder=99)
    fig.axes[i].scatter(behav[0,1],behav[0,2],c='g',marker='o',s=20,zorder=99)
    fig.axes[i].scatter(behav[-1,1],behav[-1,2],c='r',marker='o',s=20,zorder=99)
    fig.axes[i].scatter(spikes[:,1],spikes[:,2],c='blue',marker='o',s=20,zorder=99)
    label = '->'.join([codedict[j] for j in [turn['Allo-'],turn['Allo+']]])
    fig.axes[i].set_title(f"{i},{turn.name}, {label}")
for ii in range(len(fig.axes)):
    fig.axes[ii].axis("off")
    
    
#%% load cell

unit = population['TT9\\cl-maze1.3']
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
fnum=0
field = unit.fields[fnum]
perim = unit.perimeters[fnum]

#%% visualize individual turns w behavior 
plt.figure()
drawTrack()
tnum = 879
turn = turns.iloc[tnum]
ts_start, ts_end = float(turns.iloc[tnum-1]['Ts entry']), float(turns.iloc[tnum+1]['Ts exit'])
behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
plt.plot(perim[:,0], perim[:,1],c='r',zorder=99)
plt.plot(behav[:,1], behav[:,2], c='k')
plt.scatter(behav[0,1],behav[0,2],c='g',marker='o',s=20,zorder=99)
plt.scatter(behav[-1,1],behav[-1,2],c='r',marker='o',s=20,zorder=99)
plt.title(f"{unit.name} Field {fnum} Turn {turn.name}")
