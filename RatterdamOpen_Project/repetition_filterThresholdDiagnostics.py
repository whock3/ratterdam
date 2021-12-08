# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 15:46:51 2021

@author: whockei1

Script to visualize each field traversal
and label whether it will be included based on the filter params

Purpose of script is to debug/test filter thresholds and see their effect
on data

"""

import numpy as np
import utility_fx as util
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import newAlleyBounds as nab
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import repeatingPC as repPC
import matplotlib.path as path
import copy

plt.ion()

#%% Load unit and turns 
rat, day = 'R859', 'D2'
clust = 'TT7\\cl-maze1.5'
ratborders = nab.loadAlleyBounds(rat, day)
turns, _ = RepCore.loadTurns(rat,day)

ballisticTurnIdx = []
for i in range(1,turns.shape[0]-1):
   row = turns.iloc[i]
   inter = row['Inter']
   if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter:
       ballisticTurnIdx.append(i)

refturns = copy.deepcopy(turns) # keep a copy without filtering.
turns = turns.iloc[np.asarray(ballisticTurnIdx)]

unit = RepCore.loadRepeatingUnit(rat, day, clust)
RepCore.plotRoutine_RepPF_TempDyn(unit)
repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)


#%% 

## traversal params
"""
    length_thresh : float, optional
        Define a ratio of the length between the first point
        of the trajectory and the furthest pt away, to either
        dista or distb. if eithe ratio is greather than this thresh
        then True. The default is 0.4.
    dist_thresh : float, optional
        Min dist from a point of trajectory to the border.
        If the dist is greater than thresh, True,else False
        The default is 0.1.
    dist_point_thresh : int, optional
        How many points must be dist_thresh away. The default is 3.
    inside_point_thresh : int, optional
        What # of points must be inside field. The default is 3.
"""

length_thresh=0.25
dist_thresh=0.1
dist_point_thresh=2
inside_point_thresh=2

    
fnum = 1

field = unit.fields[fnum]
perim = unit.perimeters[fnum]
foverlap = overlaps[fnum]
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

traversals = [] # list of dicts. Each dict has key corresponding to traversal
                # data

for tnum, turn in turns.iterrows():
    if tnum < refturns.shape[0]-1:
    
        # Because turn definition is staggered as it moves across track,
        # pick either alley+ or alley- by convention and make the ts, labels
        # match that choice 
        
        if turn['Alley+'] in foverlap: 
            
            ts_start, ts_end = float(turn['Ts entry']), float(refturns.iloc[tnum+1]['Ts exit'])
            behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
            behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
            
            filtOutcome = RepCore.filterVisit(dista,
                                              distb,
                                              behav,
                                              perim,
                                              length_thresh=length_thresh,
                                              dist_thresh=dist_thresh,
                                              dist_point_thresh=dist_point_thresh,
                                              inside_point_thresh=inside_point_thresh)
            
            traversals.append({'tnum':tnum,'behavior':behav,'passFilter':filtOutcome})
            

ncol = 10
pad = 10
fig, _ax  = plt.subplots(int(np.ceil(len(traversals)/ncol)), ncol, figsize=(15,15))

for i,traversal in enumerate(traversals):
    
    ax = fig.axes[i]
    ax.plot(perim[:,0], perim[:,1], color='grey')
    ax.plot(traversal['behavior'][:,1], traversal['behavior'][:,2], color='k')
    ax.scatter(traversal['behavior'][0,1], traversal['behavior'][0,2],marker='o',s=30,color='k')
    ax.scatter(traversal['behavior'][-1,1], traversal['behavior'][-1,2],marker='^',s=30,color='k')
    if traversal['passFilter'] == True:
        ax.patch.set_color('green') # or whatever color you like
        ax.patch.set_alpha(0.5)
    elif traversal['passFilter'] == False:
        ax.patch.set_color('red') # or whatever color you like
        ax.patch.set_alpha(0.5)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_title(f"{i},{traversal['tnum']}")

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
prop_passing = len([1 for traversal in traversals if traversal['passFilter']==True])/len(traversals)
plt.suptitle(f"{rat} {day} {clust} Field {fnum}, {round(prop_passing,2)} traversals included")


#%% plot individual traversal
tr=87
traversal = traversals[tr]

fig, ax = plt.subplots()


ax.plot(perim[:,0], perim[:,1], color='grey',zorder=99)
ax.plot(traversal['behavior'][:,1], traversal['behavior'][:,2], color='k',zorder=99)
ax.scatter(traversal['behavior'][0,1], traversal['behavior'][0,2],marker='o',s=30,color='k',zorder=99)
ax.scatter(traversal['behavior'][-1,1], traversal['behavior'][-1,2],marker='^',s=30,color='k',zorder=99)
util.drawTrack(rat,day,ax=ax)
if traversal['passFilter'] == True:
    ax.patch.set_color('green') # or whatever color you like
    ax.patch.set_alpha(0.5)
elif traversal['passFilter'] == False:
    ax.patch.set_color('red') # or whatever color you like
    ax.patch.set_alpha(0.5)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_title(f"{tr},{traversal['tnum']}")