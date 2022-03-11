# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 12:09:45 2022

@author: whockei1

Script to look at visit filtering and how different choices of thresholds
affects the data
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
import ratterdam_DataFiltering as Filt 
import pickle 


with open("E:\\Ratterdam\\R_data_repetition\\22-02-18_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)

rat, day = 'R859', 'D2'

pop = superpop[rat][day]['units']
turns = superpop[rat][day]['turns']
refturns = superpop[rat][day]['refturns']

cellcounter, fieldcounter = 0, 0
previousDirection, currentDirection, nextDirection, prospectiveEgo, retrospectiveEgo, orientation, location, turnNums, rates = [], [], [], [], [], [], [], [], []
cellnames, fieldnums, cellIDs, fieldIDs = [], [], [], []
repeating = []
startTimes = []
rats, days, alleys = [], [], []
nfields = [] # number of fields a cell has. want to look at things by fieldedness and # fields 
traversal = [] # list of bools corresponding to whether rat went thru alley (True) or turned around (False)
visitRewards = []
allocodedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
egocodedict = {'1':'S','2':'R','3':'B','4':'L','0':'X'}

#define alleys in groups to label visits accordingly later
verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]
trackperimeter = [str(i) for i in [0,4,6,7,9,10,13,15,16,2]]
trackinterior = [str(i) for i in [1,3,14,12,5,11, 8]]

ratborders = nab.loadAlleyBounds(rat, day)
rewards = RepCore.readinRewards(rat, day)
codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}


#%% Pick a unit 
clustName = 'TT7\\cl-maze1.5'
unit = pop[clustName]
repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
#%% Pick field 
fnum = 0
field = unit.fields[fnum]
perim = unit.perimeters[fnum]
foverlap = overlaps[fnum]

foverlap = [str(i) for i in foverlap]
# Field size data, used to filter visits according to threshold pct'ages of field size
maxx,minx, maxy, miny = max(perim[:,0]), min(perim[:,0]), max(perim[:,1]), min(perim[:,1])
dista = np.sqrt(((minx-minx)**2)+((maxy-miny)**2))
distb = np.sqrt(((minx-maxx)**2)+((miny-miny)**2))

contour = path.Path(perim)


#%% Test different choices of visit thresholds 
filtPasses, all_behaviors, alleys, turnNums = [], [], [], []

for tnum, turn in refturns.iterrows():
    # first condition checks we arent at end. 
    # Second checks for discontinuities in tracking progression through alleys
    # commented second check out: and turn['Alley+'] == refturns.iloc[tnum+1]['Alley-'] 
    if tnum < refturns.shape[0]-1:
    
        # Because turn definition is staggered as it moves across track,
        # pick either alley+ or alley- by convention and make the ts, labels
        # match that choice 
        
        if turn['Alley+'] in foverlap: 
            
            ts_start, ts_end = float(turn['Ts entry']), float(refturns.iloc[tnum+1]['Ts exit'])
            behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
            behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
            
            filtOutcome = RepCore.filterVisit(dista,distb,behav,perim,
                                              length_thresh=0.01,
                                              dist_thresh=0.05,
                                              dist_point_thresh=2,
                                              inside_point_thresh=2)
            
            filtPasses.append(filtOutcome)
            all_behaviors.append(behav)
            alleys.append(turn['Alley+'])
            turnNums.append(turn.name)
            
nturns = len(all_behaviors)
ncols = 10
boundingOffset = 20
fig, ax = plt.subplots(int(np.ceil(nturns/ncols)),ncols, figsize=(12,8))

for i,(behav, alley, passBool, tnum) in enumerate(zip(all_behaviors, alleys, filtPasses, turnNums)):
    cax = fig.axes[i]
    cax.plot(behav[:,1], behav[:,2], color='k', zorder=99)
    cax.scatter(behav[0,1],behav[0,2],color='green',edgecolor='k',marker='o',s=75,zorder=99)
    cax.scatter(behav[-1,1],behav[-1,2],color='red',edgecolor='k',marker='o',s=75,zorder=99)
    
    if passBool == True:
        bgcolor = 'palegreen'
    else:
        bgcolor = 'coral'
    
    util.drawRegion(cax, ratborders.alleyInterBounds[alley], color=bgcolor)
    cax.plot(perim[:,0], perim[:,1], color='k')
    
    cax.set_xlim([minx-boundingOffset, maxx+boundingOffset])
    cax.set_ylim([miny-boundingOffset, maxy+boundingOffset])
    
    cax.set_title(tnum)
    
    
for ii in range(len(fig.axes)):
    fig.axes[ii].axis("off")

            
  








