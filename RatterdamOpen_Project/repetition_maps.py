# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:41:42 2021

@author: whockei1

Map switching across the population 
"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
import pickle
import matplotlib.path as path

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)
    

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

rat, day = 'R781', 'D3'
population, refturns = superpop[rat][day]['units'], superpop[rat][day]['refturns']

#%% Make neural df 
import matplotlib.path as path

# A dict of pandas dfs. Each key is alley. Each value is df with turns and neural responses
alleyResponses = {str(i):None for i in range(17)}

for alley in alleyResponses.keys():
    for t, turn in turns.iterrows():
        if turn['Allo+'] == alley:
            winStart, winEnd = float(turn['Ts entry']), float(turns.iloc[t+1]['Ts exit'])
            for unit in population.values():
                for field, perimeter, overlap in zip(unit.fields, unit.perimeters, unit.overlaps):
                    contour = path.Path(perimeter)
                    field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                    field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                    unit_rates = []
                    for win in time_windows:
                  
                        winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                        winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                        if winPos.shape[0] > 0:
        
        
        
        
        
        

#%%     
                    
                    
dirtime = np.asarray(dirtime)
dirbias = np.asarray(dirbias)

ncol = 1
smooth = 1
fig, _ax = plt.subplots(int(np.ceil(len(fields)/ncol)),ncol,figsize=(10,15))
for i,field in enumerate(fields):
    ax = fig.axes[i]
    ax.plot(field[:,0], util.weird_smooth(field[:,1],smooth),linewidth=2,marker='.')
    ax2 = ax.twinx()
    ax2.plot(dirtime,dirbias,color='k',linewidth=2,marker='.')
    ax.set_ylabel("Firing Rate (Hz)",fontsize=16)
    ax2.set_ylabel("Accumulated Direction Bias", fontsize=16)
    ax.set_title(f"{names[i]}")
plt.suptitle(f"Alley {alley}")
fig.subplots_adjust(hspace=0.5)

#%% plot all fields separately, look for common dynamics

rat, day = 'R781', 'D3'
population, refturns = superpop[rat][day]['units'], superpop[rat][day]['refturns']


fields = []
names = []
for unit in population.values():
    for field in unit.fields:
        fields.append(field)
        names.append(f"{unit.name} Repeating = {unit.repeating}")
        
ncols = 6
fmin = min([min(f[:,0]) for f in fields])
fmax = max([max(f[:,0]) for f in fields])
fig, _ax = plt.subplots(int(np.ceil(len(fields)/ncols)),ncols, figsize=(12,15))
for i, (name,field) in enumerate(zip(names,fields)):
    ax = fig.axes[i]
    ax.plot(field[:,0], field[:,1])
    ax.set_xlim([fmin,fmax])
    ax.set_title(name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
fig.subplots_adjust(wspace=0.1,hspace=1)

#%% Look at percent change in FR as fx of local dir bias

dirWin = 3
frWin = 3

for unit in population.items():
    for field in unit.fields:
        