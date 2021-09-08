# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 17:40:46 2021

@author: whockei1

Exploratory and prototpying analysis of winding number behavioral analysis
and how these data relate to (repeating) place field time dynamics
"""

from matplotlib import pyplot as plt
import numpy as np
import ratterdam_RepetitionCoreFx as RepCore
import utility_fx as util
from matplotlib import path
import williamDefaults as wmDef
import os, json
import more_itertools
from scipy.interpolate import PchipInterpolator as pchip



#%% plot field against winding # around blocks

clustname = "TT3\\cl-maze1.6"
f = 0
unit = population[clustname]
field = unit.fields[f]

blocknum = 1
b = w[f"Block_{blocknum}"]
ts = w['Timestamp']
sigma = 1 # smoothing sigma for the field rate vector

fig, ax = plt.subplots()
ax.plot(ts,b,'k')
ax2 = ax.twinx()
ax2.plot(field[:,0],util.weird_smooth(field[:,1],sigma),'b',marker='.',markersize=10)
ax.tick_params(axis='y',labelcolor='k')
ax.set_ylabel("Winding Number",color='k')
ax2.tick_params(axis='y',labelcolor='b')
ax2.set_ylabel("Firing Rate (Hz)",color='b')


#%%
nwins = 10
fig, ax = plt.subplots(2,3)
ts = (w['Timestamp']-unit.position[0,0])/1e6
sessionStart, sessionStop = 0, (unit.position[-1,0]-unit.position[0,0])/1e6 # better ts are in the sessionepochinfo file but this okay
windowEnds = np.linspace(sessionStart, sessionStop, num=nwins)
interpolate = False
for iax, bb in enumerate([1,2,3,4,5,6]):
    print(bb)
    blocknum = bb
    block = w[f"Block_{blocknum}"]  
    corrs = np.empty((0,2))
    
    for unitname, unit in population.items():
        for field in unit.fields:
            fmax = max(field[:,1])
            
            #monotonic spline interpolation
            interpField = pchip(field[:,0], field[:,1])
            
            for i in range(windowEnds.shape[0]-1):
                wStart, wStop = windowEnds[i], windowEnds[i+1]
                
                if interpolate:
                    xx = np.linspace(wStart, wStop, 50)
                    rateSlice = interpField(xx)
                    avgRate = np.mean(rateSlice)/fmax
                else:
                    rateSlice = field[(field[:,0] > wStart) & (field[:,0] <= wStop)]
                    avgRate = np.mean(rateSlice[:,1])/fmax
                windingSlice = block[(ts > wStart) & (ts <= wStop)]
                
                avgWinding = np.mean(windingSlice)
                corrs = np.vstack((corrs,[avgRate, avgWinding]))
                

    fig.axes[iax].scatter(corrs[:,0], corrs[:,1])
    fig.axes[iax].set_title(f"Block {blocknum}")
    fig.axes[iax].set_xlabel("Avg Rate in Window")
    fig.axes[iax].set_ylabel("Avg Winding Number in Window")
plt.suptitle(f"{rat} {day} Avg Winding Number vs Rate in Each Field, {nwins} windows, interp={interpolate}")