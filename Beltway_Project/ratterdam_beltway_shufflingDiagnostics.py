# -*- coding: utf-8 -*-
"""
Created on Wed May 26 18:53:45 2021

@author: whockei1

Visualizing Shuffling and Trial Assortment Stats
"""

import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt
import pandas as pd
from collections import OrderedDict
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


rat = 'R859'
expCode = 'BRD5'
datafile = f'E:\\Ratterdam\\{rat}\\{rat}{expCode}\\'


Def.includeRewards = 0
qualThresh = 3
alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
clustlist, clustQuals = util.getClustList(datafile) # clustList and clustQuals share same order. ith entry in each is the name and qual of same cell. 
population = OrderedDict()

#%%
clust = 'TT1\\cl-maze1.6'
alley = 1

unit = Core.UnitData(clust, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
unit.loadData_raw()
rms = np.empty((0, Def.singleAlleyBins[0]-1))
for visit in unit.alleys[alley]:
    rm = visit['ratemap1d']
    rms = np.vstack((rms, rm))
    
cmap = util.makeCustomColormap()


rms = np.asarray(rms)

stims =[]
for visit in unit.alleys[alley]:
    stims.append(visit['metadata']['stimulus'])

stims = np.asarray(stims)   

#%%
    
#Convert trials to percent-of-max and plot as a pseudo-cdf to
# find jumps.  
rmmax = np.nanmax(rms[np.nanargmax(np.nanmax(rms,axis=1)),:])
pctMax = sorted(np.nanmax(rms,axis=1)/rmmax)
plt.figure(figsize=(8,8))
plt.title(f"{clust}, Alley {alley}")
plt.plot(pctMax, '.')
plt.ylabel("Max FR Bin as Pct Max Session",fontsize=20)
plt.xlabel("Trial", fontsize=20)
plt.figure(figsize=(8,8))
plt.title(f"{clust}, Alley {alley}")
plt.imshow(rms, aspect='auto',interpolation='None',cmap=cmap,origin='lower')
plt.xlabel("Spatial Bin", fontsize=20)
plt.ylabel("Trial", fontsize=20)
# plt.figure()
# plt.plot(np.gradient(pctMax))
# plt.title(f"{clust}, Alley {alley}")
# print(f"Max grad value: {max(np.gradient(pctMax))}")

def shuffleAndView(rms, stims):
    
    stims = np.random.permutation(stims)
    dataByStim = {'rms':None, 'avg':None, 'err':None}
    
    fig = plt.figure(constrained_layout=True, figsize=(10,6))
    gs = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
    f_ax1 = fig.add_subplot(gs[0,:])
    f_ax2 = fig.add_subplot(gs[1, 0])
    f_ax3 = fig.add_subplot(gs[1, 1])
    f_ax4 = fig.add_subplot(gs[1, 2])
    
    
    
    for i, (ax, txt, c) in enumerate(zip([f_ax2, f_ax3, f_ax4], ['A','B','C'], ['r','b', 'g'])):
        dataByStim['rms'] = rms[np.where(stims==txt)[0],:]
        mask = np.ma.masked_invalid(rms[np.where(stims==txt)[0],:])
        dataByStim['avg'] = mask.mean(axis=0)
        dataByStim['err'] = np.std(mask,axis=0)/np.sqrt(rms[np.where(stims==txt)[0],:].shape[0])
        
        f_ax1.plot(dataByStim['avg'],color=c)
        f_ax1.fill_between(range(len(dataByStim['avg'])), dataByStim['avg']-dataByStim['err'], dataByStim['avg']+dataByStim['err'], color=c, alpha=0.5)
        max_for_imshow = np.nanpercentile(rms, Def.singlealley_imshow_pct_cutoff)
        ax.imshow(dataByStim['rms'], aspect='auto', interpolation='None', cmap=cmap, vmax=max_for_imshow)
        ax.set_title(txt)







