# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 13:52:32 2021

@author: whockei1

File to load processed data for figure one and plot each panel
For repetition data paper 

defaults (for e.g. plotting params) in repetition_manuscript_defaults.py
"""

import pickle, numpy as np, os, sys
from matplotlib import pyplot as plt
from scipy.stats import binom_test
import repetition_manuscript_defaults as MDef
from matplotlib.ticker import MaxNLocator

plt.ion()

#Load dict containing all recording day datasets
# structure is rat > day > unit
# unit is a Unit() class 
with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   
    
    
#%% Panel A is a track schematic, no python code used here

#%% Panel B is set of ratemaps showing repetition. no python code used here
# (ratemaps already saved w data, just load and format them)

#%% Panel C - proportion of cells within each recording day defined as repeating
# Reminder that repetition is any (pyramidal) neuron with >=2 fields in a common
# type of location, with no penalty for fields pattern breaking locations.
# locations are vertical alleys, horizontal alleys, intersections

rat_reps = {rat:{day:[] for day in superpop[rat].keys()} for rat in superpop.keys()}

for rat in superpop.keys():
    for day in superpop[rat].keys():
        reps = []
        for unit in superpop[rat][day]['units'].values():
            reps.append(unit.repeating)
        rat_reps[rat][day] = reps
            

rvals = []
rlabels = []
rsigs = [] # do binom test per day
for rat in rat_reps.keys():
    for day in rat_reps[rat].keys():
        rvals.append(sum(rat_reps[rat][day])/len(rat_reps[rat][day]))
        rlabels.append(f"{rat}{day}")
        rsigs.append(binom_test(sum(rat_reps[rat][day]),len(rat_reps[rat][day]),0.05,'greater'))
        
        
        
fig, _ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
w=0.75
ax.bar(range(len(rvals)),rvals,width=w, edgecolor='black',color='grey',linewidth=2)
for i,sig in enumerate(rsigs):
    if sig<0.05:
        ax.text(i-0.1,rvals[i]+0.015,'*',size=70,c='k')
ax.set_xticks(range(len(rlabels)))
ax.set_xticklabels(rlabels,rotation=90)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.hlines(0.05,-0.5,9,color='k',linestyle='--',linewidth=2)
ax.set_ylabel("Proportion of Repeating Cells", fontsize=MDef.ylabelsize)
plt.tight_layout()

#%% Panel D - Number of fields per cell

nfields = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if unit.repeating == True:
                nfields.append(len(unit.perimeters))
                

fig, _ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
ax.hist(nfields,facecolor='grey',edgecolor='k',linewidth=2,bins=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Number of Fields Per Repeating Cell",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))


#%% Panel E - Repetition Orientation Bias
# This is a quantity that varies between -1 and +1. It measures how aligned
# fields are within a cell in terms of whether they share an alley orientation
# Extremes are mostly H,V respectively. 0 means fields equally likely to be V,H within cell

#define V,H alleys
verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

orientationBias = []

for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if len(unit.fields)>1:
                unitbias = []
                for foverlap in unit.overlaps:
                    for regionoverlap in foverlap:
                        if regionoverlap in verticals:
                            unitbias.append(1)
                        elif regionoverlap in horizontals:
                            unitbias.append(-1)
                
                if len(unitbias)>0:
                    orientationBias.append(sum(unitbias)/len(unitbias))
                    

fig, ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
ax.hist(orientationBias,bins=np.linspace(-1,1,num=10),
        facecolor='grey',
        edgecolor='black',
        linewidth=2)
ax.set_xticks([-1,0,1])
ax.set_xticklabels(["-1 (H)", "0", "+1 (V)"])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))


#%% Correlating orientation bias with directionality
# This should not go here in fig 1, but putting in here until I figure out
# where I want to put it. 
import pandas as pd, pickle

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  

datapath  = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

#define V,H alleys
verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

orientationBias = []
meanDiffs = []

for rat in superpop.keys():
    for day in superpop[rat].keys():
        
        rdf = df[(df.Rat==rat)&(df.Day==day)]
        for unit in superpop[rat][day]['units'].values():
            if len(unit.fields)>1:
                unitbias = []
                for foverlap in unit.overlaps:
                    for regionoverlap in foverlap:
                        if regionoverlap in verticals:
                            unitbias.append(1)
                        elif regionoverlap in horizontals:
                            unitbias.append(-1)
                
                # I think this is a crude alley check (as overlaps with inters will lead to empty unitbias list)
                if len(unitbias)>0:
                    orientationBias.append(sum(unitbias)/len(unitbias))
                    celldf = rdf[rdf.CellName==unit.name]
                    meanDiff = []
                    for orien in ['V','H']:
                        odf = celldf[celldf.Orientation==orien]
                        for fname, fgroup in odf.groupby("FieldID"):
                            dirs = np.unique(fgroup.CurrDir)
                    
                            dirA = fgroup[fgroup.CurrDir==dirs[0]]
                            dirB = fgroup[fgroup.CurrDir==dirs[1]]
                            try:
                                meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
                            except:
                                meanDiff.append(0)
                    meanDiffs.append(meanDiff)
                    

fig, ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
ax.hist(orientationBias,bins=np.linspace(-1,1,num=10),
        facecolor='grey',
        edgecolor='black',
        linewidth=2)
ax.set_xticks([-1,0,1])
ax.set_xticklabels(["-1 (H)", "0", "+1 (V)"])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))


  # py history of sufficient-I-think (SIT) code to visualize correlation      
for i,ob in enumerate(orientationBias):
    for md in meanDiffs[i]:
        plt.scatter(abs(ob),md,color='k',alpha=0.7)
        
        
oriendict = {i:[] for i in np.unique([abs(i) for i in orientationBias])}
for i,ob in enumerate(orientationBias):
    for md in meanDiffs[i]:
        oriendict[abs(ob)].append(md)
    
for i,v in oriendict.items():
    plt.violinplot([v],positions=[i])