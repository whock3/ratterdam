# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 18:18:26 2021

@author: whockei1

paper fig 2: characterizing basic repetition properties
- number per day/rat
- distribution of # repeating fields per cell
- bias towards being in V,H alleys exclusively 
"""
import numpy as np, matplotlib.pyplot as plt, pickle, os
import copy
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import newAlleyBounds as nab
import repeatingPC as repPC
import utility_fx as util


#%%  Setting up superpop raw and saving 
failsafe
superpop = {}
ts = util.genTimestamp()

for rat, day in zip(['R765', 'R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886'], ['RFD5', 'D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1','D2']):
    print(rat, day)
    df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    
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


    if rat not in superpop.keys():
        superpop[rat] = {}
        
    superpop[rat][day] = {}
    superpop[rat][day]['turns'] = turns
    superpop[rat][day]['refturns'] = refturns
    superpop[rat][day]['units'] = population


with open(f"E:\\Ratterdam\\R_data_repetition\\{ts}_superpopulationRepetition.pickle","wb") as f:
    pickle.dump(superpop,f)


#%% Load superpop 
savePath = 'E:\\Ratterdam\\2021_SfNPoster_WH\\Fig2_Characterizing\\'

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
    
#%% Number / % repeating cells across rats and days 
from scipy.stats import binom_test

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
        ax.text(i-0.1,rvals[i]+0.015,'*',size=50,c='k')
ax.set_xticks(range(len(rlabels)))
ax.set_xticklabels(rlabels)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=22)
ax.hlines(0.05,-0.5,9,color='k',linestyle='--',linewidth=2)
ax.set_ylabel("Proportion of Repeating Cells", fontsize=28)


#%% How many fields per repeating cell

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
ax.tick_params(axis='both', labelsize=22)
ax.set_ylabel("Frequency", fontsize=28)
ax.set_xlabel("Number of Fields Per Repeating Cell",fontsize=28)


#%% Peak Fr repeating and nonrepeating fields, as well as field gains in repeating cells

repPeaks = []
nonrepPeaks = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if unit.repeating == True:
                for field in unit.fields:
                    repPeaks.append(max(field[:,1]))
            elif unit.repeating == False:
                #bc not all multifielded cells are repeating
                for field in unit.fields:
                    nonrepPeaks.append(max(field[:,1]))
                    
repPeaks = np.asarray(repPeaks)
nonrepPeaks = np.asarray(nonrepPeaks)
thresh=100 # hz
nHzPerBin = 1 # histogram bins should contain how many hertz, e.g. 2Hz increment per bin

repPeaks = repPeaks[repPeaks<thresh]
nonrepPeaks = nonrepPeaks[nonrepPeaks<thresh]

histMax = max([max(repPeaks),max(nonrepPeaks)])
b=np.linspace(0,histMax,int(np.ceil(histMax/nHzPerBin)))# hist bins

fig, _ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
ax.hist(repPeaks,facecolor='red',bins=b,alpha=0.7,label='Repeating Fields')
ax.hist(nonrepPeaks,facecolor='black',bins=b,alpha=0.7,label='Non-repeating Fields')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=22)
ax.set_ylabel("Frequency", fontsize=28)
ax.set_xlabel("Peak Field Rate (Hz)",fontsize=28)
plt.legend(prop={'size':25})


#%% Within repeating cells, bias towards being V or H only fields
# call V +1 and H -1. So even split would have hist centered on 0

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

orientationBias = []

for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if unit.repeating == True:
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
ax.hist(orientationBias,bins=10,facecolor='grey',edgecolor='black',linewidth=2)
ax.set_xticks([-1,0,1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=22)
ax.set_ylabel("Neuron Frequency", fontsize=28)
ax.set_xlabel("Repetition Orientation Bias",fontsize=28)