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
import pandas as pd 

plt.ion()

#Load dict containing all recording day datasets
# structure is rat > day > unit
# unit is a Unit() class 
with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   
    
alleydf = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv")
#%% Panel A is a track schematic, no python code used here

#%% Panel B is set of ratemaps showing repetition. no python code used here
# (ratemaps already saved in each rat's folder in [rat]\ratterdam_plots\[day]\overviewplots\, just load and format them)

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
                    
# next cell has this graph plus the shuffle dist results on same plot 
fig, ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
histout = ax.hist(orientationBias,bins=np.linspace(-1,1,num=10), # save bins for shuffling plot below 
        facecolor='mediumturquoise',
        edgecolor='black',
        label = "Orientation Bias Score",
        linewidth=2)
ax.set_xticks([-1,0,1])
ax.set_xticklabels(["-1 (H)", "0", "+1 (V)"])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
lgnd = plt.legend(prop={'size':MDef.legend_size})
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Repetition Orientation bias - reassignment analysis
# To quantify whether the trimodal distribution we observe is significant
# we do a reassignment analysis where all fields are pooled together and labeled with orientation
# then randomly reassigned to cells while preserving the # fields / cell in real life. Then recompute the percentiles
import copy 

nshuffs = 1000

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

all_fields = []
cell_field_nums = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if len(unit.fields) > 1:
                cellfieldcounter = 0
                for foverlap in unit.overlaps:
                    for regionoverlap in foverlap:
                        if regionoverlap in verticals:
                            all_fields.append(1)
                            cellfieldcounter += 1
                        elif regionoverlap in horizontals:
                            all_fields.append(-1)
                            cellfieldcounter += 1
                      
                if cellfieldcounter > 0:
                    cell_field_nums.append(cellfieldcounter)
                
all_fields = np.asarray(all_fields)
cell_field_nums = np.asarray(cell_field_nums)

allShuffs= []

for s in range(nshuffs):
    if s%50==0:
        print(s)
  
    field_indices = range(len(all_fields))
    
    shuffOrienBiases = []
    
    for Nfields in cell_field_nums:
        
        selected_fields_idx = np.random.choice(field_indices, Nfields, replace = False)
        
        selected_fields = all_fields[selected_fields_idx]
        
        sobias = sum(selected_fields)/len(selected_fields)
        shuffOrienBiases.append(sobias)
     
        #there has to be a better way
        for i in selected_fields_idx:
            argi = np.where(field_indices==i)
            field_indices = np.delete(field_indices, argi)
                          
    allShuffs.append(shuffOrienBiases)

shuffhists = []
for i in allShuffs:
    shuffhists.append(np.histogram(i,bins=np.linspace(-1,1,10))[0])
    

fig, ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
histout = ax.hist(orientationBias,bins=np.linspace(-1,1,num=10), # save bins for shuffling plot below 
        facecolor='mediumturquoise',
        edgecolor='black',
        label = "Orientation Bias Score",
        linewidth=2)

sbarw = histout[-1][0].get_width() # want width of bars used for real orien bias dist 
sbarx = []
for i in range(len(histout[1])-1):
    sbarx.append((histout[1][i] + histout[1][i+1])/2)
    
ax.set_xticks([-1,0,1])
ax.set_xticklabels(["-1 (H)", "0", "+1 (V)"])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
lgnd = plt.legend(prop={'size':MDef.legend_size})
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

ax.bar(sbarx,np.percentile(shuffhists,95,axis=0),
                                            width=sbarw,
                                            align='center',
                                            alpha=0.5,
                                            facecolor='gray',
                                            edgecolor='black',
                                            linewidth=3,
                                            linestyle='--',
                                            label = "$95^{th}$ Percentile Shuffle",
                                            zorder=99)
lgnd = plt.legend(prop={'size':MDef.legend_size})
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Panel F - number horizontal vs vertical fields
from collections import Counter
from scipy.stats import chisquare
oriens = []
for fid, field in alleydf.groupby("FieldID"):
    for orien, ofield in field.groupby("Orientation"):
        oriens.append(orien)
oriens = np.asarray(oriens)
c = Counter(oriens)
print(c)
print(chisquare([c['V'],c['H']]))


#%% Not a panel - number of fields which overlap multiple alleys
# this doesn't take into account fields which extend from an alley into an intersection
# which is almost certainly more common but maybe less what we mean by 'being in multiple regions'
# since theyre adjacent and a field may be normal sized but situated inbetween an alley-inter bound
from collections import Counter

numalleys = []
for fid, field in alleydf.groupby("FieldID"):
    numalleys.append(np.unique(field.Alleys).shape[0])

c = Counter(numalleys)
for i in c.keys():
    print(f"{i} fragments occured  {c[i]} times, {round(c[i]/sum(c.values()),3)*100} % of cases")
#%% Correlating orientation bias with directionality
# This should not go here in fig 1, but putting in here until I figure out
# where I want to put it. 
# import pandas as pd, pickle

# with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
#     superpop = pickle.load(f)  

# datapath  = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
# df = pd.read_csv(datapath)

# #define V,H alleys
# verticals = [2,3,5,7,16,14,11,9]
# horizontals = [0,4,6,1,12,8,15,13,10]

# orientationBias = []
# meanDiffs = []

# for rat in superpop.keys():
#     for day in superpop[rat].keys():
        
#         rdf = df[(df.Rat==rat)&(df.Day==day)]
#         for unit in superpop[rat][day]['units'].values():
#             if len(unit.fields)>1:
#                 unitbias = []
#                 for foverlap in unit.overlaps:
#                     for regionoverlap in foverlap:
#                         if regionoverlap in verticals:
#                             unitbias.append(1)
#                         elif regionoverlap in horizontals:
#                             unitbias.append(-1)
                
#                 # I think this is a crude alley check (as overlaps with inters will lead to empty unitbias list)
#                 if len(unitbias)>0:
#                     orientationBias.append(sum(unitbias)/len(unitbias))
#                     celldf = rdf[rdf.CellName==unit.name]
#                     meanDiff = []
#                     for orien in ['V','H']:
#                         odf = celldf[celldf.Orientation==orien]
#                         for fname, fgroup in odf.groupby("FieldID"):
#                             dirs = np.unique(fgroup.CurrDir)
                    
#                             dirA = fgroup[fgroup.CurrDir==dirs[0]]
#                             dirB = fgroup[fgroup.CurrDir==dirs[1]]
#                             try:
#                                 meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
#                             except:
#                                 meanDiff.append(0)
#                     meanDiffs.append(meanDiff)
                    

# fig, ax = plt.subplots(figsize=(20,15))
# ax = fig.axes[0]
# ax.hist(orientationBias,bins=np.linspace(-1,1,num=10),
#         facecolor='grey',
#         edgecolor='black',
#         linewidth=2)
# ax.set_xticks([-1,0,1])
# ax.set_xticklabels(["-1 (H)", "0", "+1 (V)"])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.tick_params(axis='both', labelsize=MDef.ticksize)
# ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
# ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))


#   # py history of sufficient-I-think (SIT) code to visualize correlation      
# for i,ob in enumerate(orientationBias):
#     for md in meanDiffs[i]:
#         plt.scatter(abs(ob),md,color='k',alpha=0.7)
        
        
# oriendict = {i:[] for i in np.unique([abs(i) for i in orientationBias])}
# for i,ob in enumerate(orientationBias):
#     for md in meanDiffs[i]:
#         oriendict[abs(ob)].append(md)
    
# for i,v in oriendict.items():
#     plt.violinplot([v],positions=[i])