# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 13:52:32 2021

@author: whockei1


Script to generate the panels for Figure 2 of repetition manuscript
Each panel in the paper corresponds to a code cell here. Code cells in python
begin with "#%%". They can be run independently. The whole script can be run to generate all panels.

Figure 2 is about basic characterization of place cell repetition as it was observed in these data
Panel A  - Percent of recorded neurons that were defined as repeating based on operational definition. This is defined
         for each recording day
Panel B - Number of place fields per repeating neuron
Panel C - Distribution of each cell's orientation bias score. This measures proportion of fields of a cell that share orientation
Panel D - Orientation alignment score. Distribution of shuffled scores, 95th percentile of shuffle, real mean OAS
"""
#%% Imports and Load Data 
import pickle, numpy as np, os, sys
from matplotlib import pyplot as plt
from scipy.stats import binom_test
import repetition_manuscript_defaults as MDef
from matplotlib.ticker import MaxNLocator
import pandas as pd
from collections import Counter

plt.ion()

%matplotlib qt5
#Load dict containing all recording day datasets
# structure is rat > day > unit
# unit is a Unit() class 
with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   
    
alleydf = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv")

#%% Panel A - Proportion of Operationally-defined Repeating Neurons across datasets

rat_reps = {rat:{day:[] for day in superpop[rat].keys()} for rat in superpop.keys()}

for rat in superpop.keys():
    for day in superpop[rat].keys():
        reps = []
        for unit in superpop[rat][day]['units'].values():
            reps.append(unit.repeating)
        rat_reps[rat][day] = reps
            

rvals = []
rlabels = []
for rat in rat_reps.keys():
    for day in rat_reps[rat].keys():
        rvals.append(sum(rat_reps[rat][day])/len(rat_reps[rat][day]))
        rlabels.append(f"{rat}{day}")        

fig, ax = plt.subplots(figsize=(20,15))
w=0.75
ax.bar(range(len(rvals)),rvals,width=w, edgecolor='black',color='grey',linewidth=2)
ax.set_xticks(range(len(rlabels)))
ax.set_xticklabels(rlabels,rotation=90, fontsize=MDef.xlabelsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.hlines(0.05,-0.5,9,color='k',linestyle='--',linewidth=1)
ax.set_ylabel("Proportion of Repeating Cells", fontsize=MDef.ylabelsize)
plt.tight_layout()
#%% Panel B - Number of fields per repeating cell

nfields = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if unit.repeating == True:
                nfields.append(len(unit.perimeters))
                

fig, _ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
labels, counts = np.unique(nfields, return_counts=True)
ax.bar(labels, counts, 
                    align='center', 
                    facecolor='grey',
                    edgecolor='k',
                    linewidth=1)
#ax.hist(nfields,facecolor='grey',edgecolor='k',linewidth=1,bins=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Number of Fields Per Repeating Cell",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))


#%% Panel C - Repetition Orientation Bias

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
                    orientationBias.append(abs(sum(unitbias)/len(unitbias)))
                    
# next cell has this graph plus the shuffle dist results on same plot 
fig, ax = plt.subplots(figsize=(20,15))
ax = fig.axes[0]
labels, counts = np.unique(orientationBias, return_counts = True)
ax.bar(labels, counts, 
                    facecolor='grey',
                    edgecolor='black',
                    label='Orientation Bias Score', 
                    linewidth=1,
                    align='center')

# histout = ax.hist(orientationBias,bins=np.linspace(0,1,num=10), # save bins for shuffling plot below 
#         facecolor='grey',
#         edgecolor='black',
#         label = "Orientation Bias Score",
#         linewidth=1)
ax.set_xticks([0,0.5,1])
ax.set_xticklabels([0,0.5,1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=MDef.ticksize)
ax.set_ylabel("Neuron Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Repetition Orientation Bias",fontsize=MDef.xlabelsize)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
lgnd = plt.legend(prop={'size':MDef.legend_size})
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)


#%% Panel D - Quantifying repetition via OAS
def calcPossibleAlignments(n):
    """
    For n = number of fields
    calculate possible alignment values
    Only considers unique types of patterns eg VVH is same as HHV
    """
    possible_oas = [(n-i)/n for i in range(int(np.floor(n/2))+1)]
    return possible_oas

def findOASRanking(alignment, n):
    pOAS = calcPossibleAlignments(n)

    #top of ratio gets you index in the ranked possible alignment scores
    # bottom normalizes by how many steps there are
    rOAS = (list(reversed(pOAS)).index(alignment)+1)/len(pOAS)
    return rOAS

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
                      
                if cellfieldcounter > 1:
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
        
        sobias = max(Counter(selected_fields).values())/len(selected_fields)

        srOAS = findOASRanking(sobias, Nfields)

        shuffOrienBiases.append(srOAS)
     
        for i in selected_fields_idx:
            argi = np.where(field_indices==i)
            field_indices = np.delete(field_indices, argi)
                          
    allShuffs.append(np.mean(shuffOrienBiases))


# Calc real OAS using above method 

field_nums = []
oas_scores = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            if len(unit.fields) > 1:
                fields = []
                for foverlap in unit.overlaps:
                    for regionoverlap in foverlap:
                        if regionoverlap in verticals:
                            fields.append(1)
                        
                        elif regionoverlap in horizontals:
                            fields.append(-1)
                      
                if len(fields) > 1:
                    field_nums.append(len(fields))
                    obias = max(Counter(fields).values())/len(fields)
                    rOAS = findOASRanking(obias, len(fields))
                    oas_scores.append(rOAS)
#  plotting 
bins = np.linspace(0.5,0.75,25)

fig, ax= plt.subplots()
ax.hist(allShuffs,
                bins=bins,
                facecolor='grey',
                edgecolor='black',
                linewidth=1,
                density=True,
                label = 'Shuffled OAS'
                )

ax.vlines(np.mean(oas_scores),0,ax.get_ylim()[1], 
                color='r')

ax.vlines(np.percentile(allShuffs, 95), 0, ax.get_ylim()[1], 
                color='k')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Frequency",fontsize=MDef.ylabelsize)
ax.set_xlabel("Orientation Alignment Score", fontsize=MDef.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)

lgnd = plt.legend(prop={'size':MDef.legend_size})


plt.show()
# %%
