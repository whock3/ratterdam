# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:50:24 2022

@author: whockei1

Analysis of overdispersion in CA1 unit responses
Attempt to relate dispersion of firing rates to behavioral variables
and specifically the distribution of trajectory choices available and/or realized

1-28-22 - attempt based on z-scoring average firing rate across individual passes
         with track broken down into groups of regions with same choice DoF
1-29-22 - attempt based on z-scoring as on 1-28, but x-axis is now number of unique
        kinds of trajectories performed at each region. Alternatively, the 
        spread, in pct, between the max sampled and min sampled trajectory. E.g.
        if there were two sampled paths in 80:20 proportion the x-axis value would be
        0.6
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt, copy
from scipy.stats import zscore

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-135920_superPopAlleyBehaviorResponse_1.5vfilt.csv"
alleydf = pd.read_csv(alleydatapath)

interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
interdf = pd.read_csv(interdatapath)


# r,d = 'R859', 'D2'

# alleydf = alleydf[(alleydf.Rat==r)&(alleydf.Day==d)]
# interdf = interdf[(interdf.Rat==r)&(interdf.Day==d)]

# alleydf = alleydf[alleydf.Repeating==True]
# interdf = interdf[interdf.Repeating==True]

# alleydf = alleydf[alleydf.NumFields==1]
# interdf = interdf[interdf.NumFields==1]

#%% 1-28-22 Overdispersion based on choice DoF at different track locations
# Group alleys and intersections, separately, into groups with the same number
# of paths available through them. Then plot zscores of samples in each group
# and compare

alley_set_mapping = {0:4,
              1:2,
              2:4,
              3:2,
              4:3,
              5:2,
              6:4,
              7:4,
              8:2,
              9:4,
              10:4,
              11:2,
              12:1,
              13:3,
              14:2,
              15:4,
              16:4}


inter_set_mapping = {'A':3,
              'B':2,
              'C':2,
              'D':3,
              'E':2,
              'F':1,
              'G':1,
              'H':2,
              'I':3,
              'J':2,
              'K':2,
              'L':3}

alley_set_zscores = {i:[] for i in range(1,5)}
inter_set_zscores = {i:[] for i in range(1,4)}


# Field Chunk Dispersion for Alleys. If field overlaps multiple
# intersections, analyze each piece separately 
for fid, field in alleydf.groupby("FieldID"):
    #field can spill over multiple alleys within perimeter or interior alleysets
    oriens = np.unique(field.Orientation)
    for o in oriens:
        ofield = field[field.Orientation==o]
        z = zscore(ofield.Rate)
        mapping = alley_set_mapping[np.unique(ofield.Alleys)[0]]
        alley_set_zscores[mapping].extend(z)
    
    
# Field Chunk Dispersion for Intersections. If field overlaps multiple
# intersections (less likely than for alleys), analyze each piece separately 
for fid, field in interdf.groupby("FieldID"):

    # most fields will overlap only 1 intersection (if any), but some are long enough to overlap 2
    inters = np.unique(field.Inters)
    for inter in inters:
        ifield = field[field.Inters==inter]
        z = zscore(ifield.Rate)
        mapping = inter_set_mapping[np.unique(ifield.Inters)[0]]
        inter_set_zscores[mapping].extend(z)
        


#%% Figures

fig, ax = plt.subplots()
labels = []
for i,(setName, zs) in enumerate(alley_set_zscores.items()):
    zs = np.asarray(zs)
    ax.violinplot([zs[~np.isnan(zs)]],positions=[i], points=len(zs[~np.isnan(zs)]), quantiles=[0.25,0.75])
    labels.append(f"Alley {setName}")
for j,(setName, zs) in enumerate(inter_set_zscores.items()):
    zs = np.asarray(zs)
    ax.violinplot([zs[~np.isnan(zs)]],positions=[i+j+1], points=len(zs[~np.isnan(zs)]), quantiles=[0.25,0.75])
    labels.append(f"Inter {setName}")   

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,fontsize=16)
ax.set_title("Firing Rate Dispersion as Function of Track Locations \n Grouped by Choice DoF", fontsize=20)
ax.set_ylabel("Z-score of Average FR over each Field Pass", fontsize=16)


#%% 1-29-22 Plotting dispersion (assayed by z-scores) versus attributes of the 
# paths actually performed. One way is zscores of average firing rates for each pass
# through field versus number of unique trajectories done through that field. 
# That is done here.
# Edit: zscores didn't show anything so I did the same thing with directionality as readout
# For alleys that is unsigned difference in mean FR in dir A vs B. For intersections
# it is first taking mean of FR for each unique traj and taking max-min of those.

# But 1/31 edit I realize that for intersections that doesn't really break direction down so will revise in new cell

zscores = []
ntrajs = []
diffs = []

traj_sample_thresh = 0

for rname, rat in alleydf.groupby("Rat"):
    for daynum, day in rat.groupby("Day"):

        for alleyNum, alley in day.groupby("Alleys"):
            
            alley_trajs = []
            
            for codename, code in alley.groupby("Code"):
                alley_trajs.append(code.shape[0])
                
            #alley_trajs = [i for i in alley_trajs if i > traj_sample_thresh]
            
            alleyNumTrajs = len(alley_trajs)
            alleySpreadTrajs = (max(alley_trajs)-min(alley_trajs))/len(alley_trajs)
            
            for fid, field in alley.groupby("FieldID"):
                
                dirs  = np.unique(field.CurrDir)
                
                if len(dirs)>1:
                    try:
                        diff = abs(field[field.CurrDir==dirs[0]].Rate.mean() - field[field.CurrDir==dirs[1]].Rate.mean())/field.Rate.max()                 
                        diffs.append(diff)
                        ntrajs.append(alleyNumTrajs)
                    except:
                        pass
                
                zscores.extend(zs)
        
     
        
#%% Inters

zscores = []
ntrajs = []
diffs = []

traj_sample_thresh = 0

inter_dir_type = 'CurrEgo' 

for rname, rat in interdf.groupby("Rat"):
    for daynum, day in rat.groupby("Day"):

        for interNum, inter in day.groupby("Inters"):
            
            inter_trajs = []
            
            for codename, code in inter.groupby("Code"):
                inter_trajs.append(code.shape[0])
                
            #alley_trajs = [i for i in alley_trajs if i > traj_sample_thresh]
            
            interNumTrajs = len(inter_trajs)
            interSpreadTrajs = (max(inter_trajs)-min(inter_trajs))/len(inter_trajs)
            
            for fid, field in inter.groupby("FieldID"):
                
                dirs  = np.unique(field[inter_dir_type])
                _diffs = []
                if len(dirs)>1:
                    for d in dirs:
                        
                        _diffs.append(field[field[inter_dir_type]==d].Rate.mean())
                
                    try:
                        diffs.append((max(_diffs)-min(_diffs))/field.Rate.max())
                        ntrajs.append(interNumTrajs)
                    except:
                        pass
                    
                zs = zscore(field.Rate)
                zscores.extend(zs)
#%% 

color = 'red'
        
ntrajdict = {i:None for i in np.unique(ntrajs)}

for nt in np.unique(ntrajs):
    
    ntrajdict[nt] = [diffs[i] for i in range(len(diffs)) if ntrajs[i]==nt]
    
    # plt.violinplot([[zscores[i] for i in range(len(zscores)) if ntrajs[i]==nt]], positions=[nt])
    

    
for k,v in ntrajdict.items():
    
    vp = plt.violinplot([v], positions=[k], widths=1)
    
    vp['bodies'][0].set_facecolor(color)
    vp['bodies'][0].set_edgecolor(color)
    vp['bodies'][0].set_linewidth(2)
    for el in ['cbars','cmaxes','cmins']:
        vp[el].set_color(color)


plt.scatter(ntrajs+np.random.normal(0,0.1,len(ntrajs)), diffs,color=color)


#%% 1-31-22 Directionality based on choice DoF at different track locations
# Group alleys and intersections, separately, into groups with the same number
# of paths available through them. Then plot directionality of samples in each group
# and compare

alley_set_mapping = {0:4,
              1:2,
              2:4,
              3:2,
              4:3,
              5:2,
              6:4,
              7:4,
              8:2,
              9:4,
              10:4,
              11:2,
              12:1,
              13:3,
              14:2,
              15:4,
              16:4}


inter_set_mapping = {'A':3,
              'B':2,
              'C':2,
              'D':3,
              'E':2,
              'F':1,
              'G':1,
              'H':2,
              'I':3,
              'J':2,
              'K':2,
              'L':3}

inter_dir_type = 'PrevDir' 


alley_set_directionality = {i:[] for i in range(1,5)}
inter_set_directionality = {i:[] for i in range(1,4)}


# Field Chunk Dispersion for Alleys. If field overlaps multiple
# intersections, analyze each piece separately 
for fid, field in alleydf.groupby("FieldID"):
    #field can spill over multiple alleys within perimeter or interior alleysets
    oriens = np.unique(field.Orientation)
    for o in oriens:
        ofield = field[field.Orientation==o]
        
        dirs = np.unique(ofield.CurrDir)
        if dirs.shape[0] > 1:
            try:
                diff = abs(field[field.CurrDir==dirs[0]].Rate.mean() - field[field.CurrDir==dirs[1]].Rate.mean())/ofield.Rate.max()
                mapping = alley_set_mapping[np.unique(ofield.Alleys)[0]]
                alley_set_directionality[mapping].append(diff)
            except:
                pass
        
    
    
# Field Chunk Dispersion for Intersections. If field overlaps multiple
# intersections (less likely than for alleys), analyze each piece separately 
for fid, field in interdf.groupby("FieldID"):

    # most fields will overlap only 1 intersection (if any), but some are long enough to overlap 2
    inters = np.unique(field.Inters)
    for inter in inters:
        ifield = field[field.Inters==inter]
        dirs = np.unique(ifield[inter_dir_type])
        if dirs.shape[0] > 1:
            dirmeans = []
            for d in dirs:
                dirmeans.append(ifield[ifield[inter_dir_type]==d].Rate.mean())
            
            mapping = inter_set_mapping[np.unique(ifield.Inters)[0]]
            try:
                inter_set_directionality[mapping].append((max(dirmeans)-min(dirmeans))/ifield.Rate.max())
            except:
                pass
            
            
#%% Figures

fig, ax = plt.subplots()
labels = []
for i,(setName, d) in enumerate(alley_set_directionality.items()):
    d = np.asarray(d)
    ax.violinplot([d[~np.isnan(d)]],positions=[i], points=len(d[~np.isnan(d)]), quantiles=[0.25,0.75])
    ax.scatter(np.asarray([i]*len(d[~np.isnan(d)]))+np.random.normal(0,0.075,len(d[~np.isnan(d)])),d[~np.isnan(d)])
    labels.append(f"Alley {setName}")
for j,(setName, d) in enumerate(inter_set_directionality.items()):
    d = np.asarray(d)
    ax.violinplot([d[~np.isnan(d)]],positions=[i+j+1], points=len(d[~np.isnan(d)]), quantiles=[0.25,0.75])
    ax.scatter(np.asarray([i+j+1]*len(d[~np.isnan(d)]))+np.random.normal(0,0.075,len(d[~np.isnan(d)])),d[~np.isnan(d)])

    labels.append(f"Inter {setName}")   

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,fontsize=16)
ax.set_title("NORMED Firing Rate Directionality as Function of Track Locations \n Grouped by Choice DoF", fontsize=20)
ax.set_ylabel("Unsigned Directional Firing, averaged over each Field Pass", fontsize=16)