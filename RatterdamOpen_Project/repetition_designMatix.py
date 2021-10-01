# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:29:40 2021

@author: whockei1

Script to create and assess statistics of design matrix 
- of cross rat/day
- of single rat/day
- viewing pop response across behavioral types 
"""
import pandas as pd, numpy as np
from matplotlib import pyplot as plt 
import ratterdam_Defaults as Def 
import utility_fx as util
from collections import Counter

prevdirlevels = ["N","E","S","W"]
nextdirlevels = ["N","E","S","W"]
currdirlevels = ["N","E","S","W"]
prospegolevels = ["S","R","B","L"]
retroegolevels = ["S","R","B","L"]

# defining the design space from the csv below, or anything like it, is wrong because it takes into account
# field locations. It's defining the FR response and labels for visits to fields. But we want ethogram of all behaviors
# (defined as three turns together) regardless of field locations. (also this overcounts bc >>1 field at some locations)
df = pd.read_csv('E:\\Ratterdam\\R_data_repetition\\20210924-145911_superPopAlleyBehaviorResponse_1.5vfilt.csv')

rat, day = 'R781', 'D3'
rdf = df[(df['Rat']==rat)&(df['Day']==day)]

designSpace = {}
for ld in prevdirlevels:
    for cd in currdirlevels:
        for nd in nextdirlevels:
            for re in retroegolevels:
                for pe in prospegolevels:
                    designSpace[f'{ld}{cd}{nd}{re}{pe}'] = 0
                    

for tnum, traversal in rdf.iterrows():
    designSpace[f"{traversal['PrevDir']}{traversal['CurrDir']}{traversal['NextDir']}{traversal['RetroEgo']}{traversal['ProspEgo']}"] += 1
    

validDesignSpace = {}
      
for k,v in designSpace.items():
    if v >= 1:
        validDesignSpace[k] = v
        
        
        
        
plt.figure()
plt.bar(range(len(validDesignSpace)), validDesignSpace.values())
plt.xticks(ticks=range(len(validDesignSpace)),labels=validDesignSpace.keys(),rotation=45)
plt.ylabel("Frequency",fontsize=20)
plt.title(f"{rat}{day} Population Obtained Design Space (Alleys)",fontsize=24)
plt.text(8,800,s="Key:\nPrevDir > CurrDir > NextDir > RetroEgo > ProspEgo",fontsize=14)
for i,(k,v) in enumerate(validDesignSpace.items()):
    plt.text(i-0.2,v+1.5,s=v,fontsize=12)
    

#%% 

#make dict mapping unique color to cellname to visualize cross plots
cmap = plt.get_cmap('gist_ncar')
colors = cmap(np.linspace(0, 1, np.unique(rdf['CellID']).shape[0]))
cdict = {}
for c, uname in zip(colors,np.unique(rdf['CellName'])):
    cdict[uname]=c
    
    
# scatter of neural response by code type, unique color per cell (see above)
fig, ax = plt.subplots(int(np.ceil(np.unique(rdf['code']).shape[0]/5)),5)
for i,(cname, codegroup) in enumerate(rdf.groupby('code')):
    for uname, ugroup in codegroup.groupby('CellName'):
        fig.axes[i].plot(ugroup.StartTime,ugroup.Rate,linestyle='',marker='o',color=cdict[uname])
        
# figure for each alley. subplot for each behavior. bar for each cell mean. +/- std
for a, agroup in rdf.groupby('Alley'):
    fig, ax = plt.subplots(int(np.ceil(agroup['code'].unique().shape[0]/5)),5,sharey=True,figsize=(15,10))
    for i,(cname, codegroup) in enumerate(agroup.groupby('code')):
        means, stds, names = [], [], []
        for uname, ugroup in codegroup.groupby('FieldID'):
            means.append(np.mean(ugroup.Rate))
            stds.append(sem(ugroup.Rate))
            names.append(uname)
        fig.axes[i].bar(range(len(means)),means,yerr=stds)
        fig.axes[i].set_xticks(range(len(names)))
        fig.axes[i].set_xticklabels(names,fontsize=8,rotation=45)
        fig.axes[i].set_ylabel("Mean +/- SEM",fontsize=12)
        fig.axes[i].set_title(cname)
        fig.axes[i].grid(True)
    plt.suptitle(f"Alley {a}")
    plt.savefig("E:\\Ratterdam\\temp\\ethogramsByField\\"+f"{rat}{day}_Alley{a}_TrajResponse.png",dpi=300)
    plt.close()
    
    
#%% Look at ratemaps with data taken from periods of different behavioral/ethographic frequency
cmap = util.makeCustomColormap()
rat, day = 'R859', 'D1'
turns = superpopulation[rat][day]['turns']
refturns = superpopulation[rat][day]['refturns']
unit = copy.deepcopy(superpopulation[rat][day]['units']['TT8\\cl-maze1.1'])
unit.spikes = unit.spikes[unit.spikes[:,2]<400]
unit.position = unit.position[unit.position[:,2]<400]

codes = []
for t, turn in refturns.iterrows():
    if t < refturns.shape[0]-2:
       # reminder data is based on 'current' alley being alley+ in the turn df 
        code = (f"{Def.allocodedict[turn['Allo-']]}"
            f"{Def.allocodedict[turn['Allo+']]}"
            f"{Def.allocodedict[refturns.iloc[t+1]['Allo+']]}"
            f"{Def.egocodedict[turn['Ego']]}"
            f"{Def.egocodedict[refturns.iloc[t+1]['Ego']]}")
    
        # code = (f"{Def.allocodedict[turn['Allo+']]}"
        #         f"{Def.allocodedict[refturns.iloc[t+1]['Allo+']]}"
        #         )
    
        codes.append(code)

behaviorSpace = Counter(codes)

freqIntervals = [0, 10, 50, 100]
ethogroups = {i:[] for i in range(len(freqIntervals)-1)}
ethorms = {i:{'spikes':np.empty((0,3)), 'occs':np.empty((0,3))} for i in range(len(freqIntervals)-1)}

# group behaviors according to their frequency, intervals coded above (by visual inspection right now)
for i in range(len(freqIntervals)-1):
    start, stop = freqIntervals[i], freqIntervals[i+1]
    for k,v in behaviorSpace.items():       
        if start < v <= stop:
            ethogroups[i].append(k)
 
# iterate over (ballistic) turns and for those that are in valid design space,
# add the the data from that turn (i.e. the alley+ part) to the right section
# of behavioral frequency dict
c = 0
for t, turn in refturns.iterrows():
    if t < refturns.shape[0]-2:
        #reminder data is based on 'current' alley being alley+ in the turn df 
        code = (f"{Def.allocodedict[turn['Allo-']]}"
            f"{Def.allocodedict[turn['Allo+']]}"
            f"{Def.allocodedict[refturns.iloc[t+1]['Allo+']]}"
            f"{Def.egocodedict[turn['Ego']]}"
            f"{Def.egocodedict[refturns.iloc[t+1]['Ego']]}")
        #ethogroups, which has the behavioral codes, and ethorms, which has the 
        # data for when the rat did each behavior, share keys so using that below
        for k, v in ethogroups.items():
            if code in v:
                c +=1
                start, stop = float(refturns.iloc[t-1]['Ts entry']), float(refturns.iloc[t+2]['Ts exit'])
                ethorms[k]['spikes'] = np.vstack((ethorms[k]['spikes'], unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)]))
                ethorms[k]['occs'] = np.vstack((ethorms[k]['occs'], unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=stop)]))
        
rms = {k:util.makeRM(v['spikes'],v['occs'],bins=[50,70]) for k,v in ethorms.items()}
fig, ax = plt.subplots(1,3,figsize=(20,9))
for k, v in ethorms.items():
    rm = rms[k]
    fig.axes[k].imshow(rm, origin='lower',aspect='auto',interpolation='None',cmap=cmap,zorder=99,vmax=np.max([np.nanpercentile(i.flatten(),98) for i in rms.values()]))    
    fig.axes[k].set_title(f"{unit.name} Interval {k}")