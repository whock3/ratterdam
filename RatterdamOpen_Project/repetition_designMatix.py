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

prevdirlevels = ["N","E","S","W"]
nextdirlevels = ["N","E","S","W"]
currdirlevels = ["N","E","S","W"]
prospegolevels = ["S","R","B","L"]
retroegolevels = ["S","R","B","L"]

df = pd.read_csv('E:\\Ratterdam\\R_data_repetition\\20210924-145911_superPopAlleyBehaviorResponse_1.5vfilt.csv')

rat, day = 'R859', 'D1'
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
        fig.axes[i].plot(ugroup.StartTimes,ugroup.Rate,linestyle='',marker='o',color=cdict[uname])
        
# figure for each alley. subplot for each behavior. bar for each cell mean. +/- std
for a, agroup in rdf.groupby('Alley'):
    fig, ax = plt.subplots(int(np.ceil(agroup['code'].unique().shape[0]/5)),5,sharey=True,figsize=(8,8))
    for i,(cname, codegroup) in enumerate(agroup.groupby('code')):
        means, stds, names = [], [], []
        for uname, ugroup in codegroup.groupby('FieldID'):
            means.append(np.mean(ugroup.Rate))
            stds.append(np.std(ugroup.Rate))
            names.append(uname)
        fig.axes[i].bar(range(len(means)),means,yerr=stds)
        fig.axes[i].set_xticks(range(len(names)))
        fig.axes[i].set_xticklabels(names,fontsize=8,rotation=45)
        fig.axes[i].set_ylabel(cname,fontsize=12)
        fig.axes[i].grid(True)
    plt.suptitle(f"Alley {a}")
    
    

