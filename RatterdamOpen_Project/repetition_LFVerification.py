# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:30:52 2021

@author: whockei1

Inspection of place fields displaying leader follower relationships with one another
Goal is to look at behavior and other variables to see if an ostensible 
LF relation is really best explained by something trivial like behavioral patterns (LF in field visits, not neural waves)
"""

#%% Imports and data readin
import numpy as np, pandas as pd, pickle, json
import matplotlib.pyplot as plt 
import utility_fx as util

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   
    
    
#%% 
    # 10-4-0 vs 15-4-0

rat, day = 'R765', 'RFD5'
#represent fields as list of ttclust,field num 
ttA, uA, fA = 10,4,0
ttB, uB, fB = 15,4,0
fieldAname, fieldBname = [f'TT{ttA}\\cl-maze1.{uA}',fA], [f'TT{ttB}\\cl-maze1.{uB}',fB]

population, refturns = superpop[rat][day]['units'], superpop[rat][day]['refturns']

fieldA,perimA,overlapsA = (population[fieldAname[0]].fields[fieldAname[1]], 
                        population[fieldAname[0]].perimeters[fieldAname[1]],
                        population[fieldAname[0]].overlaps[fieldAname[1]]
                        )

fieldB,perimB,overlapsB = (population[fieldBname[0]].fields[fieldBname[1]], 
                 population[fieldBname[0]].perimeters[fieldBname[1]],
                 population[fieldBname[0]].overlaps[fieldBname[1]]
                 )

overlapsA = [str(i) for i in overlapsA]
overlapsB = [str(i) for i in overlapsB]

fieldAVisits, fieldBVisits = [], []
# use refturns here i think?
for t, turn in refturns.iterrows():
    if turn['Alley+'] in overlapsA or turn['Inter'] in overlapsA:
        fieldAVisits.append(float(turn['Ts entry']))
    elif turn['Alley+'] in overlapsB or turn['Inter'] in overlapsB:
        fieldBVisits.append(float(turn['Ts entry']))


# create an 'order list' or 'leader list' of which field comes first
# basic idea is that alternating back and forth leaves the first field visited
#as the leader field. but any double visits switches leadering to that field
fieldAVisits = [('A',i) for i in fieldAVisits]
fieldBVisits = [('B',i) for i in fieldBVisits]
jointVisits = fieldAVisits + fieldBVisits
jointVisits = sorted(jointVisits,key=lambda x:x[1])

leaderList = [jointVisits[0]] # initialize w first field
for i,visit in enumerate(jointVisits):
    if i < len(jointVisits)-1:
        if visit[0] == jointVisits[i+1][0]:
            leaderList.append(visit)
        else:
            leaderList.append(leaderList[-1])
            
    
leaderPlot = [+1 if i[0]=='A' else -1 for i in leaderList] # convert to numbers for plotting
leaderTimes = [i[1] for i in leaderList] # get times for plotting

fieldAVisits = np.asarray([i[1] for i in fieldAVisits])
fieldBVisits = np.asarray([i[1] for i in fieldBVisits])


fig, _ax = plt.subplots(4,1,figsize=(19,15))

ax = fig.axes[0]
util.drawTrack(rat,day,ax=ax)
ax.plot(perimA[:,0], perimA[:,1],color='k',label=f"{fieldAname[0]} Field {fieldAname[1]}",zorder=99)
ax.plot(perimB[:,0], perimB[:,1],color='r',label=f"{fieldBname[0]} Field {fieldBname[1]}",zorder=99)
ax.legend()

ax = fig.axes[1]
ax.plot(fieldAVisits, [0]*fieldAVisits.shape[0], color='k', marker='o',label=f"{fieldAname[0]} Field {fieldAname[1]} (A)")
ax.plot(fieldBVisits, [0]*fieldBVisits.shape[0], color='r', marker='^', label=f"{fieldBname[0]} Field {fieldBname[1]} (B)")
ax.plot(leaderTimes, leaderPlot, color='grey', linestyle='--',label='Leader Field, +1=A,-1=B')
ax.set_title("Field Visits Over Time")
ax.legend()

ax = fig.axes[2]
ax.plot(fieldAVisits,np.cumsum([1]*fieldAVisits.shape[0]),color='k',marker='o',label=f"{fieldAname[0]} Field {fieldAname[1]}")
ax.plot(fieldBVisits,np.cumsum([1]*fieldBVisits.shape[0]),color='r', marker='^', label=f"{fieldBname[0]} Field {fieldBname[1]}")
ax.set_title("Cumulative Sum of Visits to Each Field")
ax.legend()

ax = fig.axes[3]
ax.plot(fieldA[:,0], fieldA[:,1],color='k',marker='o',label=f"{fieldAname[0]} Field {fieldAname[1]}")
ax.plot(fieldB[:,0], fieldB[:,1], color='r',marker='^',label = f"{fieldBname[0]} Field {fieldBname[1]}")
ax.set_title("Average Firing Rate on Each Visit to Field")
ax.legend()
plt.suptitle(f"Place Field Time Series Rat {rat} Day {day}, indicated fields picked up by LF analysis")


