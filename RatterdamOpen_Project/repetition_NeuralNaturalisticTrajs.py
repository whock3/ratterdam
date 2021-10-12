# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 19:02:33 2021

@author: whockei1

Script to look at what cells are doing during naturalistic trajectories
Defined series of alley traversals puntucated by a turnaround (later, or rewards too)
"""
#%% Imports and load data 
import pickle, numpy as np, matplotlib.pyplot as plt, os

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)

#%% Select day and pullout data, identify turnarounds 

rat, day = 'R781', 'D3'
  
refturns = superpop[rat][day]['refturns']
pop = superpop[rat][day]['units']
turnarounds = []
for tnum, turn in refturns.iterrows():
    if tnum < refturns.shape[0]-1:
        if turn.Ego == '3' or refturns.iloc[tnum+1].Inter == turn.Inter:
            turnarounds.append(tnum)

turnarounds = np.asarray(turnarounds)
alltrajs = []

for trA, trB in zip(turnarounds[:-1], turnarounds[1:]):
    
    tsA = float(refturns.iloc[trA]['Ts exit'])
    tsB = float(refturns.iloc[trB+1]['Ts entry'])
    rates = []
    for unitname, unit in pop.items():
        spk = unit.spikes[(unit.spikes[:,0]>tsA)&(unit.spikes[:,0]<=tsB)]
        pos = unit.position[(unit.position[:,0]>tsA)&(unit.position[:,0]<=tsB)]
        rate = spk/((pos[-1,0]-pos[0,0])/1e6)
        rates.append(rate)
    alltrajs.append(rates)
    
alltrajs = np.asarray(alltrajs)


# above isnt right. need to get firing over time/space for each.
# theyre different lengths so align them on start or end or whatever and theyll  be jagged


#%% Fields against time, annotated w turnarounds
import pickle
import numpy as np
import repeatingPC as repPC
import newAlleyBounds as nab
from matplotlib import path 
import utility_fx as util
ts = util.genTimestamp()
savepath = 'E:\\Ratterdam\\temp\\popveccorrs\\'

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
    
for rat, day in zip(['R765','R781','R781','R808','R808','R859','R859','R886','R886'],['RFD5','D3','D4','D6','D7','D1','D2','D1','D2']):
    population, turns, refturns = superpop[rat][day]['units'], superpop[rat][day]['turns'], superpop[rat][day]['refturns']
    df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    ratborders = nab.loadAlleyBounds(rat, day)
    
    with open(df+"sessionEpochInfo.txt","r") as f:
        data = f.readlines()
    session_start, session_end = [float(i) for i in data[0].split(',')]
    
    # create time windows
    seconds = (session_end - session_start)/1e6
    wnSize= 10*60*1e6 # set num mins you want in a window
    wnStep = 2*60*1e6
    time_windows = []
    stop = False
    begin = session_start
    while not stop:
        a,b = begin, begin + wnSize
        if b <= session_end:
            time_windows.append((a,b))
            begin += wnStep
        else:
            stop = True
    time_windows = np.asarray(time_windows)
    # wnSize = 0.5*60
    # time_windows = np.linspace(session_start, session_end, int(np.ceil(seconds/wnSize)))
    # turnarounds = []
    # for i in range(1,refturns.shape[0]-1):
    #     if refturns.iloc[i]['Ego'] == '3' or refturns.iloc[i+1].Inter == refturns.iloc[i].Inter:
    #         turnarounds.append(refturns.iloc[i]['Ts exit'])
            
    pop_rates = np.empty((0, time_windows.shape[0]))
    repeating = []
    includeTurnaround = []
    repeating = []
    for _,unit in population.items():
        repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
        for i in range(1): #perim in unit.perimeters:
            repeating.append(repeat)
            #contour = path.Path(perim)
            field_pos = unit.position#[contour.contains_points(unit.position[:,1:])]
            field_spikes = unit.spikes#[contour.contains_points(unit.spikes[:,1:])] 
            unit_rates = []
            for win in time_windows:
                winStart, winEnd = win[0], win[1]
                winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                if winPos.shape[0] > 0:
                    winRate = winSpikes.shape[0]/((winEnd-winStart)/1e6)
                else:
                    winRate = np.nan
                unit_rates.append(winRate)
            
            pop_rates = np.vstack((pop_rates, unit_rates))
            
    repeating=np.asarray(repeating)
    fig, ax = plt.subplots(2,1, figsize=(10,10))
    corrs = np.empty((time_windows.shape[0], time_windows.shape[0]))
    # at least one day has no repeating cells
    pop_rates_repeating = pop_rates[repeating==True,:]
    if pop_rates_repeating.shape[0] > 0:
        for i, pi in enumerate(pop_rates_repeating.T):
            for j,pj in enumerate(pop_rates_repeating.T):
                corr = ma.corrcoef(ma.masked_invalid(pi), ma.masked_invalid(pj))[0,1]
                corrs[i,j] = corr  
        fig.axes[0].set_title(f"{rat}{day} Correlating pop vector across time, 5min wins, 1 min step, Repeating fields ")
        fig.axes[0].imshow(corrs, aspect='auto', interpolation='None')
    
    corrs = np.empty((time_windows.shape[0], time_windows.shape[0]))
    pop_rates_nonrepeating = pop_rates[repeating==False,:]
    for i, pi in enumerate(pop_rates_nonrepeating.T):
        for j,pj in enumerate(pop_rates_nonrepeating.T):
            corr = ma.corrcoef(ma.masked_invalid(pi), ma.masked_invalid(pj))[0,1]
            corrs[i,j] = corr  
    fig.axes[1].set_title(f"{rat}{day} Correlating pop vector across time, Nonrepeating Fields")
    fig.axes[1].imshow(corrs, aspect='auto', interpolation='None')
    plt.savefig(savepath+f"{ts}_{rat}{day}_popcorrs.png", dpi=300)
    plt.close()

# i=0 
# for _, unit in population.items():
#     for perim in unit.perimeters:
#         contour = path.Path(perim)
#         field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
#         plt.scatter(field_spikes[:,0], [i]*field_spikes.shape[0],marker='|',c='k')
#         i += 1
# for t in turnarounds:
#     plt.vlines(float(t),0,25,color='grey')
    



        
# for each field, filter spikes/occs in border. get rates in windows
# do for all fields. annotate w turnaround times (aligned to some halfway mark)
# plot matrix fields y axis, time x axis and rate in each bin. vlines turnarounds

# 21-10-7 this was hard to interpret results. either rasters or rates (rates in small wins too)


#%% Correlating strength of directional tuning of each field with behaviorial directional bias, all rats/days
import pandas as pd
passThresh = 2
df = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv")
for rat in ['R765','R781', 'R808', 'R859', 'R886']:
    biases, ratediffs, errors, repeating = [], [], [], []
    ratdf = df[df['Rat']==rat]
    for unum, ugroup in ratdf.groupby("CellID"):
        for fnum, fgroup in ugroup.groupby("FieldNum"):
            for otype, ogroup in fgroup.groupby("Orientation"):
                directions = ogroup['CurrDir'].unique()
                if len(directions)>2:
                    print("ERROR too many directions")
                elif len(directions)==2:
                    if ogroup['Rate'].mean()>0:
                        da, db = ogroup[ogroup['CurrDir']==directions[0]], ogroup[ogroup['CurrDir']==directions[1]]
                        if da.shape[0] >= passThresh and db.shape[0] >= passThresh:
                            behavioral_bias = max(da.shape[0],db.shape[0])/ogroup.shape[0]
                            diff = abs(da['Rate'].mean()-db['Rate'].mean())
                            err = ogroup["Rate"].sem()
                            biases.append(behavioral_bias)
                            ratediffs.append(diff)
                            errors.append(err)
                            repeating.append(ogroup['Repeating'].unique()[0])
                        
    colors = ['r' if i==True else 'k' for i in repeating]
    plt.figure(figsize=(10,10))
    plt.scatter(biases, ratediffs, c=colors)
    plt.legend(['Repeating','Non-repeating'])
    plt.legend(['Non-repeating','Non-repeating'])
    plt.xlabel("Behavioral Bias",fontsize=16)
    plt.ylabel("Abs Mean Rate Difference btwn Directions",fontsize=16)
    plt.title(f"Place Fields {rat}, Directional Tuning vs Behavioral Bias",fontsize=20)
    plt.savefig(f"E:\\Ratterdam\\temp\\{rat}_dirtuningVsBias.png",dpi=300)
    plt.close()