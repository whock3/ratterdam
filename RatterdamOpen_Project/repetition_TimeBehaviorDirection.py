# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:50:07 2021

@author: whockei1

Exploratory analyses looking at possible interaction between
time, behavioral biases, and activity (specifically directionality)
"""

#%% Imports 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
import ratterdam_Defaults as Def 
import ratterdam_RepetitionCoreFx as RepCore

#%% Load data

datapath =  "E:\\Ratterdam\\R_data_repetition\\211210_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)
wnSizeMin = 10
wnStepMin = 5
wnSize= wnSizeMin*60*1e6 # set num mins you want in a window
wnStep = wnStepMin*60*1e6
passThresh=2
    
#%% Generate scatterplot for each rat and day. In each is bias against normed abs mean rate diff
# each datapoint is a timebin. 

biases,diffs = [], []
speeds = [] 
indices = []
pairs = []
    
svs = {}

for rat,day in zip(['R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886', 'R765'],['D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2','RFD5']):
   
    
    daydf = df[(df.Rat==rat)&(df.Day==day)]
    
    # get rid of turnarounds
    daydf = daydf[daydf.Traversal==True]
    
    times = daydf.StartTimes.sort_values()
    session_start, session_end = times.iloc[0], times.iloc[-1]
    
    
    _,u = RepCore.loadTurns(rat,day) # this also returns a unit class obj with position as an attribute, which is what we want
    
    ## Compute smoothed speed curve so we can average within window
    # long term this should be a separate fx, just for now its in this loop
    # input is of course already velocity filtered. This code is taken from Filt.velocity_filterting()
    
    thresh = Def.velocity_filter_thresh
    position = u.position
    gradts, gradx, grady = np.gradient(position[:,0]), np.gradient(position[:,1]), np.gradient(position[:,2])
    winsz=50
    gradx = [np.mean(gradx[0+i:winsz+i]) for i in range(len(gradx))]
    grady = [np.mean(grady[0+i:winsz+i]) for i in range(len(grady))]
    gradx = np.asarray([i/Def.ptsCm for i in gradx])
    grady = np.asarray([i/Def.ptsCm for i in grady])
    
    vx = np.asarray([1e6*(a/b) for a,b in zip(gradx,gradts)])
    vy = np.asarray([1e6*(a/b) for a,b in zip(grady,gradts)])
    v =  np.sqrt((vx**2)+(vy**2))  
    
    sv = [np.mean(v[0+i:winsz+i]) for i in range(len(v))]
    sv = np.column_stack((position[:,0], sv))
    
    #svs[f"{rat}{day}"] = {"v":v,"sv":sv}

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

    for widx, win in enumerate(time_windows):
        
        meanspeed = np.mean(sv[(sv[:,0]>win[0])&(sv[:,0]<=win[1]),1])
        
        for orien in ["V","H"]:
            
            oriendf = daydf[daydf.Orientation==orien]
            
            for fname, field in oriendf.groupby('FieldID'):
                d = np.unique(field.CurrDir)
                if d.shape[0]>1:
                    wf = field[(field.StartTimes>win[0])&(field.StartTimes<=win[1])]
                    dirA, dirB = wf[wf.CurrDir==d[0]], wf[wf.CurrDir==d[1]]
                    if dirA.shape[0] >= passThresh and dirB.shape[0] >= passThresh:
                        try:
                            
                            ## 2-2-22 note: this isnt good as np.mean returns a nan
                            # if any sample is a nan, and matplotlib will just ignore it.
                            # so plots based on this are missing (potentialy a lot) of data.
                            # pandas series.mean() ignores nans by default
                            
                            
                            bias = max(dirA.shape[0]/wf.shape[0],dirB.shape[0]/wf.shape[0])                            
                            diff = abs(np.nanmean(dirA.Rate)-np.nanmean(dirB.Rate))/np.nanmean(field.Rate)
                            biases.append(bias)
                            diffs.append(diff)
                            speeds.append(meanspeed)
                            pairs.append([f"{rat}{day}",widx,bias,diff,dirA,dirB])
                        except Exception as e:
                            print(e)

        
diffs = np.asarray(diffs)
biases = np.asarray(biases)
speeds = np.asarray(speeds)
indices = np.asarray(indices)


#%% 6 subplots showing scatter and heatmap for diffs,bias,directional tuning pairwise comps

#plotting params
_msize, _a, _c =10, 0.6, 'k'
nbins=20
dbins = np.linspace(0,max(diffs),nbins)
bbins = np.linspace(0.5,max(biases),nbins)
sbins = np.linspace(1.5,max(speeds),nbins)

fig, _ax = plt.subplots(1,3,figsize=(12,9))
plt.suptitle(f"Directional Bias vs Directionality, {wnSizeMin}min win, {wnStepMin}min step, all datasets",
             fontsize=30)
ax = fig.axes[0]
ax.scatter(biases,diffs,s=_msize,alpha=_a,c=_c)
ax.set_xlabel("Directional Bias", fontsize=32)
ax.set_ylabel("Normalized Abs Mean Difference\n in FR By Dir",fontsize=32)
ax.set_title("Behavioral Bias in Direction \n vs Directional Tuning",fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)

ax = fig.axes[1]
ax.scatter(speeds,diffs,s=_msize,alpha=_a,c=_c)
ax.set_xlabel("Average running Speed (cm/s)",fontsize=32)
ax.set_xlim([sbins[0],sbins[-1]])
ax.set_ylabel("Normalized Abs Mean Difference\n in FR by Dir",fontsize=32)
ax.set_title("Running Speed vs Directional Tuning",fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)

ax = fig.axes[2]
ax.scatter(speeds,biases,s=_msize,alpha=_a,c=_c)
ax.set_xlabel("Average running speed (cm/s)",fontsize=32)
ax.set_xlim([sbins[0],sbins[-1]])
ax.set_ylabel("Directional Bias",fontsize=32)
ax.set_title("Running Speed vs Directional Bias", fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)


# ax = fig.axes[3]
# h,_,_ = np.histogram2d(diffs,biases,bins=[dbins,bbins])
# ax.imshow(h/h.sum(axis=0),
#                         origin='lower',
#                         extent=[bbins[0],bbins[-1],dbins[0],dbins[-1]],
#                         aspect='auto',
#                         vmax=.3)
# ax.set_xlabel("Behavioral Direction Bias",fontsize=18)
# ax.set_ylabel("Absolute Mean Difference in Directions, \n normalized by field mean (%)",fontsize=18)
# ax.set_title("row-wise normed")



# ax = fig.axes[4]
# h,_,_ = np.histogram2d(diffs,speeds,bins=[dbins,sbins])
# ax.imshow(h,
#                         origin='lower',
#                         extent=[sbins[0],sbins[-1],dbins[0],dbins[-1]],
#                         aspect='auto')
# ax.set_xlabel("Average Running Speed (cm/s)",fontsize=18)
# ax.set_ylabel("Absolute Mean Difference in Directions, \n normalized by field mean (%)",fontsize=18)


# ax = fig.axes[5]
# h,_,_ = np.histogram2d(speeds,biases,bins=[sbins,bbins])
# ax.imshow(h/h.sum(axis=0),
#                         origin='lower',
#                         extent=[bbins[0],bbins[-1],sbins[0],sbins[-1]],
#                         aspect='auto'
#                         )
# ax.set_xlabel("Behavioral Direction Bias",fontsize=18)
# ax.set_ylabel("Average Running Speed",fontsize=18)
# ax.set_title("row-wise normed")

#%% troubleshooting weird envelope thing

sortedpairs = sorted(pairs,key=lambda x:x[2])
fig, ax = plt.subplots()
labels = []
ax2 = ax.twinx()
for i,p in enumerate(sortedpairs):
    a,b = p[4], p[5]
    ax.scatter([i]*a.shape[0],a.Rate,color='r')
    ax.scatter([i]*b.shape[0],b.Rate,color='b')
    ax2.plot(i,p[3],color='green',marker='^',markersize=10,alpha=0.5)
    labels.append(f"{p[0]}_{round(p[1],2)}_{round(p[2],2)}")
#ax.vlines(range(len(pairs)),0,max([max(p[4].Rate.max(),p[5].Rate.max()) for p in pairs]),
#          color='k',
#         alpha=0.5)
#_=ax.set_xticks(range(len(labels)))
#_=ax.set_xticklabels(labels,rotation=90)

