# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 18:11:05 2021

@author: whockei1
Code to generate figure panels for SfN 2021 and repetition manuscript first draft
Figure 3 code - directionality for whole population
"""
#%% Imports, applies to all panel code blocks. run once at beginning.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
from scipy.stats import mannwhitneyu
from scipy.stats import sem 

#%% Example Decoder performance
# 21-09-07_decoding. R859D2 CurrentDirection RS7.

# this wont run exactly as is, but most of the code to regen the figure panel 

runfile('E:/UserData/Documents/GitHub/ratterdam/RatterdamOpen_Project/repetition_dummyClassifier_Direction.py', wdir='E:/UserData/Documents/GitHub/ratterdam/RatterdamOpen_Project')

import json
with open("E:\\Ratterdam\\repetition_decoding\\21-09-07_decoding\\rfdata.json","rb") as f:
    realdata = json.load(f)

real = realdata[11]
plt.hist(naive_perf,bins=np.linspace(0,100,20),facecolor='dodgerblue', edgecolor='navy')
plt.hist(naive_perf,bins=np.linspace(0,1,20),facecolor='dodgerblue', edgecolor='navy')
plt.hist(naive_perf,bins=np.linspace(0,1,100),facecolor='dodgerblue', edgecolor='navy')
plt.hist(naive_perf,bins=np.linspace(0,1,50),facecolor='dodgerblue', edgecolor='navy')
plt.hist(naive_perf,bins=np.linspace(0,1,100),facecolor='dodgerblue', edgecolor='navy')
plt.hist(real[5],bins=np.linspace(0,1,100),facecolor='black', edgecolor='grey')
plt.vlines(real[4],0,250,color='r')
plt.vlines(real[4],0,250,color='r',linewidth=2)
plt.vlines(real[4],0,250,color='r',linewidth=3)
plt.xticks(fontsize=16)
plt.xticklabels(fontsize=16)
plt.xlabels(fontsize=16)
fig, ax = plt.subplots()
fig.axes.hist(naive_perf,bins=np.linspace(0,1,100),facecolor='dodgerblue', edgecolor='navy')
fig.axes[0].hist(naive_perf,bins=np.linspace(0,1,100),facecolor='dodgerblue', edgecolor='navy')
fig.axes[0].hist(real[5],bins=np.linspace(0,1,100),facecolor='black', edgecolor='grey')
fig.axes[0].vlines(real[4],0,250,color='r',linewidth=3)
fig.axes[0].tick_params(axis='both', which='major', labelsize=18)
fig.axes[0].tick_params(axis='both', which='major', labelsize=20)
fig.axes[0].set_ylabel("Frequency",fontsize=24)
fig.axes[0].set_xlabel("Classifier Performance",fontsize=24)
fig.axes[0].spikes['top'].set_visible(False)
fig.axes[0].spines['top'].set_visible(False)
fig.axes[0].spines['right'].set_visible(False)
fig.axes[0].set_xlabel("Classifier Performance",fontsize=30)
fig.axes[0].set_ylabel("Frequency",fontsize=30)


#%% Fig 3b and 3c - scatterplot of each field's directional tuning and hist of same
# v

datapath = "E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv"
df = pd.read_csv(datapath)

passThresh = 5

pvals = []
meanDiff = []
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)
        # there are two unique directions through an alley. If there are fewer,
        # then dont run bc not enough sampling. If more, there's an error
        if len(dirs) == 2:
            dirA = fgroup[fgroup.CurrDir==dirs[0]]
            dirB = fgroup[fgroup.CurrDir==dirs[1]]
            if dirA.shape[0] >= passThresh and dirB.shape[0] >= passThresh:
                pvals.append(mannwhitneyu(dirA.Rate, dirB.Rate).pvalue)
                meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            
            

meanDiff = np.asarray(meanDiff)
pvals = np.asarray(pvals)
fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.plot(meanDiff[pvals>=0.05],linestyle='',marker='.',color='k',markersize=15,label='Non-Directional')
ax.plot(meanDiff[pvals<0.05],linestyle='',marker='.',color='r',markersize=15,label='Directional')
ax.tick_params(axis='both', which='major', labelsize=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=24)
ax.set_xlabel("Field Number",fontsize=24)
plt.legend(prop={'size':20})

fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.hist(meanDiff,bins=20,color='grey',edgecolor='black')
ax.tick_params(axis='both', which='major', labelsize=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Number of Fields",fontsize=24)
ax.set_xlabel("Absolute Mean Difference (Hz)",fontsize=24)


#%% Fig 3e OOB score vs OOB shuffle and Naive classifier across datasets
with open("E:\\Ratterdam\\repetition_decoding\\21-09-07_decoding\\rfdata.json","rb") as f:
    rfdata = json.load(f)
    
oobs = [i[4] for i in rfdata]
null95 = [np.percentile(i[5],95) for i in rfdata]


naive = {'R781D3': {'RS6': 0.5337954939341422, 'RS7': 0.5367612293144207},
 'R781D4': {'RS6': 0.5310015898251192, 'RS7': 0.5694444444444444},
 'R808D6': {'RS6': 0.6726057906458798, 'RS7': 0.6341463414634146},
 'R808D7': {'RS6': 0.7523105360443623, 'RS7': 0.7068965517241379},
 'R859D1': {'RS6': 0.84375, 'RS7': 0.8634361233480177},
 'R859D2': {'RS6': 0.7798036465638148, 'RS7': 0.7804428044280443},
 'R886D1': {'RS6': 0.5480769230769231, 'RS7': 0.6162790697674418},
 'R886D2': {'RS6': 0.5841220423412203, 'RS7': 0.5549242424242424},
 'R765RFD5': {'RS6': 0.5707434052757794, 'RS7': 0.5788043478260869}}
naives = []
labels = [] # need H,V labels for plot below, use here to get them
for rat,day in zip(['R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886', 'R765'],['D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2','RFD5']):
    for rslabel in ['RS6','RS7']:
        naives.append(naive[f"{rat}{day}"][rslabel])
        if rslabel =='RS6':
            labels.append(f"{rat}{day} H")
        elif rslabel == 'RS7':
            labels.append(f"{rat}{day} V")
        
        
fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.plot(oobs,color='r',marker='o',markersize=10,label='Real OOB Score')
ax.plot(null95,color='k',marker='o',markersize=10,label='95th Percentile Shuffled OOB Score')
ax.plot(naives,color='dodgerblue',marker='o',markersize=10,label='95th Percentile Naive Classifier Score')
ax.tick_params(axis='both', which='major', labelsize=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Decoder Performance",fontsize=24)
ax.set_xlabel("Dataset",fontsize=24)
ax.set_xticks(range(len(oobs)))
plt.legend(prop={'size':20})
