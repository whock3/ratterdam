# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 18:26:47 2021

@author: whockei1

Fig 4  - Directionality in Repeating vs Nonrepeating Neurons
for Sfn 2021 and first draft paper 
"""
#%% Imports and loading data

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
from scipy.stats import mannwhitneyu
from scipy.stats import sem 
import ratterdam_Defaults as Def 

#datapath = "E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv"
datapath  = "E:\\Ratterdam\\R_data_repetition\\211210_AlleySuperpopDirVisitFiltered.csv"

df = pd.read_csv(datapath)
df = df[df.Traversal==True]

#%% Fig 4a, 4b- violins of repeating vs nonrepeating, then scatter same
meanDiff = []
repeating = []
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)
        # there are two unique directions through an alley. If there are fewer,
        # then dont run bc not enough sampling. If more, there's an error
        if len(dirs) == 2:
            dirA = fgroup[fgroup.CurrDir==dirs[0]]
            dirB = fgroup[fgroup.CurrDir==dirs[1]]
            meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            repeating.append(np.unique(fgroup.Repeating)[0])
                
meanDiff = np.asarray(meanDiff)
repeating = np.asarray(repeating)
                
fig, _ax = plt.subplots()

ax = fig.axes[0]
repviolin = ax.violinplot(meanDiff[repeating==True],[1])
repviolin['bodies'][0].set_facecolor('red')
repviolin['bodies'][0].set_edgecolor('darkred')
repviolin['bodies'][0].set_linewidth(2)
for el in ['cbars','cmaxes','cmins']:
    repviolin[el].set_color('darkred')
    
nonrepviolin = ax.violinplot(meanDiff[repeating==False],[2])
nonrepviolin['bodies'][0].set_facecolor('dodgerblue')
nonrepviolin['bodies'][0].set_edgecolor('navy')
nonrepviolin['bodies'][0].set_linewidth(2)
for el in ['cbars','cmaxes','cmins']:
    nonrepviolin[el].set_color('dodgerblue')

ax.set_xticks([1,2])
ax.set_xticklabels(["Repeating Fields", "Non-repeating Fields"],fontsize=Def.xlabelsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=Def.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)

fig, _ax  = plt.subplots()
ax = fig.axes[0]
ax.scatter(range(meanDiff[repeating==True].shape[0]),meanDiff[repeating==True],c='red',s=40,zorder=99,label='Repeating')
ax.scatter(range(meanDiff[repeating==False].shape[0]),meanDiff[repeating==False],c='dodgerblue',s=40,zorder=99,label='Non-repeating')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=Def.ylabelsize)
ax.set_xlabel("Field ID",fontsize=Def.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
plt.legend(prop={'size':30})


#%% fig 4c glm results from r, mse change in base glm vs base+CD

cdmodel = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure2\\211216_CDmodel.csv")

fig, ax = plt.subplots()
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==1)],cdmodel.m2_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==1)],
        marker='^',
        markersize=20,
        linestyle='',
        color='navy',
        alpha=0.6,
        label='Non-repeating Significant')
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==1)],cdmodel.m2_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==1)],
        marker='^',
        markersize=20,
        linestyle='',
        color='red',
        alpha=0.6,
        label='Repeating Significant')

ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==0)],cdmodel.m2_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==0)],
        marker='o',
        markersize=20,
        linestyle='',
        color='navy',
        alpha=0.6,
        label='Non-repeating Non-significant')
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==0)],cdmodel.m2_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==0)],
        marker='o',
        markersize=20,
        linestyle='',
        color='red',
        alpha=0.6,
        label='Repeating Non-significant')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=Def.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=Def.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
lgnd = plt.legend(prop={'size':40},loc='upper right')
for leg_handle in lgnd.legendHandles:
    leg_handle._legmarker.set_markersize(30)
ax.set_aspect('equal', adjustable='box')

# to do
# rep vs non scatter different color for rep status and marker for sig
# redo tuning x bias with rep or non
