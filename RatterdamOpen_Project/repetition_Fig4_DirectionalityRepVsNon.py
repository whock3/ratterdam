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

datapath = "E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv"
df = pd.read_csv(datapath)
passThresh = 5


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
            if dirA.shape[0] >= passThresh and dirB.shape[0] >= passThresh:
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
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=Def.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)

fig, _ax  = plt.subplots()
ax = fig.axes[0]
ax.scatter(range(meanDiff[repeating==True].shape[0]),meanDiff[repeating==True],c='red',s=40,zorder=99,label='Repeating')
ax.scatter(range(meanDiff[repeating==False].shape[0]),meanDiff[repeating==False],c='dodgerblue',s=40,zorder=99,label='Non-repeating')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=Def.ylabelsize)
ax.set_xlabel("Field ID",fontsize=Def.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
plt.legend(prop={'size':30})


#%% fig 4c glm results from r, mse change in base glm vs base+CD

cdmodel = pd.read_csv("E:\\Ratterdam\\2021_SfNPoster_WH\\Fig3_Directionality\\CD_model.csv")

fig, ax = plt.subplots()
ax.plot(cdmodel.m1_rmse[cdmodel.repOrNot==0],cdmodel.m2_rmse[cdmodel.repOrNot==0],
        marker='o',
        markersize=10,
        linestyle='',
        color='black',
        label='Non-repeating')
ax.plot(cdmodel.m1_rmse[cdmodel.repOrNot==1],cdmodel.m2_rmse[cdmodel.repOrNot==1],
        marker='o',
        markersize=10,
        linestyle='',
        color='red',
        label='Repeating')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=Def.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=Def.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
plt.legend(prop={'size':30})
ax.set_aspect('equal', adjustable='box')

