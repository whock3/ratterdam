# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 15:29:12 2021

@author: whockei1

Script to generate the panels for Figure 4 of repetition manuscript
Each panel in the paper corresponds to a code cell here. Code cells in python
begin with "#%%". They can be run independently. The whole script can be run to generate all panels.

Figure 4 compares the directionality tuning between repeating and nonrepeating neurons
A - violin plots showing distribution of normalized directionality scores for each group
B - scatterplot of GLM results. Each axis is RMSE of indicated model. Point coloring is whether
    the time model was associated with a significant improvement in fit as assessed by LRT
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
from scipy.stats import mannwhitneyu
from scipy.stats import sem 
import ratterdam_Defaults as Def 
import repetition_manuscript_defaults as MDef


datapath  = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)
plt.ion()
#%% Figure 4A - Violins of Repeating vs Non-repeating directionality (abs mean diff)

meanDiff = []
repeating = []
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)
        if len(dirs)==2:
            dirA = fgroup[fgroup.CurrDir==dirs[0]]
            dirB = fgroup[fgroup.CurrDir==dirs[1]]
            meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            repeating.append(np.unique(fgroup.Repeating)[0])
                    
meanDiff = np.asarray(meanDiff)
repeating = np.asarray(repeating)
                
fig, _ax = plt.subplots()

ax = fig.axes[0]
repviolin = ax.violinplot(meanDiff[repeating==True],[1])
repviolin['bodies'][0].set_facecolor('gray')
repviolin['bodies'][0].set_edgecolor('black')
repviolin['bodies'][0].set_linewidth(3)
repviolin['bodies'][0].set_alpha(0.8)
for el in ['cbars','cmaxes','cmins']:
    repviolin[el].set_color('black')
    
nonrepviolin = ax.violinplot(meanDiff[repeating==False],[2])
nonrepviolin['bodies'][0].set_facecolor('grey')
nonrepviolin['bodies'][0].set_edgecolor('black')
nonrepviolin['bodies'][0].set_linewidth(3)
nonrepviolin['bodies'][0].set_alpha(0.8)
for el in ['cbars','cmaxes','cmins']:
    nonrepviolin[el].set_color('black')

ax.set_xticks([1,2])
ax.set_xticklabels(["Repeating Fields", "Non-repeating Fields"],fontsize=MDef.xlabelsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Normalized Mean Difference (Hz)",fontsize=MDef.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)


#%% Figure 4B - GLM LRT CD, rep vs nonrep 

cdmodel = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure3_Directionality\\20220414_CDmodel.csv")

fig, ax = plt.subplots()

# nonrep sig 
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==1)],cdmodel.m2_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==1)],
        marker='o',
        markersize=20,
        linestyle='',
        color='red',
        alpha=0.6,
        label='Non-repeating Significant')

# rep sig
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==1)],cdmodel.m2_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==1)],
        marker='^',
        markersize=20,
        linestyle='',
        color='red',
        alpha=0.6,
        label='Repeating Significant')

# nonrep nonsig
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==0)],cdmodel.m2_rmse[(cdmodel.repOrNot==0)&(cdmodel.sigP==0)],
        marker='o',
        markersize=20,
        linestyle='',
        color='black',
        alpha=0.6,
        label='Non-repeating Non-significant')

# rep nonsig
ax.plot(cdmodel.m1_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==0)],cdmodel.m2_rmse[(cdmodel.repOrNot==1)&(cdmodel.sigP==0)],
        marker='^',
        markersize=20,
        linestyle='',
        color='black',
        alpha=0.6,
        label='Repeating Non-significant')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=MDef.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=MDef.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})