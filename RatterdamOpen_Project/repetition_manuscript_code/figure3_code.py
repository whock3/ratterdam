# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 15:29:12 2021

@author: whockei1

Figure 3 - Repeating versus Non-repeating Cells wrt Directionality
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
from scipy.stats import mannwhitneyu
from scipy.stats import sem 
import ratterdam_Defaults as Def 
import repetition_manuscript_defaults as MDef


datapath  = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

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
ax.set_xticklabels(["Repeating Fields", "Non-repeating Fields"],fontsize=MDef.xlabelsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=MDef.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)


#%% Figure 4B - GLM LRT CD, rep vs nonrep 

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
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=MDef.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=MDef.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)


#%% Figure 3C - schematic of GLM LRT tests. No py code used here

#%% Fig 3D - results of GLM LRTs 

#taken from R results 
rep = [31/128, 18/128, 9/128, 2/128]
nr_all = [18/75, 6/75, 6/75, 4/75]
nr_multi = [12/45, 6/45, 4/45, 2/45]
nr_single = [6/30, 0/30, 2/30, 2/30]

rows = ["+Current Dir", "+Previous Dir", "+Next Dir", "All Dirs"]

data = {'Repeating':rep,
        "Multi-field Non-repeating":nr_multi,
        "Single field Non-repeating":nr_single
       }

df = pd.DataFrame(data)

ax = df.plot(kind='bar',
        color=['red','navy','cornflowerblue','lightblue'],
        fontsize=MDef.xlabelsize
        )

ax.set_xticks([0,1,2,3])
ax.set_xticklabels(rows,rotation=0)
ax.set_ylabel("Proportion of Fields",fontsize=MDef.ylabelsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.hlines(0.05,-0.5,4,linestyle='--',color='k',linewidth=5)
plt.subplots_adjust(bottom=0.2)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    #for some reason, lgnd.legendHandles._legmarker doesnt work here? But does elsewhere?
    lhand._sizes = [MDef.legend_marker_size]
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)