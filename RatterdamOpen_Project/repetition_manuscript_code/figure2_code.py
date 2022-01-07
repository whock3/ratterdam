# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 14:54:47 2021

@author: whockei1

Figure 2 code - Directionality across whole population
"""

import pickle, numpy as np, os, sys, pandas as pd, json, pickle

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator


from scipy.stats import binom_test
from scipy.stats import mannwhitneyu
from scipy.stats import sem 

import repetition_manuscript_defaults as MDef
import ratterdam_Defaults as Def

datapath  = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

#%% Figure 3A - directional ratemaps. Taken from saved files, no py code used here

#%% Figure 3B - scatterplot of field directionality, colored by MW test result

pvals = []
meanDiff = []
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)

        dirA = fgroup[fgroup.CurrDir==dirs[0]]
        dirB = fgroup[fgroup.CurrDir==dirs[1]]
        try:
            pvals.append(mannwhitneyu(dirA.Rate, dirB.Rate).pvalue)
            meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
        except:
            pvals.append(1)
            meanDiff.append(0)
            
            

meanDiff = np.asarray(meanDiff)
pvals = np.asarray(pvals)
fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.plot(np.where(pvals>=0.05)[0],meanDiff[pvals>=0.05],
        linestyle='',
        marker='.',
        color='k',
        markersize=30,
        label='Non-directional')
ax.plot(np.where(pvals<0.05)[0],meanDiff[pvals<0.05],
        linestyle='',
        marker='.',
        color='r',
        markersize=30,
        label='Directional')

ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=MDef.ylabelsize)
ax.set_xlabel("Field Number",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Fig 3C - distribution of 3B visualized as overlaid histograms
# Need to run 3B code first, as this is another way of visualizing that

fig, _ax = plt.subplots()
ax = fig.axes[0]
c,d,e = ax.hist([meanDiff[pvals>=0.05],meanDiff[pvals<0.05]],
        stacked=True,
        density=True,
        color=['grey','lightcoral'],
        bins=np.linspace(0,12,12),
        label=['Non-directional Fields','Directional Fields'])
plt.setp(e[0],edgecolor='k',linewidth=2)
plt.setp(e[1],edgecolor='maroon',linewidth=2)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Density of Fields",fontsize=MDef.ylabelsize)
ax.set_xlabel("Absolute Mean Difference (Hz)",fontsize=MDef.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    #for some reason, lgnd.legendHandles._legmarker doesnt work here? But does above?
    lhand._sizes = [MDef.legend_marker_size]
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Fig 3D - GLM LRT analysis for CD, whole pop

cdmodel = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure2\\211216_CDmodel.csv")

fig, ax = plt.subplots()
ax.plot(cdmodel.m1_rmse[cdmodel.sigP==0],cdmodel.m2_rmse[cdmodel.sigP==0],
        marker='o',
        markersize=20,
        linestyle='',
        color='black',
        label='Non-directional')
ax.plot(cdmodel.m1_rmse[cdmodel.sigP==1],cdmodel.m2_rmse[cdmodel.sigP==1],
        marker='o',
        markersize=20,
        linestyle='',
        color='red',
        label='Directional')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=MDef.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=MDef.xlabelsize)
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Fig 3E - Example classifier run. Example dataset: 21-09-07_decoding. R859D2 CurrentDirection RS7.
# Did not save naive classifier, need to rerun.

# So either manually format this one panel or regen naive perfs
# then edit saved code from sfn 2021 poster to regen panel.

#%% Fig 3F - RF performance across datasets

