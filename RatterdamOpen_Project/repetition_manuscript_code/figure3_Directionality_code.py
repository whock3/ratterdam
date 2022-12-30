# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 14:54:47 2021

@author: whockei1

Script to generate the panels for Figure 3 of repetition manuscript
Each panel in the paper corresponds to a code cell here. Code cells in python
begin with "#%%". They can be run independently. The whole script can be run to generate all panels.

Panels:
A - ratemaps from two example neurons with directional firing fields. 
B - Distribution of normalized directional index across population, colored by whether it passes a Mann Whitney
    test 
C - Distribution of GLM RMSEs for base model versus model including direction. Significance coloring corresponds to
    outcome of LRT. 
D - Example random forest classifier run from one recoridng 
E - Random forest results for all datasets
"""

#%% Imports and load data

import pickle, numpy as np, os, sys, pandas as pd, json, pickle
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import binom_test
from scipy.stats import mannwhitneyu
from scipy.stats import sem 
import repetition_manuscript_defaults as MDef
import ratterdam_Defaults as Def

datapath = '' # Add your datapath here

df = pd.read_csv(datapath+"2022-04-05_AlleySuperpopDirVisitFiltered.csv")

plt.ion()

#%% Figure 3A - directional ratemaps. No py code used here

#%% Figure 3B - scatterplot of field directionality, colored by MW test result
pvals = []
meanDiff = []
mw_fids = [] # field IDs from MW test
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)

        dirA = fgroup[fgroup.CurrDir==dirs[0]]
        dirB = fgroup[fgroup.CurrDir==dirs[1]]
        try:
            pvals.append(mannwhitneyu(dirA.Rate, dirB.Rate).pvalue)
            meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            mw_fids.append(fname)
        except:
            # pvals.append(1)
            # meanDiff.append(0)
            pass
alpha = 0.05/2            

meanDiff = np.asarray(meanDiff)
pvals = np.asarray(pvals)
mw_fids = np.asarray(mw_fids)

fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.plot(np.where(pvals>=alpha)[0],meanDiff[pvals>=alpha],
        linestyle='',
        marker='.',
        color='k',
        markersize=30,
        label='Non-directional')
ax.plot(np.where(pvals<alpha)[0],meanDiff[pvals<alpha],
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
ax.set_ylabel("Normalized Mean Difference (Hz)",fontsize=MDef.ylabelsize)
ax.set_xlabel("Field Number",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)
ax.set_ylim([0,7.5])


#%% Fig 3C - GLM LRT analysis for CD, whole pop
cdmodel = pd.read_csv(datapath+"2022-05-25_CDmodel.csv")

fig, ax = plt.subplots()
ax.plot(cdmodel.m1_rmse[cdmodel.sigP==0],cdmodel.m2_rmse[cdmodel.sigP==0],
        marker='o',
        markersize=10,
        linestyle='',
        color='black',
        label='Non-directional')
ax.plot(cdmodel.m1_rmse[cdmodel.sigP==1],cdmodel.m2_rmse[cdmodel.sigP==1],
        marker='o',
        markersize=10,
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
ax.set_xlim([0.9,3.9])
ax.set_ylim([0.9,3.9])
ax.set_aspect(1./ax.get_data_ratio())

#%% Fig 3D - Example classifier run. Example dataset: 22-02-21. R765 DFD4 horizontal alleys

with open(datapath+"rfdata.json", "r") as f:
    rf_data = json.load(f)

with open(datapath+"naiveClassifierData.json", "r") as f:
    naive_data = json.load(f)
    
real = rf_data['R765']['DFD4']['RS6']['oobs']['Real']
shuffoob = rf_data['R765']['DFD4']['RS6']['oobs']['Shuffle']
naive = naive_data['R765DFD4']['RS6']
bins = np.linspace(0.2,0.75,35)
fig, ax = plt.subplots()
plt.hist(shuffoob,bins=bins,color='grey',linewidth=2,label='OOB Shuffle')
plt.hist(naive,bins=bins,color='blue',linewidth=2,label='Naive Classifier')
plt.vlines(real,0,200,color='r',linewidth=3,label='Real OOB Score')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Frequency",fontsize=MDef.ylabelsize)
ax.set_xlabel("Decoder Accuracy",fontsize=MDef.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Fig 3E - RF performance across datasets

# graph format was line graph (or shaded linear regions). But this evokes a 
# temporal trend so using unconnected points / error bars. 2-11-22

    
rat_list = ['R765',
            'R765',
            'R781', 
            'R781', 
            'R808', 
            'R808', 
            'R859', 
            'R859', 
            'R886', 
            'R886']

day_list = ['RFD5',
            'DFD4',
            'D3', 
            'D4',
            'D6',
            'D7',
            'D1',
            'D2',
            'D1',
            'D2']

real, oobpct5, oobpct95, naivepct5, naivepct95 = [], [], [], [], []
labels = []
for rat, day in zip(rat_list, day_list):
    
    # RS6 = horizontal alleys, RS7 = vertical alleys
    for regionset,label in zip(['RS6','RS7'],['H','V']):
        labels.append(f"{rat}{day} {label}")
        
        real.append(rf_data[rat][day][regionset]['oobs']['Real'])
        oobpct5.append(np.percentile(rf_data[rat][day][regionset]['oobs']['Shuffle'],2.5))
        oobpct95.append(np.percentile(rf_data[rat][day][regionset]['oobs']['Shuffle'],97.5))
        naivepct5.append(np.percentile(naive_data[f"{rat}{day}"][regionset],2.5))
        naivepct95.append(np.percentile(naive_data[f"{rat}{day}"][regionset],97.5))


fig, ax = plt.subplots()
ax.scatter(range(len(real)), real,
           color='r',
           marker='^',
           s=200,
           edgecolor='k',
           linewidth=0.5,
           label="Real OOB Score")
# below is the midpoint between the 5th and 95th percentiles. Has no scientific meaning, just need it for plotting errorbars
oobmids = np.asarray([(a+b)/2 for a,b in zip(oobpct5, oobpct95)])
naivemids = np.asarray([(a+b)/2 for a,b in zip(naivepct5, naivepct95)])
ax.errorbar(range(len(oobpct5)),oobmids,yerr=np.asarray([oobmids-oobpct5,oobpct95-oobmids]),
            ls='none',
            capsize=20,
            capthick=1,
            elinewidth=1,
            color='black',
            label='Shuffled OOB Score'
            )
ax.errorbar(range(len(naivepct5)),naivemids,yerr=np.asarray([naivemids-naivepct5,naivepct95-naivemids]),
            ls='none',
            capsize=20,
            capthick=1,
            elinewidth=1,
            color='blue',
            label='Naive Classifier'
            )

ax.hlines(0.5,0,18,linestyle='--',color='k',label="Nominal Chance level")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Decoder Performance",fontsize=MDef.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=90,fontsize=MDef.ticksize)
plt.legend(prop={'size':30})
plt.subplots_adjust(bottom=0.3)
# %%
