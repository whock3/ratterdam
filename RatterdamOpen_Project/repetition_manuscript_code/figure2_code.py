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

#Using 21-09-07 decoding run 

# graph format was line graph (or shaded linear regions). But this evokes a 
# temporal trend so using unconnected points / error bars. 2-11-22

with open("E:\\Ratterdam\\repetition_decoding\\21-09-07_decoding\\rfdata.json","rb") as f:
    rfdata_badorder = json.load(f)
    
rfdata = rfdata_badorder[-2:] + rfdata_badorder[:-2] # the way it's saved has R765 at the end, so move it in front so its in rat order

with open("E:\\Ratterdam\\2021_SfNPoster_WH\\Fig3_Directionality\\naivePerfs.json","rb") as f:
    naive_perfs_datasets = json.load(f)

oobpct5 = [np.nanpercentile(d[5],5) for d in rfdata]
oobpct95 = [np.nanpercentile(d[5],95) for d in rfdata]
real = [d[4] for d in rfdata]

naivepct5 = []
naivepct95 = []
labels = []
#naive perfs are in dict, rfdata in listoflists. to make sure youre grabbing in right order for both:
ratorder,dayorder,regionorder = [d[0] for d in rfdata], [d[1] for d in rfdata], [d[3] for d in rfdata]
for r,d,l in zip(ratorder, dayorder, regionorder):
    naive = naive_perfs_datasets[f"{r}{d}"][l]
    naivepct5.append(np.percentile(naive,5))
    naivepct95.append(np.percentile(naive,95))
    if l == 'RS6':
        labels.append(f"{r}{d} H")
    elif l == 'RS7':
        labels.append(f"{r}{d} V")
    
    

fig, ax = plt.subplots()
ax.scatter(range(len(real)), real,
           color='r',
           marker='^',
           s=200,
           edgecolor='k',
           linewidth=0.5,
           label="Real OOB Score")
# below is the midpoint between the 5th and 95th percentiles. Has no scientific use here, just need it for plotting errorbars
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
ax.set_ylabel("Decoder Performance",fontsize=Def.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize-10)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=90)
plt.legend(prop={'size':30})
plt.subplots_adjust(bottom=0.3)
