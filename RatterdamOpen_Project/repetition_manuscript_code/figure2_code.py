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

# most recent version using raw FR (without R765 DFD4): 211220_*

df = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv")


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
            # pvals.append(1)
            # meanDiff.append(0)
            pass
            
            

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
ax.set_ylabel("Normalized Mean Difference (Hz)",fontsize=MDef.ylabelsize)
ax.set_xlabel("Field Number",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Stats for M-W test 
from scipy.stats import binom_test
print(f"Total number field segments: {len(pvals)}")
print(f"M-W directional fields: {np.where(pvals<0.05)[0].shape[0]}")
print(binom_test(np.where(pvals<0.05)[0].shape[0],len(pvals),0.05,'greater'))

#%% Fig 3C - distribution of 3B visualized as overlaid histograms
# Need to run 3B code first, as this is another way of visualizing that

fig, _ax = plt.subplots()
ax = fig.axes[0]
c,d,e = ax.hist([meanDiff[pvals>=0.05],meanDiff[pvals<0.05]],
        stacked=True,
        density=True,
        color=['grey','lightcoral'],
        bins=np.linspace(0,9,20),
        label=['Non-directional Fields','Directional Fields'])
plt.setp(e[0],edgecolor='k',linewidth=2)
plt.setp(e[1],edgecolor='maroon',linewidth=2)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Density of Fields",fontsize=MDef.ylabelsize)

#before 220218, units were Hz. On and after this date, we are using normalized
# Firing rate. Normalizing by 2D session-avg RM. 
ax.set_xlabel("Normalized Mean Difference (Hz)",fontsize=MDef.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    #for some reason, lgnd.legendHandles._legmarker doesnt work here? But does above?
    lhand._sizes = [MDef.legend_marker_size]
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

#%% Fig 3D - GLM LRT analysis for CD, whole pop

# this 211216 model did not have its own script, I saved the CD part of a script
# that looked at all directions. 
cdmodel = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure2\\220406_CDmodel.csv")

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

#%% Fig 3E - Example classifier run. Example dataset: 22-02-21. R765 DFD4 RS6

with open("E:\\Ratterdam\\repetition_decoding\\2022-04-11_decoding\\rfdata.json", "r") as f:
    rf_data = json.load(f)

with open("E:\\Ratterdam\\repetition_decoding\\2022-04-11_decoding\\naiveClassifierData.json", "r") as f:
    naive_data = json.load(f)
    
real = rf_data['R765']['DFD4']['RS6']['oobs']['Real']
shuffoob = rf_data['R765']['DFD4']['RS6']['oobs']['Shuffle']
naive = naive_data['R765DFD4']['RS6']
    
fig, ax = plt.subplots()
plt.hist(shuffoob,bins=30,color='grey',linewidth=2,label='OOB Shuffle')
plt.hist(naive,bins=30,color='blue',linewidth=2,label='Naive Classifier')
plt.vlines(real,0,100,color='r',linewidth=3,label='Real OOB Score')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Frequency",fontsize=MDef.ylabelsize)
ax.set_xlabel("Decoder Accuracy",fontsize=MDef.xlabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)

    


#%% Fig 3F - RF performance across datasets

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
        oobpct5.append(np.percentile(rf_data[rat][day][regionset]['oobs']['Shuffle'],5))
        oobpct95.append(np.percentile(rf_data[rat][day][regionset]['oobs']['Shuffle'],95))
        naivepct5.append(np.percentile(naive_data[f"{rat}{day}"][regionset],5))
        naivepct95.append(np.percentile(naive_data[f"{rat}{day}"][regionset],95))


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
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Decoder Performance",fontsize=Def.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize-10)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=90)
plt.legend(prop={'size':30})
plt.subplots_adjust(bottom=0.3)


#%% 3G - tuning with of directional signal.
# Distributions of directionality when FR is normalized by RM dont look too different
# from each other (those MW test says are directional, vs non). Presumably the 
# effect lies in the tightness/variability of the tuning?

pvals = []
meanDiff = []
meanVar = []
meanFanos = []
for orien in ['V','H']:
    odf = df[df.Orientation==orien]
    for fname, fgroup in odf.groupby("FieldID"):
        dirs = np.unique(fgroup.CurrDir)

        dirA = fgroup[fgroup.CurrDir==dirs[0]]
        dirB = fgroup[fgroup.CurrDir==dirs[1]]
        try:
            pvals.extend([mannwhitneyu(dirA.Rate, dirB.Rate).pvalue]*2)
            meanDiff.append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            meanVar.append(np.nanmean([dirA.Rate.sem(),dirB.Rate.sem()]))
            
            fanoA = np.var(dirA.Rate)/np.nanmean(dirA.Rate)
            fanoB = np.var(dirB.Rate)/np.nanmean(dirB.Rate)
            meanFanos.extend([fanoA,fanoB])
        except:
            pass
            # pvals.append(1)
            # meanDiff.append(0)
            # meanVar.append(999)
            
bins = np.linspace(0,np.nanpercentile(meanFanos,90),50)   

meanDiff = np.asarray(meanDiff)
pvals = np.asarray(pvals)
meanVar = np.asarray(meanVar)
meanFanos = np.asarray(meanFanos)
fig, _ax = plt.subplots()
ax = fig.axes[0]
ax.hist(meanFanos[pvals>=0.05],
        bins=bins,
        color='grey',
        edgecolor='black',
        linewidth=2,
        density=True,
        stacked=True,
        alpha=0.7,
        label='Non-directional')
ax.hist(meanFanos[pvals<0.05],
        bins=bins,
        color='r',
        edgecolor='darkred',
        linewidth=2,
        density=True,
        stacked=True,
        alpha=0.7,
        label='Directional')

ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Normed Frequency",fontsize=MDef.ylabelsize)
ax.set_xlabel("Fano Factor",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)


#%% H - non-current direction (obviously this figure is getting packed and we will decide later how to unpack. lay it out first)

# manually coding numbers from 2-22 run.using allleypath <- "E:\\Ratterdam\\R_data_repetition\\220222_AlleySuperpopDirVisitFiltered.csv"
# total_current = 154
# total_previous = 167
# total_next = 163

# current_responsive = 34
# previous_responsive = 19
# next_responsive = 24

# binom_test(34,154,0.05,'greater')
# Out[31]: 2.3543593135393474e-13

# binom_test(19,167,0.05,'greater')
# Out[32]: 0.0007481893386773383

# binom_test(24,163,0.05,'greater')
# Out[33]: 2.2510912699473393e-06


#manually coding using 3/23/2022 run input data to R: "E:\\Ratterdam\\R_data_repetition\\2022-03-23_AlleySuperpopDirVisitFiltered.csv"
# running in repetition_GLMs_emmeans.R

total_current = 127
total_previous = 160
total_next = 151

current_responsive = 32
previous_responsive = 27
next_responsive = 20




fig, ax = plt.subplots()

ax.bar([0,1,2],[current_responsive/total_current,previous_responsive/total_previous,next_responsive/total_next],
       color='cornflowerblue',
       edgecolor='navy',
       linewidth=2)
ax.hlines(0.05,-0.5,2.5,linestyle='--',color='k',linewidth=3)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["Current", "Previous", "Next"],fontsize=MDef.xlabelsize)
ax.set_ylabel("Proportion of Fields Responsive",fontsize=MDef.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)