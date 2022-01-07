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
import ratterdam_Defaults as Def

#%% Example Decoder performance
# 21-09-07_decoding. R859D2 CurrentDirection RS7.

# this wont run exactly as is, but most of the code to regen the figure panel 

runfile('E:/UserData/Documents/GitHub/ratterdam/RatterdamOpen_Project/repetition_dummyClassifier_Direction.py', wdir='E:/UserData/Documents/GitHub/ratterdam/RatterdamOpen_Project')

import json
with open("E:\\Ratterdam\\repetition_decoding\\21-09-07_decoding\\rfdata.json","rb") as f:
    realdata = json.load(f)

real = realdata[11]
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
fig.axes[0].hist(real[5],bins=np.linspace(0,1,100),facecolor='b lack', edgecolor='grey')
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
from scipy.stats import mannwhitneyu
#datapath = "E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv"
datapath  = "E:\\Ratterdam\\R_data_repetition\\211210_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

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

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(5)
ax.spines['bottom'].set_linewidth(5)
ax.set_ylabel("Absolute Mean Difference (Hz)",fontsize=Def.ylabelsize)
ax.set_xlabel("Field Number",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':44})
#change the marker size manually for both lines
lgnd.legendHandles[0]._legmarker.set_markersize(30)
lgnd.legendHandles[1]._legmarker.set_markersize(30)
lgnd.get_frame().set_linewidth(4)

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
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Density of Fields",fontsize=Def.ylabelsize)
ax.set_xlabel("Absolute Mean Difference (Hz)",fontsize=Def.xlabelsize)
plt.legend(prop={'size':30})


#%% Fig 3e OOB score vs OOB shuffle and Naive classifier across datasets

# DEPRECATED - keeping for record. even though this is on github. 

# with open("E:\\Ratterdam\\repetition_decoding\\21-09-07_decoding\\rfdata.json","rb") as f:
#     rfdata = json.load(f)
    
# oobs = [i[4] for i in rfdata]
# null95 = [np.percentile(i[5],95) for i in rfdata]
# null5 = [np.percentile(i[5],5) for i in rfdata]


# # these numbers were from a run did earlier in the poster making process
# # and I only saved these 95%iles. But I want final poster to have CI bands
# # so needed to rerun and save everything so these values are obsolete and not used 

# # naive = {'R781D3': {'RS6': 0.5337954939341422, 'RS7': 0.5367612293144207},
# #  'R781D4': {'RS6': 0.5310015898251192, 'RS7': 0.5694444444444444},
# #  'R808D6': {'RS6': 0.6726057906458798, 'RS7': 0.6341463414634146},
# #  'R808D7': {'RS6': 0.7523105360443623, 'RS7': 0.7068965517241379},
# #  'R859D1': {'RS6': 0.84375, 'RS7': 0.8634361233480177},
# #  'R859D2': {'RS6': 0.7798036465638148, 'RS7': 0.7804428044280443},
# #  'R886D1': {'RS6': 0.5480769230769231, 'RS7': 0.6162790697674418},
# #  'R886D2': {'RS6': 0.5841220423412203, 'RS7': 0.5549242424242424},
# #  'R765RFD5': {'RS6': 0.5707434052757794, 'RS7': 0.5788043478260869}}
# naives = []
# labels = [] # need H,V labels for plot below, use here to get them
# for rat,day in zip(['R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886', 'R765'],['D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2','RFD5']):
#     for rslabel in ['RS6','RS7']:
#         naives.append(naive[f"{rat}{day}"][rslabel])
#         if rslabel =='RS6':
#             labels.append(f"{rat}{day} H")
#         elif rslabel == 'RS7':
#             labels.append(f"{rat}{day} V")
        
        
# fig, _ax = plt.subplots()
# ax = fig.axes[0]
# ax.plot(oobs,color='r',marker='o',markersize=20,label='Real OOB Score')
# ax.plot(null95,color='k',marker='o',markersize=20,label='95th Percentile Shuffled OOB Score')
# ax.plot(naives,color='dodgerblue',marker='o',markersize=20,label='95th Percentile Naive Classifier Score')
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.set_ylabel("Decoder Performance",fontsize=24)
# ax.set_xlabel("Dataset",fontsize=24)
# ax.set_xticks(range(len(labels)))
# ax.set_xticklabels(labels,rotation=90)
# plt.legend(prop={'size':20})

#%% fig 3e revised 

oobpct5 = [np.percentile(d[5],5) for d in rfdata]
oobpct95 = [np.percentile(d[5],95) for d in rfdata]
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
ax.plot(real,color='r',linewidth=3,markersize=20,label="Real OOB Score")
ax.fill_between(range(len(oobpct5)), oobpct5, oobpct95, 
                edgecolor='black',
                facecolor='grey',
                alpha=0.8, 
                label="Shuffled OOB Score")
ax.fill_between(range(len(naivepct5)), naivepct5, naivepct95,
                edgecolor='navy',
                facecolor='dodgerblue',
                alpha=0.8,
                label="Naive Classifier")
ax.hlines(0.5,0,18,linestyle='--',color='k',label="Nominal Chance level")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel("Decoder Performance",fontsize=Def.ylabelsize)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize-10)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=45)
plt.legend(prop={'size':30})
plt.subplots_adjust(bottom=0.3)


#%% 3f - model comparison - lrt base vs base+CD

cdmodel = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure2\\211216_CDmodel.csv")

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
ax.set_ylabel("Base Model \n+ Current Direction RMSE",fontsize=Def.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=Def.xlabelsize)
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
plt.legend(prop={'size':30})
