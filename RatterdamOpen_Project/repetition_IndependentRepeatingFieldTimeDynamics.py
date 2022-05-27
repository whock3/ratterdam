# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 19:23:22 2022

@author: whockei1

Script to analyze independence of repeating fields in terms of their temporal
dyamics. GLMs and related models will be run in R. Correlating time series
and related analyses will be done here

When mature, this will be moved to figure5_code for reptition manuscript code 

"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
from scipy.stats import sem
import repetition_manuscript_defaults as MDef
from importlib import reload
import itertools
from scipy.stats.stats import pearsonr
from scipy.interpolate import PchipInterpolator as pchip
from matplotlib import patches
from scipy.signal import correlate 
import scipy.stats

#%%

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

fieldTimeCorrs = []
x,y = [], []

i = 0

blowup_thresh = 50


for rat, rdf in alleydf.groupby("Rat"):
    for day, ddf in rdf.groupby("Day"):
        start = ddf.StartTimes.min()
        end = ddf.StartTimes.max()
        interp_ts = np.linspace(start,end, 100)

        for orien, oriendf in ddf.groupby("Orientation"):
            for cid,cell in oriendf.groupby("CellID"):
                if np.unique(cell.Repeating)[0] == True:
                    
                    field_timeseries = []
                                                            
                    for fid, field in cell.groupby("FieldID"):
                        
                        startTimes = field.StartTimes
                        rates = field.Rate
                        
                        pcf = pchip(startTimes, rates)
                        pcf_interp = pcf(interp_ts)
                        
                        pcf_interp[np.logical_or((pcf_interp>50),
                                                  (pcf_interp<0))
                        ] = 0 # this is a very rough way of doing it
                        
                        field_timeseries.append(pcf_interp)
                        
                        if i%5 == 0:
                        
                            plt.plot(startTimes,rates,linestyle='-',marker='^')
                            plt.plot(interp_ts,pcf_interp,linestyle='--')
                        
                    combs = itertools.combinations(range(len(field_timeseries)),2)
                    for pair in combs:
                        i,j = pair
                        tsA,tsB = field_timeseries[i], field_timeseries[j]
                        
                        c = correlate(tsA, tsB, mode="full", method="auto")
                        max_ind = np.argmax(np.abs(c))
                        shift = max_ind - (len(tsA)-1)
                        
                        corr = pearsonr(tsA,tsB)[0]
                        fieldTimeCorrs.append(corr)
                        
                        
                    i +=1 
                        
                        
#%%                        
fig, ax = plt.subplots()
ax.hist(fieldTimeCorrs,
        bins=25,
        color='grey',
        edgecolor='black',
        linewidth=2)
ylim = ax.get_ylim()
xlim, ylim = ax.get_xlim(), ax.get_ylim()
widths = [xlim[0],
          xlim[1]-xlim[0]]
  
heights = [ylim[1],
           ylim[1]]  

colors = ['blue','red']
alpha = 0.3
for w,h,c in zip(widths, heights, colors):    
    rect = patches.Rectangle((0,0),w,h,
                             facecolor=c,
                             alpha=alpha,
                             zorder=0
                            )

    ax.add_patch(rect)
    
ax.set_ylabel("Frequency", fontsize=40)
ax.set_xlabel("Pearson Correlation Between Interpolated Time Series", fontsize=40)


ax.tick_params(axis='both', which='major', labelsize=40, length=10)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)


                        
#%% 2022-05-23 Revised analysis to get null distribution of correlations between repeating fields

all_fields = []
cell_field_nums = []

for rat, rdf in alleydf.groupby("Rat"):
    for day, ddf in rdf.groupby("Day"):
        start = ddf.StartTimes.min()
        end = ddf.StartTimes.max()
        interp_ts = np.linspace(start,end, 100)
        for orien, oriendf in ddf.groupby("Orientation"):
            for cid,cell in oriendf.groupby("CellID"):
                if np.unique(cell.Repeating)[0] == True:

                    cell_field_nums.append(len(np.unique(cell.FieldID)))
                                                                                
                    for fid, field in cell.groupby("FieldID"):
                        
                        startTimes = field.StartTimes
                        rates = field.Rate
                        
                        pcf = pchip(startTimes, rates)
                        pcf_interp = pcf(interp_ts)
                        
                        pcf_interp[pcf_interp>blowup_thresh] = 0 # this is a very rough way of doing it
                        
                        all_fields.append(pcf_interp)

all_fields = np.asarray(all_fields)

# %%
nboots = 1000
boot_x, boot_y = [], []
bootskews = []

for i in range(nboots):

    bootstrapCorrs = []

    for cfn in cell_field_nums:

        selected_fields_idx = np.random.choice(range(len(all_fields)), cfn, replace=True)
        selected_fields = all_fields[selected_fields_idx]
        combs = itertools.combinations(range(len(selected_fields)),2)
        for pair in combs:
            i,j = pair
            tsA,tsB = selected_fields[i], selected_fields[j]
            
            # c = correlate(tsA, tsB, mode="full", method="auto")
            # max_ind = np.argmax(np.abs(c))
            # shift = max_ind - (len(tsA)-1)
            
            corr = pearsonr(tsA,tsB)[0]
            bootstrapCorrs.append(corr)
    bootstrapCorrs = np.asarray(bootstrapCorrs)
    bx, by = np.histogram(bootstrapCorrs[~np.isnan(bootstrapCorrs)], bins=np.linspace(-1,1,25))
    boot_x.append(bx)
    boot_y.append(by[:-1])
    bootskews.append(scipy.stats.skew(bootstrapCorrs,nan_policy='omit'))

bootskews_cleaned = []
for b in bootskews:
    if type(b) == float:
        bootskews_cleaned.append(b)
    elif type(b) == np.ma.core.MaskedArray:
        bootskews_cleaned.append(b.item)
# %%
