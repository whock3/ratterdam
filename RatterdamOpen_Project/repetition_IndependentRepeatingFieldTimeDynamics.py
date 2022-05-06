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



alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

fieldTimeCorrs = []
x,y = [], []

i = 0

for rat, rdf in alleydf.groupby("Rat"):
    for day, ddf in rdf.groupby("Day"):
        for orien, oriendf in ddf.groupby("Orientation"):
            for cid,cell in oriendf.groupby("CellID"):
                if np.unique(cell.Repeating)[0] == True:
                    
                    field_timeseries = []


                    start = cell.StartTimes.min()
                    end = cell.StartTimes.max()
                    
                    interp_ts = np.linspace(start,end, 100)
                                                            
                    for fid, field in cell.groupby("FieldID"):
                        
                        startTimes = field.StartTimes
                        rates = field.Rate
                        
                        pcf = pchip(startTimes, rates)
                        pcf_interp = pcf(interp_ts)
                        
                        pcf_interp[pcf_interp>50] = 0 # this is a very rough way of doing it
                        
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


                        

