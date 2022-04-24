# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:41:53 2022

@author: whockei1

Script to create a set of plots for each oriented field.
These plots relate to the hypothesized relationship between behavioral sampling and directional tuning
Three plots for each cell: 1) time versus directionality 
                            2) time versus #unique paths in that window (PCN paths)
                            3) time versus total # paths in that window (PCN paths)

Using filtered visit csv as data, not raw data


"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import itertools
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
from scipy.stats import sem
import repetition_manuscript_defaults as MDef
from importlib import reload
import numpy.ma as ma 

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

# rat_list = ['R781']
# day_list = ['D3']


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)


plotSaveTimeseries = False 

ts = util.genTimestamp()

if 'Code' not in alleydf.columns:
    codes = []
    for r, row in alleydf.iterrows():
       code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
       codes.append(code)
    alleydf = alleydf.assign(Code=codes)
    
directionality_uniquePaths_correlation = []
directionality_totalPaths_correlation = []
directionality_bias_correlation = []

avgRate_uniquePaths_correlation = []
avgRate_totalPaths_correlation = []
avgRate_bias_correlation = []


for rat,day in zip(rat_list, day_list):
    
    rdf = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]
    rdf.StartTimes = (rdf.StartTimes - rdf.StartTimes.min())/1e6 # s, ref'd to 1st sample
    
    window = 15*60
    offset = 1*60
    wins = []
    begin = 0
    stop = False
    
    while not stop:
        a,b = begin, begin + window
        if b < np.ceil(rdf.StartTimes.max()):
            wins.append((a,b))
            begin += offset
        else:
            stop = True

    for orienName, orienDf in rdf.groupby("Orientation"):
        for fid, fieldDf in orienDf.groupby("FieldID"):
            
            time = []
            directionality_time = []
            uniquePaths_time = []
            totalPaths_time = []
            bias_time = []
            avgRate_time = []
            
            for win in wins:
                windf = fieldDf[(fieldDf.StartTimes>win[0])&(fieldDf.StartTimes<=win[1])]
                 
                time.append(win[0])
                dirs = np.unique(windf.CurrDir)
                if dirs.shape[0] > 1:
                    diff = abs(windf[windf.CurrDir==dirs[0]].Rate.mean()-windf[windf.CurrDir==dirs[1]].Rate.mean())    
                    pathcount = 0
                    for cname, code in windf.groupby("Code"):
                        pathcount +=1
                    directionality_time.append(diff)
                    uniquePaths_time.append(pathcount) 
                    totalPaths_time.append(windf.shape[0])
                    ndira = windf[windf.CurrDir==dirs[0]].shape[0]
                    ndirb = windf[windf.CurrDir==dirs[1]].shape[0]
                    bias_time.append(max(ndira,ndirb)/windf.shape[0])
                    avgRate_time.append(windf.Rate.mean())
                else:
                    directionality_time.append(np.nan)
                    uniquePaths_time.append(np.nan)
                    totalPaths_time.append(np.nan)
                    bias_time.append(np.nan)
                    avgRate_time.append(np.nan)
                    

            directionality_uniquePaths_correlation.append(ma.corrcoef(ma.masked_invalid(directionality_time),
                                                                      ma.masked_invalid(uniquePaths_time)).data[0,1])
            
            directionality_totalPaths_correlation.append(ma.corrcoef(ma.masked_invalid(directionality_time),
                                                                      ma.masked_invalid(totalPaths_time)).data[0,1])
            
            directionality_bias_correlation.append(ma.corrcoef(ma.masked_invalid(directionality_time),
                                                                      ma.masked_invalid(bias_time)).data[0,1])
            
            avgRate_uniquePaths_correlation.append(ma.corrcoef(ma.masked_invalid(avgRate_time),
                                                                      ma.masked_invalid(uniquePaths_time)).data[0,1])
            
            avgRate_totalPaths_correlation.append(ma.corrcoef(ma.masked_invalid(avgRate_time),
                                                                      ma.masked_invalid(totalPaths_time)).data[0,1])
            
            avgRate_bias_correlation.append(ma.corrcoef(ma.masked_invalid(avgRate_time),
                                                                      ma.masked_invalid(bias_time)).data[0,1])
                    
                    
            if plotSaveTimeseries:
                    
            
                fig, axes = plt.subplots(5,1, figsize=(10,15))
                for ax, data, label, color in zip(fig.axes, 
                                                  [avgRate_time,
                                                   directionality_time, 
                                                   uniquePaths_time, 
                                                   totalPaths_time, 
                                                   bias_time
                                                   ],
                                                  ["Average Rate",
                                                   "Directional Firing", 
                                                   "Unique Paths", 
                                                   "Total Paths", 
                                                   "Behavioral\n Direction\n Bias"],
                                                  ['black',
                                                   'grey',
                                                   'darkblue',
                                                   'cornflowerblue',
                                                   'purple']
                                                  ):
                    
                    ax.plot(time, data, color=color,marker='^',linestyle='-',linewidth=2,markersize=10)
                    ax.set_ylabel(label, fontsize=20)
                    ax.set_xlabel("Time (s)", fontsize=20)
                    
                
                cellname = np.unique(fieldDf.CellName)[0]
                clustName = cellname.replace("\\","_")
                field = np.unique(fieldDf.FieldNum)[0]
                plt.suptitle(f"{rat}{day} {clustName} Field {field}, {orienName} Orientation", fontsize=25)
                fname = f"{ts}_{rat}{day}_{clustName}_{field}{orienName}_{window/60}min_{offset/60}step_samplingDirTuningTime"
                plt.savefig(f"E:\\Ratterdam\\temp\\PathTuningByLocation\\perField_samplingAndTuningOverTime\\{fname}.png",dpi=300)
                plt.close()
                
            
            
#%% Plotting scatter of correlations between pairs of timeseries above for all oriented fields 

fig, ax = plt.subplots(2,3)
data_vectors = [directionality_uniquePaths_correlation,
        directionality_totalPaths_correlation,
        directionality_bias_correlation,
        avgRate_uniquePaths_correlation,
        avgRate_totalPaths_correlation,
        avgRate_bias_correlation
        ]

labels = ['Directionality and Unique Paths',
          'Directionality and Total Paths',
          'Directionality and Sampling Bias',
          'Avg Rate and Unique Paths',
          'Avg Rate and Total Paths',
          'Avg Rate and Sampling Bias'
            ]

for cax, data, label in zip(fig.axes, data_vectors, labels):
    
    cax.hist(data)
    cax.set_title(label)
    
    
#%% Shuffling test of a relationship between firing (avg or directional) and behavioral sampling

nshuffles = 1000

ts = util.genTimestamp()

if 'Code' not in alleydf.columns:
    codes = []
    for r, row in alleydf.iterrows():
       code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
       codes.append(code)
    alleydf = alleydf.assign(Code=codes)

window = 10*60
offset = 2*60
wins = []
begin = 0
stop = False

data_id = 0 # keys for our dict of the time series
timeseries_dict = {}

sigUnits = []
totalModels = 0

for rat,day in zip(rat_list, day_list):

    rdf = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]
    rdf.StartTimes = (rdf.StartTimes - rdf.StartTimes.min())/1e6 # s, ref'd to 1st sample
    
    while not stop:
        a,b = begin, begin + window
        if b < np.ceil(rdf.StartTimes.max()):
            wins.append((a,b))
            begin += offset
        else:
            stop = True
    
    for orienName, orienDf in rdf.groupby("Orientation"):
        for fid, fieldDf in orienDf.groupby("FieldID"):
            
            time = []
            directionality_time = []
            uniquePaths_time = []
            totalPaths_time = []
            bias_time = []
            avgRate_time = []
            
            for win in wins:
                windf = fieldDf[(fieldDf.StartTimes>win[0])&(fieldDf.StartTimes<=win[1])]
                
                time.append(win[0])
                dirs = np.unique(windf.CurrDir)
                if dirs.shape[0] > 1:
                    diff = abs(windf[windf.CurrDir==dirs[0]].Rate.mean()-windf[windf.CurrDir==dirs[1]].Rate.mean())    
                    pathcount = 0
                    for cname, code in windf.groupby("Code"):
                        pathcount +=1
                    directionality_time.append(diff)
                    uniquePaths_time.append(pathcount) 
                    totalPaths_time.append(windf.shape[0])
                    ndira = windf[windf.CurrDir==dirs[0]].shape[0]
                    ndirb = windf[windf.CurrDir==dirs[1]].shape[0]
                    bias_time.append(max(ndira,ndirb)/windf.shape[0])
                    avgRate_time.append(windf.Rate.mean())
                else:
                    # for shuffle, we don't want nans. Want to shuffle firing rates
                    #without worrying the positions where nans are arent same in other, unshuffled vectors
                    
                    pass
                    # directionality_time.append(np.nan)
                    # uniquePaths_time.append(np.nan)
                    # totalPaths_time.append(np.nan)
                    # bias_time.append(np.nan)
                    # avgRate_time.append(np.nan)
                    
            timeseries_dict[data_id] = {'directionalFR':directionality_time,
                                        'avgRate':avgRate_time,
                                        'uniquePaths':uniquePaths_time,
                                        'totalPaths':totalPaths_time,
                                        'biases':bias_time}
            
            xdf = pd.DataFrame(data={'uniquePaths':uniquePaths_time,
                                        'totalPaths':totalPaths_time,
                                        'biases':bias_time
                                    })
            try:
                m = sm.OLS(avgRate_time, xdf)
                mfit = m.fit()
                if mfit.pvalues['uniquePaths'] < 0.05/4:
                    sigUnits.append(['uniquePaths', fid])
                elif mfit.pvalues['totalPaths'] < 0.05/4:
                    sigUnits.append(['totalPaths', fid])
                    
                totalModels += 1 
            except:
                pass
            
            data_id += 1 
 
        
shuffles_directionality_uniquePaths_correlation = []
shuffles_directionality_totalPaths_correlation = []
shuffles_directionality_bias_correlation = []

shuffles_avgRate_uniquePaths_correlation = []
shuffles_avgRate_totalPaths_correlation = []
shuffles_avgRate_bias_correlation = []

shuffled_uniquePaths_beta = []
shuffled_totalPaths_beta = [] 

shuffSigUnits = []
totalShuffModels = 0 

for s in range(nshuffles):
    
    print(s)
    
    ss_directionality_uniquePaths_correlation = []
    ss_directionality_totalPaths_correlation = []
    ss_directionality_bias_correlation = []
    
    ss_avgRate_uniquePaths_correlation = []
    ss_avgRate_totalPaths_correlation = []
    ss_avgRate_bias_correlation = []
    
    ss_sigUnits = []
    
    for dID,field_data in timeseries_dict.items():
        dirFr = np.random.permutation(field_data['directionalFR'])
        avgFr = np.random.permutation(field_data['avgRate'])
        uniquePaths = field_data['uniquePaths']
        totalPaths = field_data['totalPaths']
        bias = field_data['biases']
        
        xdf = pd.DataFrame(data={'uniquePaths':uniquePaths,
                                        'totalPaths':totalPaths,
                                        'biases':bias
                                    })
        
        try:
            m = sm.OLS(avgFr, xdf)
            mfit = m.fit()
            if mfit.pvalues['uniquePaths'] < 0.05/4:
                ss_sigUnits.append(['uniquePaths', fid])
            elif mfit.pvalues['totalPaths'] < 0.05/4:
                ss_sigUnits.append(['totalPaths', fid])
            totalShuffModels += 1 
                
        except:
            pass

                        
        ss_directionality_uniquePaths_correlation.append(ma.corrcoef(ma.masked_invalid(dirFr),
                                                                  ma.masked_invalid(uniquePaths)).data[0,1])
        
        ss_directionality_totalPaths_correlation.append(ma.corrcoef(ma.masked_invalid(dirFr),
                                                                  ma.masked_invalid(totalPaths)).data[0,1])
        
        ss_directionality_bias_correlation.append(ma.corrcoef(ma.masked_invalid(dirFr),
                                                                  ma.masked_invalid(bias)).data[0,1])
        
        ss_avgRate_uniquePaths_correlation.append(ma.corrcoef(ma.masked_invalid(avgFr),
                                                                  ma.masked_invalid(uniquePaths)).data[0,1])
        
        ss_avgRate_totalPaths_correlation.append(ma.corrcoef(ma.masked_invalid(avgFr),
                                                                  ma.masked_invalid(totalPaths)).data[0,1])
        
        ss_avgRate_bias_correlation.append(ma.corrcoef(ma.masked_invalid(avgFr),
                                                                  ma.masked_invalid(bias)).data[0,1])
        
                
          
    shuffles_directionality_uniquePaths_correlation.append(np.nanpercentile(np.abs(ss_directionality_uniquePaths_correlation),95))
    shuffles_directionality_totalPaths_correlation.append(np.nanpercentile(np.abs(ss_directionality_totalPaths_correlation),95))
    shuffles_directionality_bias_correlation.append(np.nanpercentile(np.abs(ss_directionality_bias_correlation),95))
    
    shuffles_avgRate_uniquePaths_correlation.append(np.nanpercentile(np.abs(ss_avgRate_uniquePaths_correlation),95))
    shuffles_avgRate_totalPaths_correlation.append(np.nanpercentile(np.abs(ss_avgRate_totalPaths_correlation),95))
    shuffles_avgRate_bias_correlation.append(np.nanpercentile(np.abs(ss_avgRate_bias_correlation),95))
    
    shuffSigUnits.append(len(ss_sigUnits))
                        

#%% Plotting shuffling results 

fig, ax = plt.subplots(1,3, figsize=(12,9))


fig.axes[0].scatter(np.abs(directionality_uniquePaths_correlation), np.abs(directionality_totalPaths_correlation),c='r')
fig.axes[0].scatter(shuffles_directionality_uniquePaths_correlation, shuffles_directionality_totalPaths_correlation,c='k')
uniquepaths_95 = np.nanpercentile(shuffles_directionality_uniquePaths_correlation, 97.5)
totalpaths_95 = np.nanpercentile(shuffles_directionality_totalPaths_correlation, 97.5)
fig.axes[0].vlines(uniquepaths_95, 0, 1, label = '97.5th percentile', color='k', linestyle='--')
fig.axes[0].hlines(totalpaths_95, 0, 1, color='k', linestyle='--')
fig.axes[0].set_xlabel("Abs Val of Correlation between directional rate and unique paths", fontsize=22)
fig.axes[0].set_ylabel("Abs Val of Correlation between directional firing rate and total paths", fontsize=22)
fig.axes[0].set_title("Directional Firing Rate", fontsize=30)
fig.axes[0].tick_params(axis='both', which='major', labelsize=22)

plt.legend()

fig.axes[1].scatter(np.abs(avgRate_uniquePaths_correlation), np.abs(avgRate_totalPaths_correlation),c='r')
fig.axes[1].scatter(shuffles_avgRate_uniquePaths_correlation, shuffles_avgRate_totalPaths_correlation,c='k')
uniquepaths_95 = np.nanpercentile(shuffles_avgRate_uniquePaths_correlation, 97.5)
totalpaths_95 = np.nanpercentile(shuffles_avgRate_totalPaths_correlation, 97.5)
fig.axes[1].vlines(uniquepaths_95, 0, 1, label = '97.5th percentile', color='k', linestyle='--')
fig.axes[1].hlines(totalpaths_95, 0, 1, color='k', linestyle='--')
fig.axes[1].set_xlabel("Abs Val of Correlation between average rate and unique paths", fontsize=22)
fig.axes[1].set_ylabel("Abs Val of Correlation between average firing rate and total paths", fontsize=22)
fig.axes[1].set_title('Average Firing Rate', fontsize=30)
fig.axes[1].tick_params(axis='both', which='major', labelsize=22)

plt.legend()


fig.axes[1].scatter(np.abs(avgRate_uniquePaths_correlation), np.abs(avgRate_totalPaths_correlation),c='r')
fig.axes[1].scatter(shuffles_avgRate_uniquePaths_correlation, shuffles_avgRate_totalPaths_correlation,c='k')
uniquepaths_95 = np.nanpercentile(shuffles_avgRate_uniquePaths_correlation, 97.5)
totalpaths_95 = np.nanpercentile(shuffles_avgRate_totalPaths_correlation, 97.5)
fig.axes[1].vlines(uniquepaths_95, 0, 1, label = '97.5th percentile', color='k', linestyle='--')
fig.axes[1].hlines(totalpaths_95, 0, 1, color='k', linestyle='--')
fig.axes[1].set_xlabel("Abs Val of Correlation between average rate and unique paths", fontsize=22)
fig.axes[1].set_ylabel("Abs Val of Correlation between average firing rate and total paths", fontsize=22)
fig.axes[1].set_title('Average Firing Rate', fontsize=30)
fig.axes[1].tick_params(axis='both', which='major', labelsize=22)

plt.legend()


