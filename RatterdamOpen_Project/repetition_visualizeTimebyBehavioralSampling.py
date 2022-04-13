# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 18:15:04 2022

@author: whockei1
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
import copy 
from scipy.stats import zscore


datapath  = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

#%% Load recording data with Unit() data
rat, day = 'R781', 'D3'
turns, unit = RepCore.loadTurns(rat, day)
rdf = df[(df.Rat==rat)&(df.Day==day)]

#%% add codes to dataframe of neural data from rat/day
rdf_codes = []
for r, row in rdf.iterrows():
   code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
   rdf_codes.append(code)
rdf = rdf.assign(code=rdf_codes)


#%%


# Remove turnarounds/pivots
ballisticTurnIdx = []
for i in range(1,turns.shape[0]-1):
   row = turns.iloc[i]
   inter = row['Inter']
   # edit 10/2 removing check that last turn's inter wasnt the same,
   # i.e if alley- had a turnaround. since we are looking at things
   # in terms of alley+, only remove a turn if thats where a turnaround was
   if row['Ego'] not in ['0','3'] and turns.iloc[i+1].Inter != inter:
       ballisticTurnIdx.append(i)

refturns = copy.deepcopy(turns) # keep a copy without filtering.
turns = turns.iloc[np.asarray(ballisticTurnIdx)]
    

previousDirection, currentDirection, nextDirection = [], [], [] 
currentAlleys = []

# Create categorical variable of codes
for tnum, turn in refturns.iterrows():

    if tnum < refturns.shape[0]-1:
    
        previousDirection.append(Def.allocodedict[turn['Allo-']])
        currentDirection.append(Def.allocodedict[turn['Allo+']])
        nextDirection.append(Def.allocodedict[refturns.iloc[tnum+1]['Allo+']])
        currentAlleys.append(turn['Alley+'])


codes = [f"{a}{b}{c}" for a,b,c in zip(previousDirection, currentDirection, nextDirection)]
bdf = pd.DataFrame(data={'CurrentAlley':currentAlleys,'Codes':codes})

#%% 

bias_intervals = np.linspace(0,100,21)

variability = []
zscores = []
fanos = []
means= []
rates=[]

for alleyName, adf in bdf.groupby("CurrentAlley"):
        
    abiaslookup = {}
        
    abehav = adf.value_counts()
    biases = abehav.values/abehav.values.sum()/(1/len(abehav.values))
    acodes = [i[1] for i in abehav.index]
    
    for c,b in zip(acodes,biases):
        abiaslookup[c] = b
        
    
    adata = rdf[rdf.Alleys==int(alleyName)]

    alleyzsc = pd.Series(data=zscore(adata.Rate),index=adata.index)

    for acode, acodedata in adata.groupby("code"):
        
        bias = abiaslookup[acode]
        zsc = alleyzsc.loc[acodedata.index]
        fano = np.var(acodedata.Rate)/np.nanmean(acodedata.Rate)
        mean = np.nanmean(acodedata.Rate)
        
        zscores.extend(zsc)
        fanos.append(fano)
        means.append(mean)
        rates.extend(acodedata.Rate)
        variability.extend([bias]*len(zsc))
        #variability.append(bias)
        
 
readout = rates

goodidx = ~np.isnan(readout)
readout = np.asarray(readout)
variability = np.asarray(variability)
readout = readout[goodidx]
variability = variability[goodidx]    
    

olsfit = linregress(variability,readout)
xfit = np.linspace(0, 5,100)
yfit = [(olsfit.slope*i)+olsfit.intercept for i in xfit]
plt.plot(xfit,yfit)
plt.scatter(variability, readout,c='k')
    
    
    
