# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:41:42 2021

@author: whockei1

Map switching across the population 
"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
import pickle
import matplotlib.path as path

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)
    

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]


#%% plot all fields separately, look for common dynamics

rat, day = 'R781', 'D3'
population, refturns, turns = superpop[rat][day]['units'], superpop[rat][day]['refturns'], superpop[rat][day]['turns']


fields = []
names = []
for unit in population.values():
    for field in unit.fields:
        fields.append(field)
        names.append(f"{unit.name} Repeating = {unit.repeating}")
        
ncols = 6
fmin = min([min(f[:,0]) for f in fields])
fmax = max([max(f[:,0]) for f in fields])
fig, _ax = plt.subplots(int(np.ceil(len(fields)/ncols)),ncols, figsize=(12,15))
for i, (name,field) in enumerate(zip(names,fields)):
    ax = fig.axes[i]
    ax.plot(field[:,0], field[:,1])
    ax.set_xlim([fmin,fmax])
    ax.set_title(name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
fig.subplots_adjust(wspace=0.1,hspace=1)

#%%
#######
####### Redoing looking at field FR compared to accumulated dir bias 
#######
#%% Redoing looking at fields compared to accumalated dir bias
# break field down by orientation, look at moiety on each alley

# Define alley biases
import newAlleyBounds as nab
rat, day = 'R859', 'D1'
population, refturns, turns = superpop[rat][day]['units'], superpop[rat][day]['refturns'], superpop[rat][day]['turns']

datapath =  "E:\\Ratterdam\\R_data_repetition\\211005_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)

daydf = df[(df.Rat==rat)&(df.Day==day)]

# get rid of turnarounds
daydf = daydf[daydf.Traversal==True]

datafile = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
ratborders = nab.loadAlleyBounds(rat, day)
with open(datafile+"sessionEpochInfo.txt","r") as f:
    data = f.readlines()
session_start, session_end = [float(i) for i in data[0].split(',')]

wnSizeMin = 10
wnStepMin = 1
wnSize= wnSizeMin*60*1e6 # set num mins you want in a window
wnStep = wnStepMin*60*1e6
passThresh = 1
time_windows = []


stop = False
begin = session_start
while not stop:
    a,b = begin, begin + wnSize
    if b <= session_end:
        time_windows.append((a,b))
        begin += wnStep
    else:
        stop = True
time_windows = np.asarray(time_windows)

dataset = superpop[rat][day]
turns = dataset['turns']

# alleyBiases = {str(i):[] for i in range(17)}
# for alley in alleyBiases.keys():
    
#     dirs = np.unique(turns[turns['Alley+']==alley]['Allo+'])
#     dirs = dirs[dirs!='0']
#     if len(dirs)>1:
#         for win in time_windows:
#             turnwin = turns[(turns["Alley+"]==alley)
#                             &(turns["Ts entry"].astype(float)>win[0])
#                             &(turns["Ts exit"].astype(float)<=win[1])]
            
#             dirA = turnwin[turnwin['Allo+'] == dirs[0]].shape[0]
#             dirB = turnwin[turnwin['Allo+'] == dirs[1]].shape[0]
            
#             try:
#                 alleyBiases[alley].append([(win[1]+win[0])/2,max(dirA,dirB)/(dirA+dirB)])
#             except:
#                 alleyBiases[alley].append([(win[1]+win[0])/2,np.nan])

# for key in alleyBiases.keys():
#     alleyBiases[key] = np.asarray(alleyBiases[key])    
        
    
#%calc dir bias accumulation 
field_biasaccum = {str(i):[[float(turns.iloc[0]['Ts exit']), 0]] for i in range(17)}

for t,turn in turns.iterrows():
    if t < turns.iloc[-1].name:
        if turn['Allo+'] == '2' or turn['Allo+'] == '1':
            ts = (float(turn['Ts entry']) + float(refturns.iloc[t+1]['Ts exit']))/2
            field_biasaccum[turn['Alley+']].append([ts, field_biasaccum[turn['Alley+']][-1][1]+1])
        elif turn['Allo+'] == '3' or turn['Allo+'] == '4':
            ts = (float(turn['Ts entry']) + float(refturns.iloc[t+1]['Ts exit']))/2
            field_biasaccum[turn['Alley+']].append([ts, field_biasaccum[turn['Alley+']][-1][1]-1])

for key in field_biasaccum.keys():
    field_biasaccum[key] = np.asarray(field_biasaccum[key])                                          


diffs = {i:[] for i in range(17)}

# do with filtered visits

for alley in range(17):
    
    for orien in ["V","H"]: 
        
        oriendf = daydf[(daydf.Orientation==orien)&(daydf.Alleys==alley)]
        
        for fname, field in oriendf.groupby('FieldID'):
            alleydiffs = []
    
            for widx, win in enumerate(time_windows):
                d = np.unique(field.CurrDir)
                if d.shape[0]>1:
                    wf = field[(field.StartTimes>win[0])&(field.StartTimes<=win[1])]
                    dirA, dirB = wf[wf.CurrDir==d[0]], wf[wf.CurrDir==d[1]]
                    if dirA.shape[0] >= passThresh and dirB.shape[0] >= passThresh:
                        try:
                            diff = abs(np.mean(dirA.Rate)-np.mean(dirB.Rate))/np.mean(field.Rate)
                            alleydiffs.append([(win[1]+win[0])/2,diff])
                        except Exception as e:
                            print(e)
                            
            diffs[alley].append([fname, np.asarray(alleydiffs)])



fig, ax = plt.subplots(5,4)
for i,alley in enumerate(field_biasaccum.keys()):
    ax = fig.axes[i]
    ax.set_title(f"Alley {alley}")
    ax.plot(field_biasaccum[alley][:,0], field_biasaccum[alley][:,1],color='black',linestyle='--',linewidth=2)
    ax2 = ax.twinx()
    for flabel, factivity in diffs[int(alley)]:
        try:
            ax2.plot(factivity[:,0], factivity[:,1],marker='.',linestyle='-',markersize=5)
        except:
            pass
    ax.set_ylabel("Behavioral Direction\n Bias Accumulation")
    ax2.set_ylabel("Field Directional Tuning")
plt.subplots_adjust(wspace=0.5,hspace=0.3)


## Below is code to compute FR over time and plot against
## direction bias accumulation. This doesnt pull out much, perhaps not surprisingly

#%% calc visit FR for each field broken down by alley and field bounds
field_rates = {i:[] for i in alleyBiases.keys()}
for unit in population.values():
    for i,(foverlap,perim) in enumerate(zip(unit.overlaps, unit.perimeters)):
        contour = path.Path(perim)
        field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
        field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
        
        for region in foverlap:
            if type(region)==int:
                rates = []
                for t, turn in turns.iterrows():
                    if t < turns.shape[0]-1:
                        if turn['Alley+'] == str(region):
                            
                            winStart, winEnd = float(turn['Ts entry']), float(refturns.iloc[t+1]['Ts exit'])
                            winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                            winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                            if winPos.shape[0] > 0:
                                rates.append(((winEnd+winStart)/2,winSpikes.shape[0]/((winEnd-winStart)/1e6)))
                if rates != []:              
                    field_rates[str(region)].append([f"{unit.name}_{i}_{region}", np.asarray(rates)])
                                 
# Plot alleys: directional bias, bias accum, and field rates (restricted to infield and inalley)       

fig, ax = plt.subplots(5,4)
for i,alley in enumerate(alleyBiases.keys()):
    ax = fig.axes[i]
    ax.set_title(f"Alley {alley}")
    #ax.plot(alleyBiases[alley][:,0],alleyBiases[alley][:,1],color='black',linewidth=2)
    ax.plot(field_biasaccum[alley][:,0], field_biasaccum[alley][:,1],color='black',linestyle='--',linewidth=2)
    ax2 = ax.twinx()
    for flabel, factivity in field_rates[alley]:
        ax2.plot(factivity[:,0], util.weird_smooth(factivity[:,1],2))
        #ax2.text(factivity[-1,0],factivity[-1,1],flabel)
    ax.set_ylabel("Behavioral Direction Bias Accumulation")
    ax2.set_ylabel("Firing Rate")