# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:45:07 2022

@author: whockei1
"""
import pandas as pd, matplotlib.pyplot as plt, numpy as np
import ratterdam_RepetitionCoreFx as RepCore
import newAlleyBounds as nab 
import pickle 
import copy 

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-135920_superPopAlleyBehaviorResponse_1.5vfilt.csv"
alleydf = pd.read_csv(alleydatapath)

interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
interdf = pd.read_csv(interdatapath)

r,d = 'R859','D2'
alleydf = alleydf[(alleydf.Rat==r)&(alleydf.Day==d)]
interdf = interdf[(interdf.Rat==r)&(interdf.Day==d)]

#%% behavior (still based on where there is at least one field)

ntrajs = {}
ntj = []

for rname, rat in alleydf.groupby("Rat"):
    for daynum, day in rat.groupby("Day"):

        for alleyNum, alley in day.groupby("Alleys"):
            
            alley_trajs = []
            
            for codename, code in alley.groupby("Code"):
                            
                alley_trajs.append((codename,code.shape[0]))
                
            ntrajs[alleyNum] = [len(alley_trajs)] + alley_trajs
            ntj.append(len(alley_trajs))
            
            
#%% 

ntrajs = []
diffs = []

traj_sample_thresh = 0
included_field_chunks = [] # list of tuples: [(fid,alley)]

for rname, rat in alleydf.groupby("Rat"):
    for daynum, day in rat.groupby("Day"):

        for alleyNum, alley in day.groupby("Alleys"):
            
            alley_trajs = []
            
            for codename, code in alley.groupby("Code"):
                alley_trajs.append(code.shape[0])
                
            #alley_trajs = [i for i in alley_trajs if i > traj_sample_thresh]
            
            alleyNumTrajs = len(alley_trajs)
            alleySpreadTrajs = (max(alley_trajs)-min(alley_trajs))/len(alley_trajs)
            
            for fid, field in alley.groupby("FieldID"):
                
                ntrajs.append(alleyNumTrajs)

included_field_chunks = np.asarray(included_field_chunks)   




#%% 
turns,unit = RepCore.loadTurns('R859','D2')
turns['Code'] = turns['Allo-']+turns['Ego']+turns['Allo+']

ntj = []
for alleyNum, alley in turns.groupby("Alley-"):
    alley_trajs = []
    for codename, code in alley.groupby("Code"):
        alley_trajs.append((codename,code.shape[0]))
    ntj.append(len(alley_trajs))
    
plt.figure();plt.hist(ntj)
        

#%% 

import pandas as pd, matplotlib.pyplot as plt, numpy as np

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-135920_superPopAlleyBehaviorResponse_1.5vfilt.csv"
alleydf = pd.read_csv(alleydatapath)

rat,day = 'R886','D2'
rdf = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]

rdf.StartTimes = (rdf.StartTimes - rdf.StartTimes.min())/1e6 # s, ref'd to 1st sample

window = 10*60
offset = 2*60
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
        
        
numpaths_time = []
directionality_time = []
bias_time = []

for win in wins:
    
    numpaths_win = []
    directionality_win = []
    bias_win = []
    
    windf = rdf[(rdf.StartTimes>win[0])&(rdf.StartTimes<=win[1])]    

    # for aNum, alley in windf.groupby("Alleys"):
    #     pathcount = 0
    #     for c, code in alley.groupby("Code"):
    #         pathcount += 1
    #     numpaths_win.append(pathcount)
        
    for fid, field in windf.groupby("FieldID"):
        
        for o, ofield in field.groupby("Orientation"):

            dirs = np.unique(ofield.CurrDir)
            if len(dirs)>1:
                diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                directionality_win.append(diff)
                pathcount = 0
                for cname, code in ofield.groupby("Code"):
                    pathcount +=1
                numpaths_win.append(pathcount)
                ndira = ofield[ofield.CurrDir==dirs[0]].shape[0]
                ndirb = ofield[ofield.CurrDir==dirs[1]].shape[0]
                bias_win.append(max(ndira,ndirb)/ofield.shape[0])
            
                
            
    numpaths_time.append(numpaths_win)
    directionality_time.append(directionality_win)
    bias_time.append(bias_win)
            
    
#%% plot   

#plt.figure()
# plt.hist([i for j in numpaths_time for i in j])
# plt.title("Number of paths")

    
# plt.figure()
# for i,j in zip(numpaths_time, directionality_time):
#     plt.scatter(i,j,color='k',alpha=0.3)
# plt.title("Number of paths vs directionality")
plt.figure()
for i,j in zip(bias_time, directionality_time):
    plt.scatter(i,j,color='k',alpha=0.3)
plt.title("Directional bias vs directionality")


alldir  = np.asarray([i for j in directionality_time for i in j])
allpath = np.asarray([i for j in numpaths_time for i in j])

pathdict = {}
for n in np.unique(allpath):
    pathdict[n] = [alldir[i] for i in range(len(alldir)) if allpath[i]==n]
    
fig, ax = plt.subplots(1,2,figsize=(24,10))
for k,v in pathdict.items():
    ax[0].scatter([k]*len(v)+np.random.normal(0,0.1,size=len(v)),v,color='k')
    vp = ax[0].violinplot([v],positions=[k])
    vp['bodies'][0].set_facecolor('grey')
    vp['bodies'][0].set_edgecolor('black')
    vp['bodies'][0].set_linewidth(2)
    vp['bodies'][0].set_alpha(0.7)
    for el in ['cbars','cmaxes','cmins']:
        vp[el].set_color('black')
    
ax[0].set_xlabel("Number of Unique Paths Performed in Time Window", fontsize=24)
ax[0].set_ylabel("Directionality (Abs Mean Difference, Hz) in Alleys Only", fontsize=24)
ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)

ax[1].hist([i for j in numpaths_time for i in j],color='grey',edgecolor='black',linewidth=2)
ax[1].set_xlabel("Number of Unique Paths Performed in Time Window", fontsize=24)
ax[1].set_ylabel("Frequency", fontsize=24)
ax[1].tick_params(axis='both', which='major', labelsize=20)
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)

plt.suptitle(f"{rat}{day} Relationship between Number of Unique Paths and Directionality \n \
          Windowed in Time {window/60}mins window, {offset/60}min offset",
          fontsize=32)



#%% 2-11-22 Looking into why there is a bimodal distribution in fields' peak firing rate (panel from fig 1)
import numpy as np, matplotlib.pyplot as plt, pandas as pd
import utility_fx as util 
# don't filter by passes, we're not making any comparisons or looking at directionality
alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-135920_superPopAlleyBehaviorResponse_1.5vfilt.csv"
alleydf = pd.read_csv(alleydatapath)

df = alleydf
#df = alleydf[alleydf.Repeating==False]
peakrates = []
for fid, field in df.groupby("FieldID"):
    peakrates.append(np.nanmax(field.Rate))
    
nHzPerBin = 1 # histogram bins should contain how many hertz, e.g. 2Hz increment per bin
histMax = np.max(peakrates)
b=np.linspace(0,histMax,int(np.ceil(histMax/nHzPerBin)))# hist bins

fig, ax = plt.subplots(figsize=(20,15))
ax.hist(peakrates,facecolor='blue',bins=b,alpha=0.7)



#%% 2-11-22 Speed and behavior analysis that will be folded into figure 1 
import ratterdam_RepetitionCoreFx as RepCore
import utility_fx as util 
rat, day = 'R781', 'D3'
datapath = f"E:\\Ratterdam\{rat}\\{rat}_RatterdamOpen_{day}\\"
clustList,_ = util.getClustList(datapath)
unit = RepCore.loadRepeatingUnit(rat, day, clustList[0]) # grab any unit 
pos = unit.position



#%% 2022-03-05 Trying to manually score rewards from R781 Repetition experiment
# Day 3 and D4 

#D3 
r781d3_candidateRewards = np.asarray([
                        3054990323895,
                        3055107145278,
                        3055472821503,
                        3055510024989,
                        3055569220819,
                        3055612631790,
                        3055875602556,
                        3055960489096,
                        3056206643933,
                        3056259799571,
                        3056366009509,
                        3056712933916,
                        3056824747382,
                        3056897222107,
                        3057045074765,
                        3057058589349,
                        3057216153169,
                        3057264801280,
                        3057728178538,
                        3057765684948,
                        3058045774350
                    ])


#D4 
r781d4_candidateRewards = np.asarray([3145509831365,
                    3145538126675,
                    3145583972805,
                    3145614504896,
                    3145719012407,
                    3145791686207,
                    3145884483369,
                    3145930231156,
                    3145975075603,
                    3146024526025,
                    3146072476521,
                    3146082186239,
                    3146154394280,
                    3146470952705,
                    3146485466970,
                    3146527312532,
                    3146565250144,
                    3146605291572,
                    3146748037783,
                    3146783075166,
                    3146896592500,
                    3146917113667,
                    3147141143652,
                    3147371147666,
                    3147477557649,
                    3147597877121,
                    3147650768884,
                    3147863286550,
                    3147993587545,
                    3148040336885,
                    3148059890112,
                    3148111409591,
                    3148311115888,
                    3148360233182,
                    3148805091252              
                    
                    ])


# plt.vlines(candidateRewards,0,600,'r', label='Manually-detected rewards')
# plt.vlines(rewards,0,600,'g', label = 'TTL rewards')
# plt.scatter(r781d4unit.position[:,0],r781d4unit.position[:,2],c='k',label='behavior (single axis')
# plt.legend()
# plt.xlabel("Timestamp (us)",fontsize=25)
# plt.title("Detectd Rewards for R781 D4",fontsize=32)



#%% Making new superpop code 

superpopulation = {}

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

for rat, day in zip(rat_list, day_list):
    
    if rat not in superpopulation.keys():
        superpopulation[rat] = {}
    superpopulation[rat][day] = {}
    
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    rewards = RepCore.readinRewards(rat, day)
    ratborders = nab.loadAlleyBounds(rat, day)
    
    # Code 0 means an error and a discontinuity in behavioral tracking
    # this happens infrequently enough that it's not a huge problem
    # causes are typically loss of camera tracking briefly. 
    turns = turns[turns.Ego!='0']

    codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
    
    # Remove turnarounds/pivots
    ballisticTurnIdx = []
    for i in range(1,turns.shape[0]-1):
       row = turns.iloc[i]
       inter = row['Inter']
       # edit 10/2 removing check that last turn's inter wasnt the same,
       # i.e if alley- had a turnaround. since we are looking at things
       # in terms of alley+, only remove a turn if thats where a turnaround was
       if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter:
           ballisticTurnIdx.append(i)
    
    refturns = copy.deepcopy(turns) # keep a copy without filtering.
    turns = turns.iloc[np.asarray(ballisticTurnIdx)]
    
    superpopulation[rat][day]['units'] = population
    superpopulation[rat][day]['turns'] = turns
    superpopulation[rat][day]['refturns'] = refturns
    
tstamp = util.genTimestamp()

with open(f"E:\\Ratterdam\\R_data_repetition\\{tstamp}_superPopulationRepetition.pickle", "wb") as f:
    pickle.dump(superpopulation, f)