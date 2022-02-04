# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:45:07 2022

@author: whockei1
"""
import pandas as pd, matplotlib.pyplot as plt, numpy as np
import ratterdam_RepetitionCoreFx as RepCore

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


