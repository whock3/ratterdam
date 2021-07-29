# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:31:52 2021

@author: whockei1

Repetition Project
Hierarchical Clustering of Trajectory Responses
"""

#%% Imports and defaults


import numpy as np
import utility_fx as util
import os
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import newAlleyBounds as nab
import math
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
from collections import Counter
from sklearn.manifold import Isomap
from sklearn.manifold import TSNE
from mpl_toolkits import mplot3d


region_sets = {'RS1':[12],  #decode traj and dir
                'RS2':[3,5],  # decode traj and dir
                'RS3':[14,11], # decode traj and dir
                'RS4':[0,4,6,15,13,10], #decode dir (E-W)
                'RS5':[2,16,7,9]  # decode dir (N-S)
                }

codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}


#%% Load Data

rat, day = 'R859', 'D1'
population, turns = RepCore.loadRecordingSessionData(rat, day)

 #%% Create data arrays
#Here's the logic. If you want the acivity from when the animal was on a given
# alley on a certain pass, you take the ts entry of turn n to the ts exit of
# turn n+1 and that corresponds to time spent on alley+. 
currentDir, previousDir, nextDir = [], [], []
ego = []
currentAlley = []
traj = []
X = np.empty((0, len(population)))

for t in range(1,turns.shape[0]-1):
    
    turn_nm1 = turns.iloc[t-1]
    turn = turns.iloc[t]
    turn_np1 = turns.iloc[t+1]
            
    # cD= turn['Ego'] # currentDir value this turn
    # pD = turn_nm1['Ego'] # previous ''
    # nD = turn_np1['Ego'] # next ''
    cD = turn['Allo+']
    pD = turn['Allo-']
    nD = turn_np1['Allo+']
    
    ego.append(turn['Ego'])
    
    currentAlley.append(turn['Alley+']) # use this later to get visits to regions in a certain set 
    
    start, stop = float(turn['Ts entry']), float(turn_np1['Ts exit'])
    duration = (stop-start)/1e6
    
    popvector = []
    for unitname, unit in population.items():
        #basic stuff, get the firing rate (# spikes / time on alley ) for each unit and append
            spike_count = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]
            rate = spike_count / duration
            popvector.append(rate)

    popvector = np.asarray(popvector)
    X = np.vstack((X, popvector))
    currentDir.append(cD)
    previousDir.append(pD)
    nextDir.append(nD)
    traj.append(f"{pD}{cD}{nD}")
    
currentDir = np.asarray(currentDir)
nextDir = np.asarray(nextDir)
previousDir = np.asarray(previousDir)
traj = np.asarray(traj)
currentAlley = np.asarray(currentAlley)
ego = np.asarray(ego)

#convert numeric codes to cardinal directions 
traj = [f"{codedict[t[0]]}{codedict[t[1]]}{codedict[t[2]]}" for t in traj]
traj = np.asarray(traj)

#%% Dimensionality Reduction 3d (Isomap) and visualize 
rsl = "RS3"
rs = region_sets[rsl]
Xsub,targetsub = np.empty((0, len(population))), []
for r in rs:
    Xsub = np.vstack((Xsub, X[currentAlley==str(r),:]))
    targetsub.extend(traj[currentAlley==str(r)])
    
targetsub = np.asarray(targetsub)
n=3
iso = Isomap(n_components=n)

Xtransf = iso.fit_transform(Xsub)

ax = plt.figure(figsize=(20,20)).add_subplot(projection='3d')
ax.scatter3D(Xtransf[:,0], Xtransf[:,1], Xtransf[:,2])

for i,txt in enumerate(targetsub):
    ax.text(Xtransf[i,0],Xtransf[i,1], Xtransf[i,2],txt)

plt.title(f"{rat}{day} {rsl}, Isomap dim={n}", fontsize=18)


#%% Dimensionality Reduction 2d (Isomap) and visualize

rsl = "RS1"
rs = region_sets[rsl]
Xsub,targetsub = np.empty((0, len(population))), []
for r in rs:
    Xsub = np.vstack((Xsub, X[currentAlley==str(r),:]))
    targetsub.extend(currentDir[currentAlley==str(r)])
    
targetsub = np.asarray(targetsub)
n=2
iso = Isomap(n_components=n)

Xtransf = iso.fit_transform(Xsub)
plt.scatter(Xtransf[:,0], Xtransf[:,1],s=1)

for i,txt in enumerate(targetsub):
    plt.annotate(txt,(Xtransf[i,0], Xtransf[i,1]))

# #color code the letters
# for i,txt in enumerate(targetsub):
#     for letter,offset in zip(txt,[0,0.1,0.2]):
#         c = {'N':'r','S':'b','E':'g','W':'y','X':'k'}[letter]
#         plt.annotate(letter,(Xtransf[i,0]+offset,Xtransf[i,1]),color=c)


    
plt.title(f"{rat}{day} {rsl}, Isomap n={n}", fontsize=18)

plt.xlabel("Dim 1",fontsize=16)

plt.ylabel("Dim 2", fontsize=16)
    
#%% Dimensionality Reduction 2d (isomap) and Visualize 3d with time
rsl = "RS1"
rs = region_sets[rsl]
Xsub,targetsub = np.empty((0, len(population))), []
for r in rs:
    Xsub = np.vstack((Xsub, X[currentAlley==str(r),:]))
    targetsub.extend(traj[currentAlley==str(r)])
    
targetsub = np.asarray(targetsub)

iso = Isomap(n_components=2)

Xtransf = iso.fit_transform(Xsub)

ax = plt.figure(figsize=(20,20)).add_subplot(projection='3d')
ax.scatter3D(Xtransf[:,0], Xtransf[:,1], list(range(Xtransf.shape[0])))

for i,txt in enumerate(targetsub):
    ax.text(Xtransf[i,0],Xtransf[i,1],i,txt)

plt.title(f"{rat}{day} {rsl}, Isomap n=2, z=time", fontsize=18)

plt.xlabel("Dim 1",fontsize=16)

plt.ylabel("Dim 2", fontsize=16)

#%% Dim Reduction (tSNE) 2d 

rsl = "RS2"
rs = region_sets[rsl]
Xsub,targetsub = np.empty((0, len(population))), []
for r in rs:
    Xsub = np.vstack((Xsub, X[currentAlley==str(r),:]))
    targetsub.extend(traj[currentAlley==str(r)])
    
targetsub = np.asarray(targetsub)
n=2

Xtransf = TSNE(n_components=n).fit_transform(Xsub)
plt.scatter(Xtransf[:,0], Xtransf[:,1],s=1)

for i,txt in enumerate(targetsub):
    plt.annotate(txt,(Xtransf[i,0], Xtransf[i,1]))

# #color code the letters
# for i,txt in enumerate(targetsub):
#     for letter,offset in zip(txt,[0,0.2,0.4]):
#         c = {'N':'r','S':'b','E':'g','W':'y','X':'k'}[letter]
#         plt.annotate(letter,(Xtransf[i,0]+offset,Xtransf[i,1]),color=c)


    
plt.title(f"{rat}{day} {rsl}, tSNE n={n}", fontsize=18)

plt.xlabel("Dim 1",fontsize=16)

plt.ylabel("Dim 2", fontsize=16)

#%% Dimensionality Reduction 3d (tSNE) and visualize 

rsl = "RS2"
rs = region_sets[rsl]
Xsub,targetsub = np.empty((0, len(population))), []
for r in rs:
    Xsub = np.vstack((Xsub, X[currentAlley==str(r),:]))
    targetsub.extend(traj[currentAlley==str(r)])
    
targetsub = np.asarray(targetsub)
n=3

Xtransf = TSNE(n_components=n).fit_transform(Xsub)


ax = plt.figure(figsize=(20,20)).add_subplot(projection='3d')
ax.scatter3D(Xtransf[:,0], Xtransf[:,1], Xtransf[:,2])

for i,txt in enumerate(targetsub):
    ax.text(Xtransf[i,0],Xtransf[i,1], Xtransf[i,2],txt)

plt.title(f"{rat}{day} {rsl}, tSNE dim={n}", fontsize=18)
