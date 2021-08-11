# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 19:10:36 2021

@author: whockei1

Identifying Naturalistic Trajectories
Use schematic track to construct conditional transition probabilites
A trajectory is defined as two or more consecutive significantly biased
conditional transition probability tables for an alley

"""

#%% Imports 
import numpy as np
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import confounds as conf
from collections import Counter
from matplotlib.patches import Rectangle
import newAlleyBounds as nab
#%% Load data
rat, day = 'R859', 'D1'
population, turns = RepCore.loadRecordingSessionData(rat, day)  
ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859, 'R765':nab.R765, 'R886':nab.R886}[rat]


#%% Connectivity lookup and trajElement init
directConnections = {'0':np.array(['2','3','4']),
                     '4':np.array(['0','3','5','6']),
                     '6':np.array(['4','5','7']),
                     '2':np.array(['0','16']),
                     '3':np.array(['0','4','1','12','14']),
                     '5':np.array(['4','6','12','8','11']),
                     '7':np.array(['6','8','9']),
                     '1':np.array(['2','16','3','14','12']),
                     '12':np.array(['1','3','14','5','11','8']),
                     '8':np.array(['12','5','11','7','9']),
                     '16':np.array(['2','1','15']),
                     '14':np.array(['1','3','12','15','13']),
                     '11':np.array(['12','5','8','13','10']),
                     '9':np.array(['8','7','10']),
                     '15':np.array(['16','14','13']),
                     '13':np.array(['15','14','11','10']),
                     '10':np.array(['13','11','9'])
                     
                     }

trajElements = {str(i):None for i in range(17)}
#%% Create Conditional Transition Probability Table 

alley = '6'

alleyTurnIdx = np.where(turns['Alley+']==alley)[0]
turnsOnAlley = turns.iloc[alleyTurnIdx]

prevAlleys = list(set(directConnections[alley]).intersection(np.unique(turnsOnAlley['Alley-'])))



# Create conditional transition probability table
# Each row i (outer dict) is the previous alley and each column j (inner dict)
# is the next alley. The value i,j is the probability (over all turns coming
# from i) that the rat goes to j. 

conditionalTransitionProbabilities = {i:{j:0 for j in directConnections[alley]} for i in directConnections[alley]}

for pA in prevAlleys:
    prevTurns = turnsOnAlley[turnsOnAlley['Alley-']==pA]
    ptidx = np.array([i for i,_ in prevTurns.iterrows()])
    c = Counter(turns.iloc[ptidx+1]['Alley+'])
    
    for nA,count in c.items():
        if count >= 3:
            conditionalTransitionProbabilities[pA][nA] = count/sum(c.values())
    
cpt = pd.DataFrame.from_dict(conditionalTransitionProbabilities,orient='index')

#%% Identify single trajectory element 
# (i.e. a prev alley, current, next that occur frequently together. 
# Naturalistic trajectories will be made from finding multiple, overlapping such elements 

transThresh = 0.5 # 0.5 is a simple choice for now 8-11-21. much work will be needed to pick right thresh
ti = np.where(cpt>transThresh)

trajElements[alley] = list(zip(directConnections[alley][ti[0]], [alley]*len(ti[0]), directConnections[alley][ti[1]]))


#%% Visualize Trajectory

def drawRegion(ax, bounds, color, edge):
    """
    Bounds are the corners of a region on the track
    Format is [[xmin, xmax], [ymin, ymax]]
    Ax is an axis to which to add the region
    Color is the face color, edge is the edgecolor 
    """
    x0,y0 = bounds[0][0], bounds[1][0]
    w,h = bounds[0][1] - bounds[0][0], bounds[1][1] - bounds[1][0]
    ax.add_patch(Rectangle((x0,y0),w,h,facecolor=color, edgecolor=edge,linewidth=4))
    ax.autoscale_view() # for some reason the axes dont update automatically, so run this

ncol=2
fig, axes = plt.subplots(int(np.ceil(len(trajElements[alley])/ncol)),ncol,figsize=(10,10))
plt.suptitle(f"Trajectory Elements For Alley {alley}")
for i,trajectory in enumerate(trajElements[alley]):
    ax = fig.axes[i]
    ax.set_title(f"{trajectory[0]} -> {trajectory[1]} -> {trajectory[2]}")
    for a in range(17):
        
       if str(a) in trajectory: 
           if trajectory.index(str(a)) == 0:
               drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','green')
           elif trajectory.index(str(a)) == len(trajectory)-1:
               drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','red')
           else:
               drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','lightgrey')
    
    
       else:
           drawRegion(ax, ratborders.alleyInterBounds[str(a)],'lightgrey','lightgrey')
           
#%% prototype - Viewing trajectory elements as sequences if they exist
for a in range(17):
    a = str(a)
    te = trajElements[a]
    if te is not None:
        for el in te:
            for b in range(17):
                b = str(b)
                te_comp = trajElements[b]
                if te_comp is not None:
                    for el_comp in te_comp:
                        if el[1] == el_comp[0] and el[2] == el_comp[1]:
                            print(a,b,el,el_comp)
    

#%% prototype - Assembling trajectory elements into sequences if they exist

fullTrajectories = []
for a in range(17):
    a = str(a)
    te = trajElements[a]
    if te is not None:
        for el in te:
            for b in range(17):
                b = str(b)
                te_comp = trajElements[b]
                if te_comp is not None:
                    for el_comp in te_comp:
                        if el[1] == el_comp[0] and el[2] == el_comp[1]:
                            fullTrajectories.append(f"{el[0]}-{el[1]}-{el[2]}-{el_comp[2]}")
    
    
    
#%% Assemble trajectories iteratively. 

def assemble(trA, trB):
    """ Trajectories are built up by first taking two traj elements and iteratively
    checking to see if the interiors have overlapping alleys. check that here
    Assumes inputs will be string of alleys separated by '-'s """
    trlistA, trlistB = trA.split('-'), trB.split('-')
    lenoverlap = len(set(trlistA).intersection(set(trlistB)))    
    if lenoverlap == len(trlistA)-1:
        
        #logic is seeing which should go first. the earlier seq has its first 
        # element not in the common domain while the second seq has its last
        # element not in the common domain
        if trlistA[0] not in trlistB:
            newtraj = trlistA + [trlistB[-1]] #this works bc the new traj just adds the ends of the old ones to their overlap
        else:
            newtraj = trlistB + [trlistA[-1]]
        return '-'.join(newtraj)
    else:
        return None
   
###############################
#%% Putting Things Together ###
###############################

def drawRegion(ax, bounds, color, edge):
    """
    Bounds are the corners of a region on the track
    Format is [[xmin, xmax], [ymin, ymax]]
    Ax is an axis to which to add the region
    Color is the face color, edge is the edgecolor 
    """
    x0,y0 = bounds[0][0], bounds[1][0]
    w,h = bounds[0][1] - bounds[0][0], bounds[1][1] - bounds[1][0]
    ax.add_patch(Rectangle((x0,y0),w,h,facecolor=color, edgecolor=edge,linewidth=4))
    ax.autoscale_view() # for some reason the axes dont update automatically, so run this


def createCTPT(alley):
    alleyTurnIdx = np.where(turns['Alley+']==alley)[0]
    turnsOnAlley = turns.iloc[alleyTurnIdx]
    
    prevAlleys = list(set(directConnections[alley]).intersection(np.unique(turnsOnAlley['Alley-'])))
    
    
    
    # Create conditional transition probability table
    # Each row i (outer dict) is the previous alley and each column j (inner dict)
    # is the next alley. The value i,j is the probability (over all turns coming
    # from i) that the rat goes to j. 
    
    conditionalTransitionProbabilities = {i:{j:0 for j in directConnections[alley]} for i in directConnections[alley]}
    
    for pA in prevAlleys:
        prevTurns = turnsOnAlley[turnsOnAlley['Alley-']==pA]
        ptidx = np.array([i for i,_ in prevTurns.iterrows()])
        c = Counter(turns.iloc[ptidx+1]['Alley+'])
        
        for nA,count in c.items():
            if count >= 3:
                conditionalTransitionProbabilities[pA][nA] = count/sum(c.values())
        
    cpt = pd.DataFrame.from_dict(conditionalTransitionProbabilities,orient='index')
    
    return cpt


def findTrajElements(trajElements, alley, cpt, transThresh=0.5):
    
    transThresh = 0.5 # 0.5 is a simple choice for now 8-11-21. much work will be needed to pick right thresh
    ti = np.where(cpt>transThresh)
    trajElements[alley] = list(zip(directConnections[alley][ti[0]], [alley]*len(ti[0]), directConnections[alley][ti[1]]))

    return trajElements



def drawTrajElements(trajElements,alley):
    """Takes dict of trajectory elements and key, each value is len 3 list
    with prev,current,next alley of the t element. draws them on subplots
    annotated with green (start) and red (end)
    """
    ncol=2
    fig, axes = plt.subplots(int(np.ceil(len(trajElements[alley])/ncol)),ncol,figsize=(4,3.5))
    plt.suptitle(f"Trajectory Elements For Alley {alley}")
    for i,trajectory in enumerate(trajElements[alley]):
        ax = fig.axes[i]
        ax.set_title(f"{trajectory[0]} -> {trajectory[1]} -> {trajectory[2]}")
        for a in range(17):
            
           if str(a) in trajectory: 
               if trajectory.index(str(a)) == 0:
                   drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','green')
               elif trajectory.index(str(a)) == len(trajectory)-1:
                   drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','red')
               else:
                   drawRegion(ax, ratborders.alleyInterBounds[str(a)],'navy','lightgrey')
        
        
           else:
               drawRegion(ax, ratborders.alleyInterBounds[str(a)],'lightgrey','lightgrey')
               
               
               

for a in range(17):
    try:
        a = str(a)
        cpt = createCTPT(a)
        trajElements = findTrajElements(trajElements, a, cpt)
        drawTrajElements(trajElements,a)
    except:
        print(a)
