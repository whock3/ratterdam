# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 13:36:04 2021

@author: whockei1

Aggregating data based on turns to visualize intuitively any relation between
firing and direction/turn/trajectory.

"""

#%% Imports
import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
import alleyTransitions as alleyTrans
import math
import bisect
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
#%% Setup
pos, turns = alleyTrans.alleyTransitions(unit.position, ratborders, graph=False)
turns = pd.DataFrame(turns)
turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']

turns = pd.DataFrame(data=turns)
turns.dropna(inplace=True) 


#%% Helper fx 

def drawRegion(ax, bounds,color):
    """
    Bounds are the corners of a region on the track
    Format is [[xmin, xmax], [ymin, ymax]]
    Ax is an axis to which to add the region
    Color can be anything in practice it will be the rate 
    """
    x0,y0 = bounds[0][0], bounds[1][0]
    w,h = bounds[0][1] - bounds[0][0], bounds[1][1] - bounds[1][0]
    ax.add_patch(Rectangle((x0,y0),w,h,color=color))
    ax.autoscale_view() # for some reason the axes dont update automatically, so run this


#%% Turn based rate maps - setup
# Group data by allo dir pre (4 plots) and allo dir post (4 plots) +/- 1.5 turns
#(the half turn is bc the turns end at the intersection so punch out a bit more
# in time to get the full +/1 1 turn)
#Reminder about code order: N,E,S,W,  F,R,B,L


turnRMS = {'Pre':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
                  }, 
           
           'Post':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
                  '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
                  }
           }

for field in unit.fields:
    for visit in field:
        turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
        try:
            # first/last turn has no previous/next turn, so just ignore it
            if turnIdx > 2 and turnIdx < turns.shape[0]-2:
                spikesPre = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
                spikesPost = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

                turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'],spikesPre))
                turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'],spikesPost))
    
                occPre = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
                occPost = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

                turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'],occPre))
                turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'],occPost))
        except:
            print("Likely invalid code found")

for d in ['Pre', 'Post']:
    fig, ax = plt.subplots(2,2)
    for i,code in enumerate([('1','N'),('2','E'),('3','S'),('4','W')]):
        s,o = turnRMS[d][code[0]]['Spikes'], turnRMS[d][code[0]]['Pos']
        rm = util.makeRM(s,o)
        fig.axes[i].imshow(rm,aspect='auto',interpolation='None',cmap=cmap,vmax=7,origin='lower')
        fig.axes[i].set_title(code[1])
    plt.suptitle(f"{d}-turn Bearing")
    
    


#%% Extract and plot each trajectory (2d and schematic alleys visited) color-coded by rate 

savepath = 'E:\\Ratterdam\\repetition_decoding\\R859D2_trajectories\\'
ncols=10
window=1

codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
for unitname, unit in population.items():
    u = unitname.split('\\')[0]+unitname.split('\\')[1]
    with PdfPages(savepath+f"{u}_fieldTrajectories.pdf") as pdf:
        
        for f,field in enumerate(unit.fields):
            _vmax = max([i[1] for i in field])
            fig, ax = plt.subplots(int(np.ceil(len(field)/ncols)),ncols,figsize=(15,12))
            norm = mpl.colors.Normalize(vmin=0,vmax=_vmax)
            for i,visit in enumerate(field):
                
                turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
                turndata = [np.nan]
                if turnIdx > 1 and turnIdx < turns.shape[0]-2:
                    regions = []
                    dirs = []
                    trajturns = turns.iloc[turnIdx-window:turnIdx+window]
                    ts_start, ts_end = float(trajturns.iloc[0]['Ts exit']), float(trajturns.iloc[-1]['Ts exit'])
                    behavtraj = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end),1:]
                    behavtraj = behavtraj[(behavtraj[:,0]>0)&(behavtraj[:,1]>0)]
            
                    for t in trajturns.iterrows():
                        regions.append(t[1]['Alley-']) # iterrows gives a tuple and the second entry is the pd series we want.
                        regions.append(t[1]['Inter'])
                        regions.append(t[1]['Alley+'])
                        dirs.append(codedict[t[1]['Allo-']])
                        dirs.append(codedict[t[1]['Allo+']])
                        regions = list(set(regions)) # there is redundancy based on overlapping def of turns, so get unique regions
                        
                        
                    for region in regions:
                        drawRegion(fig.axes[i],ratborders.alleyInterBounds[region],cmap(norm(visit[1])))
                    
                    fig.axes[i].plot(unit.perimeters[f][:,0], unit.perimeters[f][:,1],linestyle='--')
                    fig.axes[i].plot(behavtraj[:,0],behavtraj[:,1],color='k')
                    fig.axes[i].scatter(behavtraj[-1,0],behavtraj[-1,1],color='green',marker='^',s=40,zorder=99)
                    fig.axes[i].axis("off")
                    fig.axes[i].set_title(f"{','.join(dirs)}", fontsize=9)
                    plt.suptitle(f"Unit {unitname} Field {f}, Max rate {_vmax}")
                            
            try:
                pdf.savefig()
            except:
                print(f"Could not save {unitname} field {f}, moving on...")
            plt.close()
                    
            
    
    
    
    
    












