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
import newAlleyBounds as nab
import math
import bisect
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
#%% Setup
rat = 'R781'
day = 'D3'
ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859}[rat]
savepath = "E:\\Ratterdam\\repetition_decoding\\"
datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(datapath)
population = {}
qualThresh = 3

for i,clust in enumerate(clustList):
    
    if clustQuals[i] >= qualThresh:
   
        print(clust)
        unit = RepCore.loadRepeatingUnit(datapath, clust, smoothing=1)                                   
        rm = util.makeRM(unit.spikes, unit.position)
        if np.nanpercentile(rm, 95) > 1.:
            population[clust] = unit
            print(f"{clust} included")
        else:
            print(f"{clust} is not included")
        
        
# Session endpoints data
with open(datapath+"sessionEpochInfo.txt","r") as f:
    lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    
nepoch=3
intervals = np.linspace(start,end,nepoch+1)
pos, turns = alleyTrans.alleyTransitions(unit.position, ratborders, graph=False)
turns = pd.DataFrame(turns)
turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']

turns = pd.DataFrame(data=turns)
turns.dropna(inplace=True) 

cmap = util.makeCustomColormap()


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



    


#%% Extract and plot each trajectory (2d and schematic alleys visited) color-coded by rate 
savepath = f'E:\\Ratterdam\\repetition_decoding\\{rat}{day}_trajectories\\'
window=1


timestamp = util.genTimestamp()

codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
for unitname, unit in population.items():
        
    u = unitname.split('\\')[0]+unitname.split('\\')[1]
    print(u)
    
    with PdfPages(f"{savepath}_{timestamp}_{u}_Trajectories.pdf") as pdf:

        
        for f,field in enumerate(unit.fields):
            allregions, alltrajs, alllabels, allrates = [], [], [], []
            _vmax = max([i[1] for i in field])
            norm = mpl.colors.Normalize(vmin=0,vmax=_vmax)
            for i,visit in enumerate(field):
                            
                turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
                turndata = [np.nan]
                if turnIdx > 1 and turnIdx < turns.shape[0]-2:
                    regions = []
                    dirs = []
                    trajturns = turns.iloc[turnIdx-window:turnIdx+window+1]
                    ts_start, ts_end = float(turns.iloc[turnIdx-window-1]['Ts exit']), float(turns.iloc[turnIdx+window+1]['Ts exit'])
                    behavtraj = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end),1:]
                    behavtraj = behavtraj[(behavtraj[:,0]>0)&(behavtraj[:,1]>0)]
                    
                   
                    it=0
                    for _, t in trajturns.iterrows():
                        regions.append(t['Alley-']) 
                        regions.append(t['Inter'])
                        regions.append(t['Alley+'])
                        dirs.append(codedict[t['Allo-']])
                        if it == trajturns.shape[0]-1:                   
                            dirs.append(codedict[t['Allo+']])
                        regions = list(set(regions)) # there is redundancy based on overlapping def of turns, so get unique regions
                        it+=1 # pd.iterrows() gives the absolute row number in the overall df so do a manual increment for the size of the sub-df were using here
                    
                    
                    
                    allregions.append(regions)
                    subplotlabel = f"{','.join(dirs)}"
                    alllabels.append(subplotlabel)
                    allrates.append(visit[1])
                    alltrajs.append(behavtraj)
            allregions = np.asarray(allregions)
            alllabels = np.asarray(alllabels)
            allrates = np.asarray(allrates)
            alltrajs = np.asarray(alltrajs)
            
            plotOrder = np.flip(np.argsort(allrates))
                    

            ncols=10 
            fig, ax = plt.subplots(int(np.ceil(len(field)/ncols)),ncols,figsize=(15,15))
    
            for i,pid in enumerate(plotOrder):           
    
                for region in allregions[pid]:
                    drawRegion(fig.axes[i],ratborders.alleyInterBounds[region],cmap(norm(allrates[pid])))
                    
                fig.axes[i].plot(unit.perimeters[f][:,0], unit.perimeters[f][:,1],linestyle='--',linewidth=3)
                fig.axes[i].plot(alltrajs[pid][:,0],alltrajs[pid][:,1],color='black',linewidth=2)
                fig.axes[i].scatter(alltrajs[pid][-1,0],alltrajs[pid][-1,1],color='lime',marker='^',s=40,zorder=99)
                fig.axes[i].scatter(alltrajs[pid][0,0],alltrajs[pid][0,1],color='red',marker='o',s=40,zorder=99)
                
                
                fig.axes[i].set_title(str(pid)+" "+alllabels[pid], fontsize=12)
                
            plt.suptitle(f"Unit {unitname} Field {f}, Max rate {round(_vmax,3)} Hz")
            plt.subplots_adjust(left=0.01, right=0.89,top=0.95,bottom=0.05,wspace=0.25,hspace=0.25)
            for ii in range(len(fig.axes)):
                fig.axes[ii].axis("off")
                fig.axes[ii].set_aspect('equal', adjustable='box')
            cax = fig.add_axes([0.9,0.2,0.01,0.5])
            cbar = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
            cbar.set_label("Hz")
    
    
                                    
            try:
                pdf.savefig()
                print("Saved")
            except:
                print(f"Could not save {unitname} field {f}, moving on...")
            plt.close()
                        
                
        
        
        
        
        
    











