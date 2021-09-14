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
import repeatingPC as repPC
import copy
#%% Setup

for rat,day in zip(['R781', 'R781', 'R808', 'R808', 'R859', 'R859', 'R886', 'R886', 'R765','R765'],['D3', 'D4', 'D6', 'D7', 'D1', 'D2', 'D1', 'D2','RFD5','DFD8']):

    try:
        ratborders = nab.loadAlleyBounds(rat, day)
        savepath = f"E:\\Ratterdam\\raw_data_for_JK\\{rat}{day}_trajectories\\"
        datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
        
        population, turns = RepCore.loadRecordingSessionData(rat, day)
        
        ballisticTurnIdx = []
        for i in range(1,turns.shape[0]-1):
           row = turns.iloc[i]
           inter = row['Inter']
           if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter and turns.iloc[i-1].Inter != inter:
               ballisticTurnIdx.append(i)
        
        refturns = copy.deepcopy(turns) # keep a copy without filtering.
        turns = turns.iloc[np.asarray(ballisticTurnIdx)]
        
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
        
        
        def drawTrack(ax=None):
            if ax == None:
                ax = plt.gca()
            for i in range(17):
                drawRegion(ax, ratborders.alleyInterBounds[str(i)],'lightgrey')
                  
                    
        
        
        #%% Extract and plot each trajectory (2d and schematic alleys visited) color-coded by rate 
        window=1
        
        
        timestamp = util.genTimestamp()
        
        codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
        for unitname, unit in population.items():
                
            u = unitname.split('\\')[0]+unitname.split('\\')[1]
            print(u)
            repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)

            
            with PdfPages(f"{savepath}{timestamp}_{u}_Trajectories.pdf") as pdf:
                
                fig, ax = plt.subplots(figsize=(8,6))
                ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                              cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                       extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
                ax.set_title(f"{unitname}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
                ax.axis('equal')
                ax.set_ylim([0,480])
                ax.set_xlim([0, 640])
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                
                # fields over time (+ a line to add field borders to overall ratemap above)
                # 6/17/21 this code is adapted from the ratterdam_RepetitionCore.plotRoutine... fx but
                # have removed 1) lines dealing w if field dynamic traces are smoothed or not 2) time vs visits on xaxis
                # and 3) whether youre manually removing bad fields at end bc that issue was partially solved
                for i, field in enumerate(unit.smoothedFields):       
                       ax.plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i], label=f"Field {i}")
                ax.legend()
                
                
                pdf.savefig()
                plt.close()
                
                for f,field in enumerate(unit.fields):
                    
                    foverlap = overlaps[f]
                    
                    allregions, alltrajs, alllabels, allrates = [], [], [], []
                    _vmax = max([i[1] for i in field])
                    norm = mpl.colors.Normalize(vmin=0,vmax=_vmax)
                    for i,visit in enumerate(field):
                                    
                        turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
                        turndata = [np.nan]
                        if turnIdx >= 1 and turnIdx < turns.shape[0]-1:
                                                
                            regions = []
                            dirs = []
                            trajturns = turns.iloc[turnIdx-1:turnIdx+1+1] # want to be explicit: +1 to get next turn and another +1 bc python isnt right inclusive
                            
                            if trajturns.iloc[0]['Alley+'] == trajturns.iloc[1]['Alley-'] and trajturns.iloc[1]['Alley+'] == trajturns.iloc[2]['Alley-']:
                            
                                if int(trajturns.iloc[1]['Alley+']) in foverlap and int(trajturns.iloc[1]['Alley-']) in foverlap:
                                    predir, postdir = trajturns.iloc[0]['Alley-'], trajturns.iloc[2]['Alley+']
                                    ts_start, ts_end = float(refturns.iloc[turnIdx-1-1]['Ts entry']), float(refturns.iloc[turnIdx+1+1+1]['Ts exit'])
                                    regions = [trajturns.iloc[0]['Alley-'], 
                                               trajturns.iloc[1]['Alley-'], 
                                               trajturns.iloc[1]['Alley+'], 
                                               trajturns.iloc[2]['Alley+'],
                                               trajturns.iloc[0]['Inter'],
                                               trajturns.iloc[1]['Inter'],
                                               trajturns.iloc[2]['Inter']
                                               ]
                                    dirs = [trajturns.iloc[0]['Allo-'], trajturns.iloc[0]['Allo+'], trajturns.iloc[1]['Allo+'], trajturns.iloc[2]['Allo+']]
    
                                elif int(trajturns.iloc[1]['Alley+']) in foverlap:
                                    predir, postdir = trajturns.iloc[1]['Alley-'], trajturns.iloc[2]['Alley+']
                                    ts_start, ts_end = float(trajturns.iloc[0]['Ts entry']), float(refturns.iloc[turnIdx+1+1+1]['Ts exit'])
                                    regions = [trajturns.iloc[1]['Alley-'], 
                                               trajturns.iloc[1]['Alley+'], 
                                               trajturns.iloc[2]['Alley+'],
                                               trajturns.iloc[1]['Inter'],
                                               trajturns.iloc[2]['Inter']
                                               ]
                                    dirs = [trajturns.iloc[1]['Allo-'], trajturns.iloc[1]['Allo+'], trajturns.iloc[2]['Allo+']]
    
                                elif int(trajturns.iloc[1]['Alley-']) in foverlap:
                                    predir, postdir = trajturns.iloc[0]['Alley-'], trajturns.iloc[1]['Alley+']
                                    ts_start, ts_end = float(refturns.iloc[turnIdx-1-1]['Ts entry']), float(trajturns.iloc[1]['Ts exit'])
                                    regions = [trajturns.iloc[0]['Alley-'], 
                                               trajturns.iloc[1]['Alley-'], 
                                               trajturns.iloc[1]['Alley+'],
                                               trajturns.iloc[0]['Inter'],
                                               trajturns.iloc[1]['Inter']
                                               ]
                                    dirs = [trajturns.iloc[0]['Allo-'], trajturns.iloc[1]['Allo-'], trajturns.iloc[1]['Allo+']] 
                  
                                elif trajturns.iloc[1]['Inter'] in foverlap:
                                    predir, postdir = trajturns.iloc[1]['Alley-'], trajturns.iloc[1]['Alley+']
                                    ts_start, ts_end = float(trajturns.iloc[0]['Ts entry']), float(trajturns.iloc[2]['Ts exit'])
                                    regions = [trajturns.iloc[1]['Alley-'], 
                                               trajturns.iloc[1]['Alley+'], 
                                               trajturns.iloc[1]['Inter']]
                                    dirs = [trajturns.iloc[1]['Allo-'], trajturns.iloc[1]['Allo+']]
                            
                                else:
                                    print("Error - no field overlap detected")
                                
                                behavtraj = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end),1:]
                                dirs = [codedict[i] for i in dirs]
                                
                                # gets rid of -99s (lost data)
                                behavtraj = behavtraj[(behavtraj[:,0]>0)&(behavtraj[:,1]>0)]
               
                     
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
                        print(f"Field {f} Saved")
                    except:
                        print(f"Could not save {unitname} field {f}, moving on...")
                    plt.close()
                                
                        
    except:
        print(f"Error with {rat} {day} turn visualization")
            
            
            
    











