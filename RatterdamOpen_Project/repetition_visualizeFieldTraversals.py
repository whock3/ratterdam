# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 16:08:31 2021

@author: whockei1

Revised script to plot each trajectory through the field and the associated spikes
Including turn right before and after field traversal
Revision is based on defining field traversal as a (potentially) multi-turn event
with (potentially) multiple different headings. Here, the different headings thru
field (if there are >1) are not segmented, just plotted w spikes for viewer to analyze
visually.

Code is based on repetition_revisedTurnAssignment.py to get traversals
and repetition_visualizeDir to organize and visualize (but this was buggy due to how traversals were defined)
"""


#%% Imports
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
from matplotlib import path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import repeatingPC as repPC
import copy
import more_itertools
import pickle


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
              
                    
with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
        
for rat, day in zip(['R765'],['RFD5']):
    df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    savepath = f"E:\\Ratterdam\\raw_data_for_JK\\{rat}{day}_trajectories\\"
    
    population, turns, refturns = superpop[rat][day]['units'], superpop[rat][day]['turns'], superpop[rat][day]['refturns']
    rewards = RepCore.readinRewards(rat, day)
    ratborders = nab.loadAlleyBounds(rat, day)
    codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
    cmap = util.makeCustomColormap()
    
    # Remove turnarounds/pivots
    # ballisticTurnIdx = []
    # for i in range(1,turns.shape[0]-1):
    #    row = turns.iloc[i]
    #    inter = row['Inter']
    #    if row['Ego'] != '3' and turns.iloc[i+1].Inter != inter and turns.iloc[i-1].Inter != inter:
    #        ballisticTurnIdx.append(i)
    
    # refturns = copy.deepcopy(turns) # keep a copy without filtering.
    # turns = turns.iloc[np.asarray(ballisticTurnIdx)]
    
    #%% Run routine
    
    timestamp = util.genTimestamp()
    
    for clutName, unit in population.items():    
        u = unit.name.split('\\')[0]+unit.name.split('\\')[1]
        print(u)
        try:
            repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
        
            with PdfPages(f"{savepath}{timestamp}_{u}_Trajectories.pdf") as pdf:
                
                fig, ax = plt.subplots(figsize=(10,8))
                ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                              cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                       extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
                ax.set_title(f"{unit.name}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
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
            
                for fnum, (field, perim, foverlap) in enumerate(zip(unit.fields, unit.perimeters, overlaps)):
                
                    
                    contour = path.Path(perim) # get a path around the field which we will use to filter spikes in-field latert
                    
                    foverlap = [str(i) for i in foverlap]
                    # Field size data, used to filter visits according to threshold pct'ages of field size
                    maxx,minx, maxy, miny = max(perim[:,0]), min(perim[:,0]), max(perim[:,1]), min(perim[:,1])
                    dista = np.sqrt(((minx-minx)**2)+((maxy-miny)**2))
                    distb = np.sqrt(((minx-maxx)**2)+((miny-miny)**2))
                    
                    # Iterate over turns and assign to fields, get rates and labels
                    
                    turnNums = []
                    
                    for tnum, turn in refturns.iterrows():
                        # trim the ends bc were looking ahead and behind as we go and dont want an out of index error
                        if tnum > 2 and tnum < refturns.shape[0]-2:
                            if turn['Alley-'] in foverlap or turn['Alley+'] in foverlap or turn['Inter'] in foverlap:
                                
                                # Get behavior and spikes using refturns which has all the turns bc we need the flanking turns
                                # because here we are looking for turn into/out of the field. so +/-2 gets us the next/prev alleys 
                                ts_start, ts_end = float(refturns.iloc[tnum-2]['Ts entry']), float(refturns.iloc[tnum+2]['Ts exit'])
                                behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
                                behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
                                s = unit.spikes[(unit.spikes[:,0]>ts_start)&(unit.spikes[:,0]<=ts_end)]
                                
                                filtOutcome = RepCore.filterVisit(dista,distb,behav,perim,length_thresh=0.3,dist_thresh=0.1,dist_point_thresh=3,inside_point_thresh=3)
                                if filtOutcome == True:
                                    turnNums.append(tnum)
            
            
                    turnNums = np.asarray(turnNums)
                    
                        
                    # Goal here is to take all turns involving a field and group them by
                    # the axis of travel.
                    
                    # group turns into sets of consecutive turns, each group being the
                    # turns that involve the field for a given visit to that field
                    visitTurns = [list(group) for group in more_itertools.consecutive_groups(turnNums)]
                    
            
                        
                    ###############################################
                    # Loop over traversals and get data for field #
                    ###############################################
                       
                    traversalRates = []
                    traversalTrajectories = []
                    traversalSpikes = []
                    traversalRegions = []
                    traversalRewards = []
                    traversalTurnNums = []
                    
                    for group in visitTurns:
                    
                        visitData = []
                        regionData = []
                        
                        
                                
                        # get start, end timestamps for each alley along the traversal through the field
                        # The alleys must overlap w the field, i.e we are not getting the alleys pre/post traversal
                        for idx in group:
                            turn = refturns.iloc[idx]
                            if turn['Alley-'] in foverlap:
                                if len(visitData)>=1 and refturns.iloc[idx-2]['Ts entry'] not in visitData[-1]:
                                    visitData.append((refturns.iloc[idx-2]['Ts entry'], turn['Ts exit']))
                                    
                                    # add the track regions involved so we can plot them later
                                    _ = [regionData.append(i) for i in  [turn['Alley-'], 
                                                                         turn['Inter'],
                                                                         turn['Alley+'],
                                                                         refturns.iloc[idx-2]['Alley+'],
                                                                         refturns.iloc[idx-2]['Inter'],
                                                                         refturns.iloc[idx-1]['Inter']
                                                                         ]
                                                                         ]
                                elif visitData == []:
                                    visitData.append((refturns.iloc[idx-2]['Ts entry'], turn['Ts exit']))
                                    
                                    # add the track regions involved so we can plot them later
                                    _ = [regionData.append(i) for i in  [turn['Alley-'], 
                                                                         turn['Inter'],
                                                                         turn['Alley+'],
                                                                         refturns.iloc[idx-2]['Alley+'],
                                                                         refturns.iloc[idx-2]['Inter'],
                                                                         refturns.iloc[idx-1]['Inter']
                                                                         ]
                                                                         ]
                                    
                            if turn['Alley+'] in foverlap:
                                if len(visitData)>=1 and refturns.iloc[idx-1]['Ts entry'] not in visitData[-1]:
                                    visitData.append((refturns.iloc[idx-1]['Ts entry'], refturns.iloc[idx+2]['Ts exit']))
                                    
                                    # add the track regions involved so we can plot them later
                                    _ = [regionData.append(i) for i in  [turn['Alley-'], 
                                                                         turn['Inter'],
                                                                         turn['Alley+'],
                                                                         refturns.iloc[idx+2]['Alley-'],
                                                                         refturns.iloc[idx+2]['Inter'],
                                                                         refturns.iloc[idx+1]['Inter']
                                                                         ]
                                                                         ]
                                    
                                elif visitData == []:
                                    visitData.append((refturns.iloc[idx-1]['Ts entry'], refturns.iloc[idx+2]['Ts exit']))
                                    _ = [regionData.append(i) for i in  [turn['Alley-'], 
                                                                         turn['Inter'],
                                                                         turn['Alley+'],
                                                                         refturns.iloc[idx+2]['Alley-'],
                                                                         refturns.iloc[idx+2]['Inter'],
                                                                         refturns.iloc[idx+1]['Inter']
                                                                         ]
                                                                         ]
                                    
                            # check passes through inter when field doesnt overlap involved alleys. Wont have double counting here. 
                            if turn['Alley-'] not in foverlap and turn['Alley+'] not in foverlap and turn['Inter'] in foverlap:
                                # intersection traversals with a 90deg turn are ambiguous as to their direction without more involved analysis. 
                                if turn['Allo-'] == turn['Allo+']:
                                    visitData.append((refturns.iloc[idx-1]['Ts entry'], refturns.iloc[idx+1]['Ts exit']))
                                    _ = [regionData.append(i) for i in  [turn['Alley-'], 
                                                                         turn['Inter'],
                                                                         turn['Alley+']
                                                                         ]
                                                                         ]
                        #visitData might be empty if there is a right-angle turn that only involes an intersection
                        # that overlaps the field. 
                        if visitData != []:
                        
                            regionData = np.unique(regionData) # we double count often above, so get unique track regions
                            alltimes = []
                            for v in visitData:
                                for i in v:
                                    alltimes.append(i)
                                    
                            alltimes = sorted([float(i) for i in alltimes])
                            
                                        
                            ts_start, ts_end = alltimes[0], alltimes[-1]
                            spikes = unit.spikes[(unit.spikes[:,0]>ts_start)&(unit.spikes[:,0]<=ts_end)]
                            spikes = spikes[contour.contains_points(spikes[:,1:])]
                            trajectory = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
                            trajectory = trajectory[(trajectory[:,1]>0)&(trajectory[:,2]>0)&(trajectory[:,1]<640)&(trajectory[:,2]<480)]
                            rate = spikes.shape[0]/((ts_end-ts_start)/1e6)
                            #rewards is a list of reward times, so youre generating bool list of whether reward i is in the current time window
                            isReward = np.where(np.asarray([(ts_start < i < ts_end) for i in rewards])==True)[0]
                            if isReward.shape[0]>0:
                                traversalRewards.append(True)
                            else:
                                traversalRewards.append(False)
                            
                            traversalRates.append(rate)
                            traversalTrajectories.append(trajectory)
                            traversalSpikes.append(spikes)
                            traversalRegions.append(regionData)
                            traversalTurnNums.append(group[0])
                        
                    ##############################################
                    # Now plot traversals in order of their rate #
                    ##############################################
                    
                    plotOrder = np.flip(np.argsort(traversalRates))
                    numPlots = len(plotOrder)
                    norm = mpl.colors.Normalize(vmin=0,vmax=max(traversalRates))
                    plotsPerPage = 25
                    nrows = 5
                    plotIntervals = np.ceil(np.linspace(0,numPlots-1,int(np.ceil(numPlots/plotsPerPage))+1)).astype(int)
                    
                    for pI in range(len(plotIntervals)-1):
                        pStart, pEnd = plotIntervals[pI], plotIntervals[pI+1]
                        fig, ax = plt.subplots(nrows, int(np.ceil(plotsPerPage/nrows)),figsize=(12,10))
                        
                        for aix, pid in enumerate(plotOrder[pStart:pEnd]):
                        
                            trajectory = traversalTrajectories[pid]
                            spikes = traversalSpikes[pid]
                            rate = traversalRates[pid]
                            regions = traversalRegions[pid]
                        
                            drawTrack(ax=fig.axes[aix])
                            for region in regions:
                                drawRegion(ax=fig.axes[aix], bounds=ratborders.alleyInterBounds[region], color=cmap(norm(traversalRates[pid])))
                            
                            fig.axes[aix].plot(trajectory[:,1], trajectory[:,2], color='k')
                            fig.axes[aix].scatter(trajectory[0,1], trajectory[0,2],color='k',marker='o',s=50,zorder=99)
                            fig.axes[aix].scatter(trajectory[-1,1], trajectory[-1,2], color='k',marker='^',s=50,zorder=99)
                            fig.axes[aix].scatter(spikes[:,1], spikes[:,2], marker='+',c='g',s=25, zorder=99)
                            fig.axes[aix].set_title(f"{pid},{traversalTurnNums[pid]},{round(rate,2)}Hz,{int(traversalRewards[pid])}",fontsize=8)
                            fig.axes[aix].plot(perim[:,0],perim[:,1],c='k',linestyle='--')
                                
                            
                        plt.suptitle(f"{rat}{day} {unit.name} Field {fnum} (p{pI+1})")
                        for ii in range(len(fig.axes)):
                                        fig.axes[ii].axis("off")
                                       
                        cax = fig.add_axes([0.9,0.2,0.01,0.5])
                        cbar = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
                        cbar.set_label("Hz")
                        
                        pdf.savefig()
                        plt.close()
        except:
            print(f"Error with {u}")
        