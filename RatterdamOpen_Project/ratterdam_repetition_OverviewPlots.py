# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 13:59:01 2021

@author: whockei1

Ratterdam Repetition Project

Overview graphs - 1 PDF per cell
page 1-
page 2-
page 3-

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



#%% Load Data and 
rat = "R781"
day = "D4"
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\overviewPlots\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(df)
population = {}
qualThresh = 3
date = util.genTimestamp()
cmap = util.makeCustomColormap()

for i,clust in enumerate(clustList):
    
    if clustQuals[i] >= qualThresh:
   
        try:
            print(clust)
            unit = RepCore.loadRepeatingUnit(df, clust, smoothing=2)
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
        except:
            pass
        
#%% Make pdf for each unit

pagewidth = 12

for unitname, unit in population.items():
    
    print(f"Creating overview plots for Rat {rat} day {day} cell {unitname}...")
    
    unitname = unitname.split('\\')[0]+unitname.split('\\')[1]
    with PdfPages(savepath+date+"_"+unitname+"overview.pdf") as pdf:
        
        # page 1 - ratemap, field dynamics over time, IFD plots   
        
        p1gs = gridspec.GridSpec(3,2)
        p1gs.update(wspace=0.025, hspace=0.3) 
        page1fig = plt.figure(figsize=(pagewidth,15))
        p1ax1, p1ax2, p1ax3, p1ax4 = page1fig.add_subplot(p1gs[0,:]), page1fig.add_subplot(p1gs[1,:]), page1fig.add_subplot(p1gs[2,1]), page1fig.add_subplot(p1gs[2,0])
        
        # overall ratemap
        p1ax1.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                       cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
        p1ax1.set_title(f"{clust}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
        p1ax1.axis('equal')
        p1ax1.set_ylim([0,480])
        p1ax1.set_xlim([0, 640])
        p1ax1.set_xticks([])
        p1ax1.set_yticks([])
        p1ax1.spines['top'].set_visible(False)
        p1ax1.spines['right'].set_visible(False)
        p1ax1.spines['bottom'].set_visible(False)
        p1ax1.spines['left'].set_visible(False)
        
        # fields over time (+ a line to add field borders to overall ratemap above)
        # 6/17/21 this code is adapted from the ratterdam_RepetitionCore.plotRoutine... fx but
        # have removed 1) lines dealing w if field dynamic traces are smoothed or not 2) time vs visits on xaxis
        # and 3) whether youre manually removing bad fields at end bc that issue was partially solved
        for i, field in enumerate(unit.smoothedFields):
            
            #xvals in time since session beginning, converted to mins
            xval = field[:,0]
            xval = [((i-field[0,0])/1e6)/60 for i in xval]
            
            p1ax2.plot(xval, field[:,1], color=unit.colors[i], marker='.',alpha=0.8, label=i)
            p1ax1.plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i])
            p1ax2.text(xval[0]-0.1,field[0,1]-0.1,i)
            p1ax2.tick_params(axis='y', labelsize=14)
            p1ax2.tick_params(axis='x', labelsize=14)
            p1ax2.set_xlabel("Time in session (mins)", fontsize=12)
            p1ax2.set_ylabel(f"Firing Rate Hz (sigma = {unit.smoothing})", fontsize=14)
            p1ax2.spines['right'].set_visible(False)
            p1ax2.spines['top'].set_visible(False)
            p1ax2.set_title("Place Field Dynamics", fontsize=20)
        p1ax2.legend()
                  
            
        # IFD plots and autocorr of them  if >1 field
        
        if len(unit.fields)>1:
            
            diffmats, wins = RepCore.makeSemaphores(unit.fields)
            dm_vmax = np.percentile([np.nanmax(dm.flatten()) for dm in diffmats],95)
            
            
            #nested gridspec is a bit tedious, see: https://stackoverflow.com/questions/28933233/embedding-multiple-gridspec-layouts-on-a-single-matplotlib-figure
            nrow = int(np.ceil(len(diffmats)/7))
            ncol = 7
            inner_grid = gridspec.GridSpecFromSubplotSpec(nrow,ncol, subplot_spec=p1gs[2,1])
            
            c = 0# super hacky but get counters of how many r/c pairs were thru to index scalar diffmat in list
            for row in range(nrow):
                for col in range(ncol):   
                    try:
                        if row == 0 and col == 0:
                            ax = plt.subplot(inner_grid[row,col])
                            ax.set_title(f"IFD Matrices, 95%ile = {round(dm_vmax,2)} Hz")
                            im = ax.imshow(diffmats[c],interpolation='None',vmax=dm_vmax,aspect='auto',cmap=cmap)
                            c += 1
                        else:
                            ax = plt.subplot(inner_grid[row,col])
                            ax.imshow(diffmats[c],interpolation='None',vmax=dm_vmax,aspect='auto',cmap=cmap)                        
                            c += 1  
                            
                        #idk why i needed to pull this out of the if/else, otherwise all subplots are the first diffmat.
                        ax.set_xticks([])
                        ax.set_yticks([])
                        ax.set_xticklabels([])
                        ax.setyticklabels([])
                    except:
                        pass
        
                        
            # Autocorr of IFD   
            corrdm = RepCore.corrOfMats(diffmats,plot=False,ret=True)
            p1ax4.imshow(corrdm, aspect='auto',interpolation='None',origin='lower')
            p1ax4.set_title("Autocorrelation of IFD Matrices")
        
        
        try:
            pdf.savefig()
        except:
            print("Error with page 1, continuing on without saving it")
        plt.close()
        
        
        #page 2 - directionality
        
        
        spos,sspikes = conf.dirFiltWindow(unit.position,unit.spikes,1,10)
        fig,ax = plt.subplots(2,2,figsize=(pagewidth,12))
        
            
        allrmsflat, allrms = [],[]
        for i, (p,s) in enumerate(zip(spos,sspikes)):
            rm = util.makeRM(s,p,bins=[30,50])
            allrmsflat.extend(rm.flatten())
            allrms.append(rm)
            vmax = np.nanpercentile(allrmsflat, 99)
            
        for i, (r,d) in enumerate(zip(allrms,['N','E','S','W'])):
            fig.axes[i].imshow(r,aspect='auto',interpolation='None',origin='lower',cmap=cmap,vmax=vmax)
            fig.axes[i].set_title(d)
            fig.axes[i].set_xticks([])
            fig.axes[i].set_yticks([])
            fig.axes[i].axis("equal") 
            plt.suptitle(f"98%ile {round(vmax,2)}Hz")
            
        try:
            pdf.savefig()
        except:
            print("Error with page 2, continuing on without saving it")
        plt.close()
        
    

    
    
    
    
    

