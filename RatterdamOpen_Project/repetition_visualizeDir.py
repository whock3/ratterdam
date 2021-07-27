# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 16:26:06 2021

@author: whockei1

Script to visualize directionality within CA1 neurons
1. Look at whole population across each alley grouped by direction. Direction
is defined using travel through schematic track, not pt-by-pt movement
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

#%% Load data, define defaults


for rat,day in zip(['R859','R859','R781','R781','R808','R808'],['D1','D2','D3','D4','D6','D7']):
    print(rat,day)
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    codedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
    ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859}[rat]
    
    
    #%% Create data arrays by region
    # Each region has an (v,n) shaped array where v are the visits to the alley
    # and n are the neurons. Each row is a pop response vector, each row is a unit's
    # response across trials. 
    
    alleyArrays = {str(i):{'neural_data':np.empty((0,len(population))),'labels':np.empty((0))} for i in range(17)}
    
    for t in range(1,turns.shape[0]-1):
        turn_nm1 = turns.iloc[t-1]
        turn = turns.iloc[t]
        turn_np1 = turns.iloc[t+1]
                
        cD = codedict[turn['Allo+']]
        
        if cD != 'X':
        
            start, stop = float(turn['Ts entry']), float(turn_np1['Ts exit'])
            duration = (stop-start)/1e6
            
            popvector = []
            for i, (unitname, unit) in enumerate(population.items()):
                #basic stuff, get the firing rate (# spikes / time on alley ) for each unit and append
                spike_count = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)].shape[0]
                rate = spike_count / duration
                popvector.append(rate)
        
            popvector = np.asarray(popvector)
            alleyArrays[turn['Alley+']]['neural_data'] = np.vstack((alleyArrays[turn['Alley+']]['neural_data'], popvector))
            alleyArrays[turn['Alley+']]['labels'] = np.append(alleyArrays[turn['Alley+']]['labels'], cD)
        
        
    #%% Calc rate by direction for each cell in region and visualize
    timestamp = util.genTimestamp()
    fig,ax = plt.subplots(5,4, sharey=False, figsize=(17,15))
    
    for i,alley in enumerate([0,4,6,2,3,5,7,1,12,8,16,14,11,9,15,13,10]):
        
        if alley in [0,4,6,2,16,15,13,10,9,7]:
            loc = 'Perim'
        else:
            loc = 'Interior'
        
        alley = str(alley)
    
        X,Y = alleyArrays[alley]['neural_data'], alleyArrays[alley]['labels']
        dirs = sorted(np.unique(Y))
        if len(dirs) == 2:
            Xa, Xb = X[Y==dirs[0],:], X[Y==dirs[1],:] # assuming two unique directions on current alley, filter whole array X into 2 subarrays
            if Xa.shape[0] >=3 and Xb.shape[0] >=3:

                width=0.4
                fig.axes[i].bar(np.arange(Xa.shape[1])-(width/2), Xa.mean(axis=0), width, yerr=Xa.std(axis=0)/np.sqrt(Xa.shape[0]),color='red')
                fig.axes[i].bar(np.arange(Xb.shape[1])+(width/2), Xb.mean(axis=0), width, yerr=Xb.std(axis=0)/np.sqrt(Xb.shape[0]),color='blue')
     
                fig.axes[i].set_title(alley + f" {loc}", fontsize=14)
            else:
                print(f"Insufficient sampling for {rat}{day} alley ")
                fig.axes[i].set_title(f"{alley} {loc} Insufficient Sampling")
        else:
            print(f"Number unique dirs = {len(dirs)} for {rat}{day} {alley}, so not running")
            fig.axes[i].set_title(f"{alley} {loc} Total bias", fontsize=14)
    
    plt.suptitle(f"{rat}{day} Population Directionality by Alley. Mean +/- SEM", fontsize=20)
    plt.subplots_adjust(wspace=0.25,hspace=0.75)
    plt.savefig("E:\\Ratterdam\\repetition_decoding\\raw_directionality\\"+f"{timestamp}_{rat}{day}_directionByAlley.png",dpi=300)
    plt.close()
        
        
    










