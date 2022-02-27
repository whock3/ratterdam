# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:08:14 2022

@author: whockei1

Repetition project behavioral characterization
Purpose is to quantify and visualize things like running speed thru alley
to put in the "basic characterization" figure (Fig 1) in repetition manuscript 


"""
import numpy as np
import utility_fx as util
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import pandas as pdo
import ratterdam_DataFiltering as Filt 
import utility_fx as util 
import newAlleyBounds as nab

rat, day  = 'R765', 'DFD4'
turns, unit = RepCore.loadTurns(rat,day)
turns = turns[turns.Ego!='0']

ratborders = nab.loadAlleyBounds(rat, day)
rewards = RepCore.readinRewards(rat, day)




verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]

speeds = []
speedbins = np.linspace(0,60,13)
alleybins = np.linspace(0,15,15)
for tnum, turn in turns.iterrows():
    if tnum < turns.shape[0]-1:
    
        # Because turn definition is staggered as it moves across track,
        # pick either alley+ or alley- by convention and make the ts, labels
        # match that choice 
        alley = turn['Alley+']
        ts_start, ts_end = float(turn['Ts entry']), float(turns.iloc[tnum+1]['Ts exit'])
        behav = unit.position[(unit.position[:,0]>ts_start)&(unit.position[:,0]<=ts_end)]
        behav = behav[(behav[:,1]>0)&(behav[:,2]>0)]
        
        if alley in verticals:
            lintraj = behav[:,np.r_[0,2]] # ts,y
        elif alley in horizontals:
            lintraj = behav[:,:-1] # ts, x
           
        
        diff = np.diff(lintraj,axis=0)
        diff[:,0] = diff[:,0]/1e6
        diff[:,1] = diff[:,1]/Def.ptsCm_macaulay
        speed = np.asarray([abs((lintraj[1:,1]-lintraj[0,1])/Def.ptsCm_macaulay), abs(diff[:,1]/diff[:,0])])
        h = np.histogram2d(speed[1,:],speed[0,:],bins=[speedbins,alleybins])
        speeds.append(h[0]/h[0].sum())
            
        
       