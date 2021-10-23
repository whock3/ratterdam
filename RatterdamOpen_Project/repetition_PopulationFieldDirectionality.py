# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 14:37:03 2021

@author: whockei1

Script to analyze field directionality for each field in a recording day
Additionally annotated with whether cell is repeating or not
Will also visualize the results
"""

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
import repeatingPC as repPC
import math
import bisect
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib import path
import copy
import traceback
import numpy.ma as ma
from scipy.stats import mannwhitneyu


datadir = 'E:\\Ratterdam\\R_data_repetition\\'
data, replabels, pvalues = [], [], []

metadata = [('R765','RFD5','20210916-151948_R765RFD5_segmentedDirections_1.5vfilt.csv'),
            ('R781','D3','20210903-102901_R781D3_segmentedDirections_1.5vfilt.csv'),
            ('R781','D4','20210916-150730_R781D4_segmentedDirections_1.5vfilt.csv'),
            ('R808','D6','20210903-113101_R808D6_segmentedDirections_1.5vfilt.csv'),
            ('R808','D7','20210903-113320_R808D7_segmentedDirections_1.5vfilt.csv'),
            ('R859','D1','20210902-170828_R859D1_segmentedDirections_1.5vfilt.csv'),
            ('R859','D2','20210903-112710_R859D2_segmentedDirections_1.5vfilt.csv'),
            ('R886','D1','20210916-151110_R886D1_segmentedDirections_1.5vfilt.csv'),
            ('R886','D2','20210916-151358_R886D2_segmentedDirections_1.5vfilt.csv')
            ]


for rat, day, fname in metadata:

    df = pd.read_csv(datadir+fname)
    ratborders = nab.loadAlleyBounds(rat, day)
    
    population, turns = RepCore.loadRecordingSessionData(rat, day)
    unitList = np.unique(df['Unit'])
    
    #%% Compute average field directional difference 
    for unitname in unitList:
        
        unit = population[unitname]
        repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)
    
        
        fields = np.unique(df[df['Unit']==unitname]['Field'])
        
        for f in fields:
            
            replabels.append(repeat)
            
            subdf = df[(df['Unit']==unitname)&(df['Field']==f)]
            directions = np.unique(subdf['Direction'])
            
            if 'E' in directions and 'W' in directions and 'N' in directions and 'S' in directions:
                bf = fields.shape[0] + 2
            else:
                bf = fields.shape[0]
            
            if 'E' in directions and 'W' in directions:           
                horizMeanDiff = np.abs(subdf[subdf['Direction']=='E']['Rate'].mean()-
                                       subdf[subdf['Direction']=='W']['Rate'].mean())   
                xval = horizMeanDiff
                
                # Mann whitney will error out if all values are the same. 
                # obviously in that case there is no sig difference so manually code pval=1
                try:
                    xpval = mannwhitneyu(subdf[subdf['Direction']=='E']['Rate'],
                                         subdf[subdf['Direction']=='W']['Rate']).pvalue/bf
                except:
                    xpval = 1
            else:
                xval = 0
                xpval = 1 # code as arbitrary nonsig
                
            if 'N' in directions and 'S' in directions:
                vertMeanDiff = np.abs(subdf[subdf['Direction']=='N']['Rate'].mean()-
                                       subdf[subdf['Direction']=='S']['Rate'].mean())        
                yval = vertMeanDiff
                try:
                    ypval = mannwhitneyu(subdf[subdf['Direction']=='N']['Rate'],
                                         subdf[subdf['Direction']=='S']['Rate']).pvalue/bf
                except:
                    ypval = 1
            
            else:
                yval = 0
                ypval = 1
                
                
            data.append([xval,yval])
            pvalues.append([xpval, ypval])
        
data = np.asarray(data)
pvalues = np.asarray(pvalues)