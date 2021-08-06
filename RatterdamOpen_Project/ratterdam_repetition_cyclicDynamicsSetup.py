# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 21:29:04 2021

@author: whockei1

Set up data for cyclic dynamics analysis Ivan A and earlier work done by Vivek K.

"""
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
import math
import bisect
import pandas as pd
from scipy.interpolate import PchipInterpolator as pchip


#%% define fx 

def removeFieldNans(unit):
    for i,field in enumerate(unit.fields):
         fieldnn = field[~np.isnan(field).any(axis=1)]
         unit.fields[i] = fieldnn
    return unit

#%% read in data
    
rat = "R859"
day = "D1"
savepath = f'E:\\Ratterdam\\leader_follower_data\\{rat}\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(df)
population = {}
qualThresh = 3

print(clustList)

for i,clust in enumerate(clustList):
    
    if clustQuals[i] >= qualThresh:
        
        try:
            print(clust)
            unit = RepCore.loadRepeatingUnit(df, clust, smoothing=1)
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
        except:
            pass


#%% Create dataframe
df_data = pd.DataFrame()
df_info = pd.DataFrame()


# 8-6-21: before today, code would loop over fields and use pchip interpolation
# to interpolate fields and equate all the sample numbers across fields.
# Baryshinkov group actually needs fields to be raw and they deal with uneven data
# sizes. This issue was noted a long time ago but they have been using one day R859D2 
# which i either created the csv for manualy or did in a jupyter notebook and dont feel
#like finding since its pretty trivial. 

for clustname, unit in population.items():  
    
    for i,field in enumerate(unit.fields):
        df_data = df_data.append(pd.Series(field[:,0]),ignore_index=True)
        df_info = df_info.append(pd.Series([clustname,f"field {i}",'timestamps']),ignore_index=True)
        
        df_data = df_data.append(pd.Series(field[:,1]),ignore_index=True)
        df_info = df_info.append(pd.Series([clustname,f"field {i}",'firingrates']),ignore_index=True)
        
            
df = pd.concat([df_info, df_data], axis=1)
df.to_csv(savepath+f"{rat}{day}_spline.csv",header=True,index=False)

        
        
        
        