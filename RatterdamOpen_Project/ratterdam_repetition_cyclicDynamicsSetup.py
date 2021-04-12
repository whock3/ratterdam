# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 21:29:04 2021

@author: whockei1

Set up data for cyclic dynamics analysis (done by Vivek K.)

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
day = "D2"
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList, clustQuals = util.getClustList(df)
population = {}
qualThresh = 3

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

start, stop = 2912779894,6780860303
npoints = 60*10
interpTime = np.linspace(start, stop, npoints)
st = pd.Series(interpTime)

df_data = pd.DataFrame()
df_info = pd.DataFrame()


for clustname, unit in population.items():
    
    unit = removeFieldNans(unit)
    
    pchip_fields = [pchip(field[:,0], field[:,1]) for field in unit.fields]
       
    interpFields = []
    for pc in pchip_fields:
        interpFields.append(pc(interpTime)[30:-30])
            
    
    field_series_mult = [pd.Series(f) for f in interpFields]
    
    for i,ser in enumerate(field_series_mult):
        df_data = df_data.append(ser,ignore_index=True)
        df_info = df_info.append(pd.Series([clustname,i]),ignore_index=True)
        

            
df_data.columns = st[30:-30]
df_info.columns = ["Neuron", "Field"]

df = pd.concat([df_info, df_data], axis=1)
df.to_csv("R859D2_spline.csv",header=True,index=False)

        
        
        
        