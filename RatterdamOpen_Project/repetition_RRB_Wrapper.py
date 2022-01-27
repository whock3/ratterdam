# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 17:12:04 2022

@author: whockei1

Wrapper for repetition_Directionality_RegionRepeatingBreakdown.py

Interface to filter data in different ways and compare the results 
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt, copy
import repetition_Directionality_RegionRepeatingBreakdown as RRB

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
interdf = pd.read_csv(interdatapath)


#%% Helper fx 

def convert_toggles_to_string(toggles):
    """toggles is a dict"""
    
    
    params = []
    for togKey in ['repeating_type', 'region_type', 'region_location']:
        params.append(toggles[togKey])
    
    label = ','.join(params)
    return label
        

#%% 

all_rep = RRB.calculate_subgroup_directionality(alleydf, interdf, 
                                                toggles = {
                                                'repeating_type':'R',
                                                'region_type':'All',
                                                'region_location':'All'
                                                            }
                                                )


all_sf = RRB.calculate_subgroup_directionality(alleydf, interdf, 
                                                toggles = {
                                                'repeating_type':'SF',
                                                'region_type':'All',
                                                'region_location':'All'
                                                            }
                                                )



#%% 

toggle_list = [
    {
    'repeating_type':'All',
    'region_type':'All',
    'region_location':'I'
                        },
    {
    'repeating_type':'All',
    'region_type':'All',
    'region_location':'P'
                        },
        


      
     
    ]

labels = []
diffs = []

for tog in toggle_list:
    
    diffs.append(RRB.calculate_subgroup_directionality(alleydf, interdf, tog))
    labels.append(convert_toggles_to_string(tog))
    
    
fig, ax = plt.subplots()
ax.violinplot(diffs, positions=range(len(diffs)))
ax.set_title("Mean Directional Difference Across Fields", fontsize=20)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=45,fontsize=14)
ax.set_ylabel("Mean Directional Firing Rate (Hz)", fontsize=14)
    