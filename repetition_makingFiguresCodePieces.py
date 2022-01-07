# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 00:45:58 2021

@author: whockei1
"""

import ratterdam_RepetitionCoreFx as RepCore
import ratterdam_Defaults as Def
import utility_fx as util 




#%% 

ratborders = nab.loadAlleyBounds(rat, day)
letters = [i for i in ascii_uppercase[:12]]
fig, ax = plt.subplots()
txt = ['S','E','N']
selected_alleys = [3,12,5]
selected_inters = ['F','G']
selected_color = 'navy'
unselected_color = 'grey'
txtcolor = 'white'
for i in range(17):
    if i in selected_alleys:
        idx = selected_alleys.index(i)
        util.drawRegion(ax, ratborders.alleyInterBounds[str(i)],selected_color,txt=txt[idx],txtcolor=txtcolor) 
    else:
        util.drawRegion(ax, ratborders.alleyInterBounds[str(i)],unselected_color)
for j in letters:
    if j in selected_inters:
        util.drawRegion(ax, ratborders.alleyInterBounds[j],selected_color)
    else:
        util.drawRegion(ax, ratborders.alleyInterBounds[j],unselected_color)




#%% GLM LRT grouped bars 
import pandas as pd, numpy as np, matplotlib.pyplot as plt 
# based on alleypath <- "E:\\Ratterdam\\R_data_repetition\\211210_AlleySuperpopDirVisitFiltered.csv"
# rewards removed, turnarounds excluded. 30% overlap. 
# filtOutcome = RepCore.filterVisit(dista,distb,behav,perim,
# length_thresh=0.2,
# dist_thresh=0.1,
# dist_point_thresh=2,
# inside_point_thresh=2)
# 3 passes each dir
rep = [31/128, 18/128, 9/128, 2/128]
nr_all = [18/75, 6/75, 6/75, 4/75]
nr_multi = [12/45, 6/45, 4/45, 2/45]
nr_single = [6/30, 0/30, 2/30, 2/30]
rows = ["+Current Dir", "+Previous Dir", "+Next Dir", "All Dirs"]

data = {'Repeating':rep,
        "All Non-repeating":nr_all,
        "Multi-field Non-repeating":nr_multi,
        "Single field Non-repeating":nr_single
       }

df = pd.DataFrame(data)

ax = df.plot(kind='bar',
        color=['red','navy','cornflowerblue','lightblue'],
        fontsize=42
        )

ax.set_xticks([0,1,2,3])
ax.set_xticklabels(rows,rotation=0)
ax.set_ylabel("Proportion of Fields",fontsize=54)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.hlines(0.05,-0.5,4,linestyle='--',color='k',linewidth=5)
plt.subplots_adjust(bottom=0.2)



#%% Examples of temporal dynamics in repeating fields
# Plotting fields on separate subplots (RM just taken from a pdf)

rat, day  ='R859', 'D3'
clust = 'TT3\\cl-maze1.3'

unit = RepCore.loadRepeatingUnit(rat, day, clust)