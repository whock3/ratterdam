# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 17:44:53 2021

@author: whockei1

Fig 5 - Noncurrent direction coding
Sfn2021 poster and paper draft
"""
#%% imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ratterdam_Defaults as Def
from scipy.stats import mannwhitneyu

#%% Fig 5b - bars of best models by LRT (5a is a schematic of the testing)

# taken from r, proportion of fields tested with best model being:
#,c,p,n, pcn
replrts = [0.238,0.137,0.119,0.045]
nonreplrts = [0.238,0.119,0.071,0.047]
labels = ['+C','+C+P','+C+N','+P+C+N']

fig, _ax = plt.subplots()
ax = fig.axes[0]
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

rects1 = ax.bar(x - width/2, replrts, width,
                label='Repeating Fields',
                facecolor='red',
                edgecolor='darkred',
                linewidth=2)
rects2 = ax.bar(x + width/2, nonreplrts, width,
                label='Non-repeating Fields',
                facecolor='dodgerblue',
                edgecolor='navy',
                linewidth=2)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Proportion of Fields \n in Each Group \nBest Fit By Model',fontsize=Def.ylabelsize)
ax.set_xticks(x)
ax.set_xticklabels(labels,rotation=0)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.hlines(0.05,-0.5,len(labels),linestyle='--',color='k',linewidth=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(prop={'size':30})
plt.subplots_adjust(left=0.2)


#%% 5c - violins of rmse changes

df = pd.read_csv("E:\\Ratterdam\\2021_SfNPoster_WH\\Fig5_NoncurrentDirectionality\\lrt_models_RMSE.csv")

fig, ax = plt.subplots()
v = ax.violinplot([df.m2_rmse-df.m1_rmse,
               df.m3_rmse-df.m2_rmse,
               df.m4_rmse-df.m2_rmse,
               df.m5_rmse-((df.m3_rmse+df.m4_rmse)/2)],
              positions=range(4),
              )
for vpart in v['bodies']:
    vpart.set_facecolor('grey')
    vpart.set_edgecolor('black')
    vpart.set_linewidth(2)
    
for el in ['cbars','cmaxes','cmins']:
    v[el].set_color('black')
               
ax.set_xticks(range(4))
ax.set_xticklabels(['+C','+C+P','+C+N','P+C+N'],rotation=0)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.set_ylabel("Change in RMSE Compared to \nNext Simplest Model",fontsize=Def.ylabelsize)  
ax.hlines(0,-0.5,3.5,color='k',linestyle='--',linewidth=2) 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False) 
plt.subplots_adjust(left=0.15)


#%% 5d - violins of aic changes

df = pd.read_csv("E:\\Ratterdam\\2021_SfNPoster_WH\\Fig5_NoncurrentDirectionality\\lrt_models_aic.csv")
df.dropna(inplace=True)
fig, ax = plt.subplots()
v = ax.violinplot([df.m2_aic-df.m1_aic,
               df.m3_aic-df.m2_aic,
               df.m4_aic-df.m2_aic,
               df.m5_aic-((df.m3_aic+df.m4_aic)/2)],
              positions=range(4),
              )

for vpart in v['bodies']:
    vpart.set_facecolor('grey')
    vpart.set_edgecolor('black')
    vpart.set_linewidth(2)
    
for el in ['cbars','cmaxes','cmins']:
    v[el].set_color('black')
               
ax.set_xticks(range(4))
ax.set_xticklabels(['+C','+C+P','+C+N','P+C+N'],rotation=0)
ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.set_ylabel("Change in AIC Compared to \nNext Simplest Model",fontsize=Def.ylabelsize)  
ax.hlines(0,-0.5,3.5,color='k',linestyle='--',linewidth=2) 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False) 
plt.subplots_adjust(left=0.15)
