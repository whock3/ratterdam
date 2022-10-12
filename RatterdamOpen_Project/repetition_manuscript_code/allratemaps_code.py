# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:49:52 2022

@author: whockei1

All ratemaps figure 
Self-explanatory, have a ratemap for each unit

"""
#%%
%matplotlib qt5

import matplotlib.pyplot as plt, numpy as np
import utility_fx as util, williamDefaults as wmDef
import pickle


with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  

#%%
    
# should be 178 neurons 
cmap = util.makeCustomColormap()
r,c = 10,5
fig, ax = plt.subplots(r,c,figsize=(20,20))
i=0
j=0
for rat in superpop.keys():
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():

            if i%(r*c) == 0:

                for _ax in fig.axes:
                    _ax.set_xticks([])
                    _ax.set_yticks([])
                    for spineloc in ['top','bottom','left','right']:
                        _ax.spines[spineloc].set_visible(False)

                fig, ax = plt.subplots(r,c,figsize=(20,20))
                fig.subplots_adjust(wspace=0, hspace=0.2)

                j=0


            ax = fig.axes[j]
            ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                       cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])

            if len(unit.fields) > 0:
                color = 'black'
            else:
                color = 'red'
            ax.set_title(f"{rat}{day} {unit.name}",fontsize=8,y=0.9,color=color)
            for spineloc in ['top','bottom','left','right']:
                ax.spines[spineloc].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect('equal')
            
            i += 1
            j += 1
            
# plt.subplots_adjust(left=0.1,
#                     right=0.9,
#                     bottom=0.15,
#                     top=0.8,
#                     wspace=0.1,
#                     hspace=0.1   
#     )


#%% Getting specific neurons (i.e. SF from one day) for presentation purposes

cmap = util.makeCustomColormap()
fig, ax = plt.subplots(5,3,figsize=(25,20))
i=0
rat, day = 'R859','D2'

allSpikes = np.empty((0,3))
for unit in superpop[rat][day]['units'].values():
    if len(unit.fields)==1:
        allSpikes = np.vstack((allSpikes, unit.spikes))
        ax = fig.axes[i]
        ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                   cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
            extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
        ax.set_title(f"{rat}{day} {unit.name}",fontsize=5,y=0.9)
        for spineloc in ['top','bottom','left','right']:
            ax.spines[spineloc].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        
        i += 1
            
            

fig, ax= plt.subplots()
rm = util.makeRM(allSpikes,unit.position)
plt.imshow(rm,aspect='auto',interpolation='None',origin='lower',cmap=cmap,vmax=10)