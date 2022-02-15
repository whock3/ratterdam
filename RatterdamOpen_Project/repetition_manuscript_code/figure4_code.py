# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:33:35 2022

@author: whockei1

Figure 4 code - repeating field directionality. These fields are not all identically
tuned within the same repeating cell

2 example units
scatter of pairwise differences between fields. expressed as signed difference in FR
hist of product of the pairwise signed differences 

"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
from scipy.stats import sem
import repetition_manuscript_defaults as MDef
from importlib import reload


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)


#%% Example units. Plot ratemap and for each field plot bars of average firing
# for each direction +/- SEM


cmap = util.makeCustomColormap()
example_cell_list = [["R781","D4","TT3\\cl-maze1.2","V"]
                     
                     ]

for celldata in example_cell_list:
    cdf = alleydf[(alleydf['Rat']==celldata[0])
                  &(alleydf['Day']==celldata[1])
                  &(alleydf['CellName']==celldata[2])
                  &(alleydf['Orientation']==celldata[3])            
                  ]
    fieldgroups = cdf.groupby("FieldNum")
    
    ncols=3
    fig, ax_ = plt.subplots(int(np.ceil(len(fieldgroups)/ncols)+1),ncols, figsize=(8,8))
    for i, (fname, fg) in enumerate(fieldgroups):
         i=i+1 # shift plots, want rm in first 
         dgroups = fg.groupby("CurrDir")['Rate'].mean()
         bar = dgroups.plot(kind='bar',
                      yerr=fg.groupby("CurrDir")['Rate'].sem(),
                      ax=fig.axes[i],
                      fontsize=24,
                      edgecolor='k',
                      linewidth=2
                      )
         bar.set_ylabel("Mean Rate +/- SEM",fontsize=36)
         bar.set_title(f"Field {i-1}",fontsize=36)
         bar.spines['right'].set_visible(False)
         bar.spines['top'].set_visible(False)
         bar.spines['left'].set_linewidth(MDef.spine_width)
         bar.spines['bottom'].set_linewidth(MDef.spine_width)

    rat = celldata[0]
    day = celldata[1]
    clustName= celldata[2]
        
    unit = RepCore.loadRepeatingUnit(rat, day, clustName)
    ax = fig.axes[0]
    ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                  cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
           extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
    ax.set_title(f"{unit.name}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
    ax.axis('equal')
    ax.set_ylim([0,480])
    ax.set_xlim([0, 640])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for i, field in enumerate(unit.smoothedFields):       
        ax.plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i], label=f"Field {i}")
    lgnd = ax.legend(prop={'size':26})
    for lhand in lgnd.legendHandles:
        lhand._legmarker.set_linewidth(5)
    lgnd.get_frame().set_linewidth(MDef.legend_frame_width)
    
    clustName = clustName.replace("\\","_")
    
    #rescale ylim of bar graphs, cant share y bc ratemap is on different scale
    ymax = max([ax.get_ylim() for ax in fig.axes][1:])[1]
    for ax in fig.axes[1:]:
        ax.set_ylim([0,ymax])
    