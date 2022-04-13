# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 12:08:34 2022

@author: whockei1

Repetition - Independence of repeating field directionality

Ran 2-way ANOVAs in repetition_RepeatingFieldDirectionality.py 
Spent a while in fall 2021 figuring out best processing pipeline for it

Now wrapping that in script here and doing further analysis to explore result


As of 2/10/22 using this file for Figure 4 - will copy into a file
called Fig 4... when I think the analysis is stable. 
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
import itertools
from scipy.stats import linregress 


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)


def fd_anova(ocell):
    """
    Runs a two-way anova with following model:
    'Rate ~ C(FieldNum)*C(CurrDir)'
    We want to know if there's an interaction between firing rate and field
    ie. field tuning is not shared across fields 
    
    Input is orientation-filtered field.
    We check if there's > 1 field inside
    
    Returns
    -------
    p value of interaction term in model.
    or None if orientation-filtered cell only has 1 (filtered) field

    """
    if np.unique(ocell.FieldID).shape[0] > 1:
    
        mod = ols('Rate ~ C(FieldNum)*C(CurrDir)',data=ocell).fit()
        aov = sm.stats.anova_lm(mod,typ=2)
        pint = aov['PR(>F)']['C(FieldNum):C(CurrDir)']
    
    else:
        pint = None
        
    return pint

#%% 
x = []
y = []
sigs = []
prods = []
all_diffs = []

for cellid, cell in alleydf.groupby("CellID"):
    for oname, ocell in cell.groupby("Orientation"):
        
        if np.unique(ocell.FieldID).shape[0] > 1:
        
            pint = fd_anova(ocell)
            signed_diffs = []
            
            if pint != None and pint < 0.025:
                sig = True
            else:
                sig = False
            sigs.append(sig)
            for fid, field in ocell.groupby("FieldID"):
            
                dirs = sorted(np.unique(field.CurrDir)) # sort alphabetically to keep signs consistent below 
                if len(dirs) > 2:
                    print(f"dir error {fid} {oname}")
                dirA = field[field.CurrDir==dirs[0]]
                dirB = field[field.CurrDir==dirs[1]]
                
              
                #this fails if all trials are 0
                try:
                    diff = (dirA.Rate.mean()-dirB.Rate.mean())/field.Rate.mean()
                except:
                    diff = 0
                
                signed_diffs.append(diff)
                all_diffs.append(diff)
                
            # get pairwise comparisons and products
            combs = itertools.combinations(signed_diffs,2)
            for pair in combs:
                i,j = pair
                    
                x.append(i)
                y.append(j)
                
                prods.append(i*j)
                        
#%% Plot scatter signed dir A vs B

from scipy.stats import linregress
import matplotlib.patches as patches

reload(MDef)

colors = ['grey' if s==False else 'red' for s in sigs]


fig, ax = plt.subplots()
ax.scatter(x,y, s=400,c='grey',edgecolor='k',linewidth=1)
ax.set_ylabel("Signed Normalized\n Directionality Field A", fontsize=MDef.ylabelsize)
ax.set_xlabel("Signed Normalized\n Directionality Field B", fontsize=MDef.xlabelsize)
ax.set_title("Pairwise Directionality Comparison of Aligned Repeating Fields \nwithin each Repeating Neuron",
             fontsize=MDef.titlesize)

xlim, ylim = ax.get_xlim(), ax.get_ylim()


ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.hlines(0, xlim[0], xlim[1], linestyle='--',color='k',linewidth=2)
ax.vlines(0, ylim[0], ylim[1], linestyle='--',color='k',linewidth=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)

## run ols and plot line 
olsfit = linregress(x,y)
xfit = np.linspace(xlim[0], xlim[1],100)
yfit = [(olsfit.slope*i)+olsfit.intercept for i in xfit]
ax.plot(xfit,yfit,color='black',linewidth=3)
ax.fill_between(xfit,yfit-olsfit.stderr,yfit+olsfit.stderr,color='grey',alpha=0.7)
plt.text(2,2.5,f"$R^2$ = {round(olsfit.rvalue**2,3)}",fontsize=32) # this position is hardcoded, if data changes a lot need to move

xlim, ylim = ax.get_xlim(), ax.get_ylim()

alpha = 0.3

widths = [xlim[1]-xlim[0],
     xlim[0],
     xlim[0],
     xlim[1]-xlim[0]
     ]

heights = [ylim[1]-ylim[0],
     ylim[1]-ylim[0],
     ylim[0],
     ylim[0]  
     ]

colors = ['red','blue','red','blue']
     
     
for w,h,c in zip(widths, heights, colors):    
    rect = patches.Rectangle((0,0),w,h,
                             facecolor=c,
                             alpha=alpha,
                             zorder=0
                            )

    ax.add_patch(rect)

#%%

# plot product 
fig, ax = plt.subplots()
#ax.scatter(prods,range(len(prods)),s=80,c='grey',edgecolor='black',linewidth=2)
ax.hist(prods,bins=50,color='grey',edgecolor='k',linewidth=3)
xlim, ylim = ax.get_xlim(), ax.get_ylim()
ax.vlines(0, ylim[0], ylim[1], linestyle='--',color='k',linewidth=5)
ax.set_ylabel("Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Directionality Product", fontsize=MDef.xlabelsize)
ax.set_title("Product of Signed Directionality of Aligned Repeating Fields \nwithin each Repeating Neuron",
             fontsize=MDef.titlesize)

alpha = 0.3

ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)


xlim, ylim = ax.get_xlim(), ax.get_ylim()
widths = [xlim[0],
          xlim[1]-xlim[0]]
  
heights = [ylim[1],
           ylim[1]]  

colors = ['blue','red']
        
for w,h,c in zip(widths, heights, colors):    
    rect = patches.Rectangle((0,0),w,h,
                             facecolor=c,
                             alpha=alpha,
                             zorder=0
                            )

    ax.add_patch(rect)
    
    
#%% Look at pairwise directional tuning relationships only within common
# 'arms' of the maze. I.e. the same linear segment running along each axis. 4V, 3H

maze_arm_organization = {'V':{'V1':[2,16],
                         'V2':[3,14],
                         'V3':[5,11],
                         'V4':[7,9],
                         },
                         'H':{
                         'H1':[0,4,6],
                         'H2':[1,12,8],
                         'H3':[15,13,10]
                         }
                         }


x = []
y = []
sigs = []
prods = []
field_arm_info = []


for cellid, cell in alleydf.groupby("CellID"):
    for oname, ocell in cell.groupby("Orientation"):
        
        orientation_structure = maze_arm_organization[oname]
        
        if np.unique(ocell.FieldID).shape[0] > 1:
        
            signed_diffs = []
            alleys = []
            fids = []

            for fid, field in ocell.groupby("FieldID"):
            
                dirs = sorted(np.unique(field.CurrDir)) # sort alphabetically to keep signs consistent below 
                if len(dirs) > 2:
                    print(f"dir error {fid} {oname}")
                dirA = field[field.CurrDir==dirs[0]]
                dirB = field[field.CurrDir==dirs[1]]
                #this fails if all trials are 0
                try:
                    diff = (dirA.Rate.mean()-dirB.Rate.mean())/field.Rate.mean()
                except:
                    diff = 0
                
                signed_diffs.append(diff)
                alleys.append(np.unique(field.Alleys)[0])
                fids.append(fid)
                
            combs = itertools.combinations(range(len(signed_diffs)),2)
            
            for pair in combs:
                i,j = pair
                alleyA, alleyB = alleys[i], alleys[j]
                
                # this is clumsy but getting each arm each field is on
                # assumes one and only one hit when checking orientation_structure 
                for arm, armAlleys in orientation_structure.items():
                    if alleyA in armAlleys:
                        armA = arm        
                for arm, armAlleys in orientation_structure.items():
                    if alleyB in armAlleys:
                        armB = arm
                if armA != armB:
                    diffA, diffB = signed_diffs[i], signed_diffs[j]
                    x.append(diffA)
                    y.append(diffB)
                    
                    fidA, fidB = fids[i], fids[j]
                    sign = np.sign(diffA*diffB)
                    field_arm_info.append([getFieldInfoFromID(fidA,alleydf), getFieldInfoFromID(fidB, alleydf),sign])
                    
            
#%% plot scatter for x,y. Same code as above, but up there it's all
# pairwise comparisons within repeating cell. here its pairwise comparisons
# within "arm" of maze

from scipy.stats import linregress
import matplotlib.patches as patches

reload(MDef)

fig, ax = plt.subplots()
ax.scatter(x,y, s=400,c='grey',edgecolor='k',linewidth=1)
ax.set_ylabel("Signed Normalized\n Directionality Field A", fontsize=MDef.ylabelsize)
ax.set_xlabel("Signed Normalized\n Directionality Field B", fontsize=MDef.xlabelsize)
ax.set_title("Pairwise Directionality Comparison of Aligned Repeating Fields \nwithin Different Arms within Repeating Neuron",
             fontsize=MDef.titlesize)

xlim, ylim = ax.get_xlim(), ax.get_ylim()


ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.hlines(0, xlim[0], xlim[1], linestyle='--',color='k',linewidth=2)
ax.vlines(0, ylim[0], ylim[1], linestyle='--',color='k',linewidth=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)

## run ols and plot line 
olsfit = linregress(x,y)
xfit = np.linspace(xlim[0], xlim[1],100)
yfit = [(olsfit.slope*i)+olsfit.intercept for i in xfit]
ax.plot(xfit,yfit,color='black',linewidth=3)
ax.fill_between(xfit,yfit-olsfit.stderr,yfit+olsfit.stderr,color='grey',alpha=0.7)
plt.text(1,2.5,f"$R^2$ = {round(olsfit.rvalue**2,3)}",fontsize=32) # this position is hardcoded, if data changes a lot need to move

xlim, ylim = ax.get_xlim(), ax.get_ylim()

alpha = 0.3

widths = [xlim[1]-xlim[0],
     xlim[0],
     xlim[0],
     xlim[1]-xlim[0]
     ]

heights = [ylim[1]-ylim[0],
     ylim[1]-ylim[0],
     ylim[0],
     ylim[0]  
     ]

colors = ['red','blue','red','blue']
     
     
for w,h,c in zip(widths, heights, colors):    
    rect = patches.Rectangle((0,0),w,h,
                             facecolor=c,
                             alpha=alpha,
                             zorder=0
                            )

    ax.add_patch(rect)