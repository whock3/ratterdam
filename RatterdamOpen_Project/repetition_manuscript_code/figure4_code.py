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
import itertools
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
from scipy.stats import sem
import repetition_manuscript_defaults as MDef
from importlib import reload


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
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
    
    
    
#%% 
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
    
#%% Run anovas 
datapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
df = pd.read_csv(datapath)
totalAnovas = 0
passThresh = 2
fdirIntsAlley = []
fdirAlley = []
orientationDirs = {'V':['N','S'],'H':['E','W']}
for cellid in df['CellID'].unique():
    u = df[df['CellID']==cellid]
    
    oriens = np.unique(u.Orientation)
    if oriens.shape[0] == 1:
        padj = 0.05
    elif oriens.shape[0] == 2:
        padj = 0.025
    else:
        print("error wrong number orientation")
        
    for orien in oriens:
        uOriented = u[u['Orientation']==orien]
        # get num passes by dir per field
        pint = fd_anova(uOriented)
        if pint is not None:
            if pint < padj:
                fdirIntsAlley.append(cellid)
            
            totalAnovas += 1
            
            
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
    
    
    
