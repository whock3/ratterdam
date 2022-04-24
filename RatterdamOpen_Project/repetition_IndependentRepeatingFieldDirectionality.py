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
import pickle 


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
arm_shareness = []


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
                if armA == armB:
                    diffA, diffB = signed_diffs[i], signed_diffs[j]
                    x.append(diffA)
                    y.append(diffB)
                    arm_shareness.append(1)
                elif armA != armB:
                    diffA, diffB = signed_diffs[i], signed_diffs[j]
                    x.append(diffA)
                    y.append(diffB)
                    arm_shareness.append(0)
                    
arm_shareness = np.asarray(arm_shareness)
x = np.asarray(x)
y = np.asarray(y)                    
                   
            
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
    
    
#%% Looking at relationship between whether pair of fields is on same arm/different arm
# and the path length / path stereotype between those locations.

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
samearm_alleys = []
diffarm_alleys = []

for cellid, cell in alleydf.groupby("CellID"):
    
    rat = np.unique(cell.Rat)[0]
    day = np.unique(cell.Day)[0]
    
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
                  
                if armA == armB:
                    samearm_alleys.append([rat, day, alleyA, alleyB])
                else:
                    diffarm_alleys.append([rat, day, alleyA, alleyB])
                    
#%% go through alley pairs and check avg path length and dist (d-l dist)

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  

samearm_interPairLength = []
diffarm_interPairLength = []

for entry in samearm_alleys:
    rat, day, alleyA, alleyB = entry
    
    refturns = superpop[rat][day]['refturns']
    
    # logic is we are using turn df that is unfiltered so indices are each visit # in order
    # so for each visit to Alley A, find min diff of that index to all the indices of alley B
    # i.e. the number of visits (to any alley) between them
    for visit in refturns[refturns['Alley+']==str(alleyA)].index:
       lengthToPartner = np.min(abs(visit-refturns[refturns['Alley+']==str(alleyB)].index))
       samearm_interPairLength.append(lengthToPartner)
       
       
for entry in diffarm_alleys:
    rat, day, alleyA, alleyB = entry
    
    refturns = superpop[rat][day]['refturns']
    
    # logic is we are using turn df that is unfiltered so indices are each visit # in order
    # so for each visit to Alley A, find min diff of that index to all the indices of alley B
    # i.e. the number of visits (to any alley) between them
    for visit in refturns[refturns['Alley+']==str(alleyA)].index:
       lengthToPartner = np.min(abs(visit-refturns[refturns['Alley+']==str(alleyB)].index))
       diffarm_interPairLength.append(lengthToPartner)
       
    
#%% plot 

fig, ax = plt.subplots()

mymax = max(max(samearm_interPairLength), max(diffarm_interPairLength))

ax.hist(samearm_interPairLength,
        density=True, 
        bins=np.linspace(0,mymax,100),
        color='red',
        label='Field pair on same arm'
        )

ax.hist(diffarm_interPairLength,
        density=True,
        stacked=True,
        bins=np.linspace(0,mymax,100),
        color='k',
        alpha=0.6,
        label='Field pair on different arms'
        )


ax.tick_params(axis='both', which='major', labelsize=22)

ax.set_ylabel("Normalized, Stacked Frequency", fontsize=25)

ax.set_xlabel("Distance in alley visits between visits to alleys in pair", fontsize=25)

ax.set_title("Relationship between arms alleys are on and the path lengths between alleys", fontsize=30)


#%% Same arm vs different arm, shuffling to see R2 under null hypothesis
# that they contribute equally to overall R2 (of about 4.3%)

import copy

xfit = np.linspace(min(x), max(x),100)

real_olsfit = linregress(x,y)
yfit = [(real_olsfit.slope*i)+real_olsfit.intercept for i in np.linspace(min(x), max(x),100)]
real_olsfit_same = linregress(x[arm_shareness==1], y[arm_shareness==1])
real_olsfit_diff = linregress(x[arm_shareness==0], y[arm_shareness==0])
real_yfit_same = [(real_olsfit_same.slope*i)+real_olsfit_same.intercept for i in np.linspace(min(x), max(x),100)]
real_yfit_diff = [(real_olsfit_diff.slope*i)+real_olsfit_diff.intercept for i in np.linspace(min(x), max(x),100)]

arm_shareness_copy = copy.deepcopy(arm_shareness)

nshuff = 1000

shuff_r2_same = []
shuff_lines_same = []
shuff_r2_diff = []
shuff_lines_diff = []
shuff_r2_differences = [] # diff above means diff arms, here differences means the arthimetic difference: r2_same - r2_diff
nsame = len(np.where(arm_shareness==1)[0])
ndiff = len(np.where(arm_shareness==0)[0])

# same always < diff, not by definition but bc its harder to be on same arm. but code it flexibly 

ndownsample = min(nsame, ndiff)

for s in range(nshuff):
    
    if s%10 == 0:
        print(s)
    
    arm_shareness_copy = np.random.permutation(arm_shareness_copy)
    
    same_shuff_x = x[np.random.choice(np.where(arm_shareness_copy==1)[0], ndownsample, replace=False)]
    same_shuff_y = y[np.random.choice(np.where(arm_shareness_copy==1)[0], ndownsample, replace=False)]

    diff_shuff_x = x[np.random.choice(np.where(arm_shareness_copy==0)[0], ndownsample, replace=False)] 
    diff_shuff_y = y[np.random.choice(np.where(arm_shareness_copy==0)[0], ndownsample, replace=False)]
    
    
    olsfit_same = linregress(same_shuff_x, same_shuff_y)
    yfit_same = [(olsfit_same.slope*i)+olsfit_same.intercept for i in np.linspace(min(x), max(x),100)]
    
    olsfit_diff = linregress(diff_shuff_x, diff_shuff_y)
    yfit_diff = [(olsfit_diff.slope*i)+olsfit_diff.intercept for i in np.linspace(min(x), max(x),100)]
    
    
    shuff_r2_same.append(olsfit_same.rvalue**2)
    shuff_lines_same.append(yfit_same)
    shuff_r2_diff.append(olsfit_diff.rvalue**2)
    shuff_lines_diff.append(yfit_diff)
    
    shuff_r2_differences.append(olsfit_same.rvalue**2 - olsfit_diff.rvalue**2)

shuff_lines_same = np.asarray(shuff_lines_same)
shuff_lines_diff = np.asarray(shuff_lines_diff)    
    
# shuffle regression lines 
fig, ax = plt.subplots()

ax.plot(xfit, real_yfit_same, color='green', label = 'Same arm')
ax.plot(xfit, real_yfit_diff, color='orange', label = 'Different arm')
ax.plot(xfit, yfit, color='black',linewidth=3,zorder=99, label = 'All data')

ax.plot(xfit,np.percentile(shuff_lines_diff,97.5,axis=0),
                            color='orange',
                            linewidth=2,
                            linestyle='--', 
                            label = "Different arm $2.5^{th}$ / $97.5^{th}$ confidence band")

ax.plot(xfit,np.percentile(shuff_lines_same,97.5,axis=0),
                            color='green',
                            linewidth=2,
                            linestyle='--',
                            label = 'Same arm $2.5^{th}$ / $97.5^{th}$ confidence band')

ax.plot(xfit,np.percentile(shuff_lines_diff,2.5,axis=0),
                            color='orange',
                            linewidth=2,
                            linestyle='--')

ax.plot(xfit,np.percentile(shuff_lines_same,2.5,axis=0),
                            color='green',
                            linewidth=2,
                            linestyle='--')

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Normalized Directional Rate (Hz)",fontsize=Def.ylabelsize)
ax.set_xlabel("Normalized Directional Rate (Hz)",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':44})
#change the marker size manually for both lines
lgnd.get_frame().set_linewidth(2)

# shuffled difference in r2 
fig, ax = plt.subplots()
realr2diff = real_olsfit_same.rvalue**2 - real_olsfit_diff.rvalue**2
shuffpercentile = np.percentile(shuff_r2_differences, 95)
plt.hist(shuff_r2_differences,bins=25, label="Shuffle")
plt.vlines(shuffpercentile,0,400,color='k', label="$95^{th}$ percentile")
plt.vlines(realr2diff,0,400,color='r', label= 'Observed')

ax.tick_params(axis='both', which='major', labelsize=Def.ticksize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.set_ylabel("Frequency",fontsize=Def.ylabelsize)
ax.set_xlabel("Difference in $R^2$",fontsize=Def.xlabelsize)
lgnd = plt.legend(prop={'size':44})
#change the marker size manually for both lines
lgnd.get_frame().set_linewidth(2)
    
    
    
#%% Looking at polarity scatter by PD, ND

opposite_dir_lookup = {'N':'S',
                       'E':'W',
                       'W':'E',
                       'S':'N'
                        }

turn_type = field.NextDir

for cellid, cell in alleydf.groupby("CellID"):
    
    rat = np.unique(cell.Rat)[0]
    day = np.unique(cell.Day)[0]
    
    for oname, ocell in cell.groupby("Orientation"):
        
        orientation_structure = maze_arm_organization[oname]
        
        if np.unique(ocell.FieldID).shape[0] > 1:
        
            signed_diffs = []
            alleys = []
            fids = []

            for fid, field in ocell.groupby("FieldID"):
            
            # 2022-04-20 come back to this
            # for one, it doesnt distinguish between which side of the alley
            # youre on and that could make a difference 
            
                dirs = sorted(np.unique(turn_type)) # sort alphabetically to keep signs consistent below 
                if len(dirs)>1:
                    meanDirDiffs = {}
                    for direction in dirs:
                        meanDirDiffs[direction] = field[turn_type==direction].Rate.mean()
                    maxdir = max(meanDirDiffs, key = lambda x: meanDirDiffs[x])
                    maxdirrate = meanDirDiffs[maxdir]
                    
                    #try to get opposite dir FR, if that doesnt exist just pick randomly
                    try:
                        mindirrate = meanDirDiffs[opposite_dir_lookup[maxdir]]
                    except: 
                        # this convoluted line randomly selects a direction to be 
                        # compared with the max FR direction to get a directionality value
                        # the reshaping business is bc choice needs 1d and argwhere returns 2d
                        mindir = dirs[np.random.choice(np.argwhere(dirs!=maxdir).reshape(np.argwhere(dirs!=maxdir).shape[0],))]
                        mindirrate = meanDirDiffs[mindir]
                        

                diff = maxdirrate - mindirrate
                interdiffs.append(diff)

                
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
                  
                if armA == armB:
                    samearm_alleys.append([rat, day, alleyA, alleyB])
                else:
                    diffarm_alleys.append([rat, day, alleyA, alleyB])