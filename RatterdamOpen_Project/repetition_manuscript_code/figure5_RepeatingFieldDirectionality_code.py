# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:33:35 2022

@author: whockei1


Script to generate the panels for Figure 5 of repetition manuscript
Each panel in the paper corresponds to a code cell here. Code cells in python
begin with "#%%". They can be run independently. The whole script can be run to generate all panels.

Figure 5 analyzes the distribution of directional tuning across repeating fields.
Pairs of fields are further divided into whether they are on the same corridor or not

A - two example repeating neurons. Ratemaps and bar graphs of each field's directional tuning.
B - Scatterplot showing comparison of directionality index for each pair of repeating fields
C - B, but for pairs on the same corridor
D - B, but for pairs on different corridor
E - shuffling test showing significance of difference in correlation
F - 

"""

import numpy as np, matplotlib.pyplot as plt, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import itertools
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import repetition_manuscript_defaults as MDef
from scipy.stats import linregress 
import copy
from collections import Counter 

from matplotlib import cm
from matplotlib import colors as mpl_colors
import pickle 


plt.ion()
alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  


#%% Panel 5A Example units. Plot ratemap and for each field plot bars of average firing
# for each direction +/- SEM


cmap = util.makeCustomColormap()
example_cell_list = [["R781","D4","TT3\\cl-maze1.2","V"],
                     ["R859","D1","TT1\\cl-maze1.1","V"]
                     
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
    
                        

#%%  Cell to calc pairwise field polarity 

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
arm_shareness = []
pair_info = []


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
                
                if np.unique(field.Alleys).shape[0] == 1:
           
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
                fidA, fidB = fids[i], fids[j]
                diffA, diffB = signed_diffs[i], signed_diffs[j]

                x.append(diffA)
                y.append(diffB)
                
                pair_info.append({'rat':rat, 
                                  'day':day, 
                                  'alleyA':alleyA, 
                                  'alleyB':alleyB, 
                                  'fidA':fidA, 
                                  'fidB':fidB, 
                                  'orien':oname,
                                  'prod':diffA*diffB})
                
                # this is clumsy but getting each arm each field is on
                # assumes one and only one hit when checking orientation_structure 
                for arm, armAlleys in orientation_structure.items():
                    if alleyA in armAlleys:
                        armA = arm        
                for arm, armAlleys in orientation_structure.items():
                    if alleyB in armAlleys:
                        armB = arm
                arm_shareness.append(int(armA==armB))
                    
arm_shareness = np.asarray(arm_shareness)
pair_info = np.asarray(pair_info)
x = np.asarray(x)
y = np.asarray(y)                    

olsAll = linregress(x,y)
olsSame_shuff = linregress(x[arm_shareness==1], y[arm_shareness==1])
olsDiff_shuff = linregress(x[arm_shareness==0], y[arm_shareness==0])

#%% Panel B 

size=250
fig, ax = plt.subplots()
ax.set_aspect('equal')


ax.scatter(x,y, 
            facecolor='grey',
            edgecolor='black',
            label='All Data',
            linewidth=2)

ax.text(2,-2,"$R^2_{all}$ ="+ f"{round(olsAll.rvalue**2,3)}",fontsize=45)


ax.set_ylabel("Signed Normalized\n Directionality Field A", fontsize=MDef.ylabelsize)
ax.set_xlabel("Signed Normalized\n Directionality Field B", fontsize=MDef.xlabelsize)
ax.set_title("Pairwise Directionality Comparison of Aligned Repeating Fields \nwithin each Repeating Neuron",
             fontsize=MDef.titlesize)

xlim, ylim = ax.get_xlim(), ax.get_ylim()
xfit = np.linspace(xlim[0], xlim[1],100)

for fit, linecolor, errcolor,label in zip([olsAll], ['black'], ['grey'], ['All Data']):

    yfit = [(fit.slope*i)+fit.intercept for i in xfit]
    ax.plot(xfit,yfit,color=linecolor,linewidth=3)
    ax.fill_between(xfit,yfit-fit.stderr,yfit+fit.stderr,color=errcolor,alpha=0.7,label=label)


ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.hlines(0, xlim[0], xlim[1], linestyle='--',color='k',linewidth=2)
ax.vlines(0, ylim[0], ylim[1], linestyle='--',color='k',linewidth=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylim([-2.5,3])
ax.set_xlim([-3.5,4])


#%% Panels C, D
size=250
tick_spacing = 0.5

import matplotlib.ticker as ticker 

for shareness, colors, label, fit in zip([0,1], 
                                        [['cornflowerblue', 'navy'], ['lightcoral', 'firebrick']],
                                        ['Different Corridor', 'Same Corridor'],
                                        [olsDiff_shuff, olsSame_shuff]):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    ax.scatter(x[arm_shareness==shareness], y[arm_shareness==shareness], color=colors[0], edgecolor=colors[1], s=size)
    
    ax.set_ylabel("Signed Normalized\n Directionality Field A", fontsize=MDef.ylabelsize)
    ax.set_xlabel("Signed Normalized\n Directionality Field B", fontsize=MDef.xlabelsize)
    
    ax.set_ylim([-2.5,3])
    ax.set_xlim([-3.5,4])
    
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xfit = np.linspace(xlim[0], xlim[1],100)

    
    yfit = [(fit.slope*i)+fit.intercept for i in xfit]
    ax.plot(xfit,yfit,color=colors[1],linewidth=3)
    ax.fill_between(xfit,yfit-fit.stderr,yfit+fit.stderr,color=colors[0],alpha=0.7,label=label)
    
    
    ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
    ax.hlines(0, xlim[0], xlim[1], linestyle='--',color='k',linewidth=2)
    ax.vlines(0, ylim[0], ylim[1], linestyle='--',color='k',linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(MDef.spine_width)
    ax.spines['bottom'].set_linewidth(MDef.spine_width)
    ax.set_title(label,fontsize=40)
    
    xticks = np.arange(xlim[0], xlim[1], 0.5)
    yticks = np.arange(ylim[0], ylim[1], 0.5)


#%% Panels E, F

arm_shareness_copy = copy.deepcopy(arm_shareness)

real_olsfit_same = linregress(x[arm_shareness==1], y[arm_shareness==1])
real_olsfit_diff = linregress(x[arm_shareness==0], y[arm_shareness==0])
real_yfit_same = [(real_olsfit_same.slope*i)+real_olsfit_same.intercept for i in np.linspace(min(x), max(x),100)]
real_yfit_diff = [(real_olsfit_diff.slope*i)+real_olsfit_diff.intercept for i in np.linspace(min(x), max(x),100)]

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
    
    
    olsSame_shuff = linregress(same_shuff_x, same_shuff_y)
    yfit_same = [(olsSame_shuff.slope*i)+olsSame_shuff.intercept for i in np.linspace(min(x), max(x),100)]
    
    olsDiff_shuff = linregress(diff_shuff_x, diff_shuff_y)
    yfit_diff = [(olsDiff_shuff.slope*i)+olsDiff_shuff.intercept for i in np.linspace(min(x), max(x),100)]
    
    
    shuff_r2_same.append(olsSame_shuff.rvalue**2)
    shuff_lines_same.append(yfit_same)
    shuff_r2_diff.append(olsDiff_shuff.rvalue**2)
    shuff_lines_diff.append(yfit_diff)
    
    shuff_r2_differences.append(olsSame_shuff.rvalue**2 - olsDiff_shuff.rvalue**2)

shuff_lines_same = np.asarray(shuff_lines_same)
shuff_lines_diff = np.asarray(shuff_lines_diff)    
    
# shuffle regression lines 
fig, ax = plt.subplots()

ax.plot(xfit, real_yfit_same, 
        color='red', 
        linewidth=3,
        label = 'Same segment'
        )
ax.plot(xfit, real_yfit_diff, 
        color='navy', 
        linewidth=3,
        label = 'Different segment'
        )
ax.plot(xfit, yfit, color='black',linewidth=3,zorder=99, label = 'All data')


diff_lower = np.percentile(shuff_lines_diff, 2.5,axis=0)
diff_upper = np.percentile(shuff_lines_diff, 97.5,axis=0)
same_lower = np.percentile(shuff_lines_same, 2.5,axis=0)
same_upper = np.percentile(shuff_lines_same, 97.5,axis=0)

ax.fill_between(xfit, diff_lower, diff_upper,
                color='cornflowerblue',
                alpha=0.6,
                label="Different arm $2.5^{th}$ / $97.5^{th}$ confidence band")
ax.plot(xfit,diff_lower,
        color='navy',
        linewidth=3)
ax.plot(xfit,diff_upper,
        color='navy',
        linewidth=3
        )

ax.fill_between(xfit, same_lower, same_upper,
                color='lightcoral',
                alpha=0.6,
                label = 'Same arm $2.5^{th}$ / $97.5^{th}$ confidence band')
ax.plot(xfit,same_lower,
        color='firebrick',
        linewidth=3)
ax.plot(xfit,same_upper,
        color='navy', 
        linewidth=3
        )

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

#%% Panel G 
allcolors, allx, ally = [], [], []
for shareness in [0,1]:
    
    sharex, sharey = [], []
    colors = []
    
    for pair in pair_info[arm_shareness==shareness]:
        
        colors.append(pair['prod'])
        
        alleyA, alleyB = pair['alleyA'], pair['alleyB']
       
        turns = superpop[pair['rat']][pair['day']]['refturns']
        
        # important to keep the ordering consistent: North, East are the positive
        # South, West are the negative. IOW, it's alphabetical. And the numeric codes
        # from Def.allocodedict (the coding scheme in the turns df) is consistent
        # with this N=1, S=3, 2=E, 4=W
        if pair['orien'] == 'V':
            dirs = ['1', '3']
        elif pair['orien'] == 'H':
            dirs = ['2', '4']
        
        biases = []
        for alley, l in zip([alleyA, alleyB], [sharex, sharey]): 
        
            alleyDirSampling = turns[turns['Alley+']==str(alley)]['Allo+']
            c = Counter(alleyDirSampling)
            l.append(np.diff(np.asarray([c[i] for i in dirs])/sum(c.values())))
            
    allx.append(sharex)
    ally.append(sharey)
    allcolors.append(colors)
    
    
    
fig, ax = plt.subplots(1,2)


mymin = min([min(i) for i in allcolors])
mymax = max([max(i) for i in allcolors])

for i, (cax,label,xs,ys,col) in enumerate(zip(fig.axes,
                                                 ['Different Arm', 'Same Arm'],
                                                 allx,
                                                 ally,
                                                 allcolors
                                             )    
                                             ):

    
    divnorm = mpl_colors.TwoSlopeNorm(vmin=mymin,vmax=mymax,vcenter=0)
    sc = cax.scatter(xs, ys,
                c=np.asarray(col),
                cmap = 'seismic',
                norm=divnorm,
                s=200,
                edgecolor='k',
                linewidth=0.5)
    
    if i == 1:

        fig.colorbar(sc, ax=cax)
    
    cax.set_xlabel("Field A Sampling Bias",fontsize=35)
    cax.set_ylabel("Field B Sampling Bias",fontsize=35)
    # cax.set_ylim([-1,1])
    # cax.set_xlim([-1,1])
    ymin, ymax = cax.get_ylim()
    xmin, xmax = cax.get_xlim()
    cax.hlines(0,xmin, xmax,color='k',linestyle='--')
    cax.vlines(0,ymin, ymax,color='k',linestyle='--')
    cax.set_title(label,fontsize=30)
    cax.spines['top'].set_visible(False)
    cax.spines['right'].set_visible(False)
    cax.spines['left'].set_linewidth(1.5)
    cax.spines['bottom'].set_linewidth(1.5)
    cax.tick_params(axis='both', which='major', labelsize=35, length=8)
    cax.set_aspect('equal')
    

    
plt.suptitle("Directional Sampling Bias For Repeating Field Pairs",fontsize=30)
            
plt.show()