# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 13:15:33 2022

@author: whockei1

Followup analyses for the same arm / different arm analysis

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
from scipy.stats import linregress 
from matplotlib import patches
import pickle 
from scipy.interpolate import PchipInterpolator as pchip
from collections import Counter 
import scipy.stats


alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   


#%% 

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

ii = 0

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
                    ii += 1
                    if np.unique(field.Alleys).shape[0] > 1:
                        print(oname, fid)
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


#%% Correlate difference in sampling bias with product of directional tuning
# this is hard to interpret - keeping here in case i want to do something with it or spend more time looking at it

prods = x*y
allBiasDiffs = []
for pair in pair_info:
    
    alleyA, alleyB = pair['alleyA'], pair['alleyB']
    turns = superpop[pair['rat']][pair['day']]['refturns']
    
    if pair['orien'] == 'V':
        dirs = ['1', '3']
    elif pair['orien'] == 'H':
        dirs = ['2', '4']
    
    biases = []
    for alley in [alleyA, alleyB]: 
    
        alleyDirSampling = turns[turns['Alley+']==alley]['Allo+']
        c = Counter(alleyDirSampling)
        biases.append(np.diff(np.asarray([c[i] for i in dirs])/sum(c.values())))
        
    biasDiff = biases[0] - biases[1]
    allBiasDiffs.append(biasDiff)
    
    
#%% Make two graphs, one for same arm pairs, the other for different arm pairs
# And within each plot, plot the directional sampling bias for field A against B

from matplotlib import cm
from matplotlib import colors as mpl_colors


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
                linewidth=2)
    
    if i == 1:

        fig.colorbar(sc, ax=cax)
    
    cax.set_xlabel("Field A Sampling Bias",fontsize=25)
    cax.set_ylabel("Field B Sampling Bias",fontsize=25)
    ymin, ymax = cax.get_ylim()
    xmin, xmax = cax.get_xlim()
    cax.hlines(0,xmin, xmax,color='k',linestyle='--')
    cax.vlines(0,ymin, ymax,color='k',linestyle='--')
    cax.set_title(label,fontsize=30)
    cax.spines['top'].set_visible(False)
    cax.spines['right'].set_visible(False)
    cax.spines['left'].set_linewidth(3)
    cax.spines['bottom'].set_linewidth(3)
    cax.tick_params(axis='both', which='major', labelsize=20)

    
plt.suptitle("Directional Sampling Bias For Repeating Field Pairs",fontsize=30)
            
            
            
#%% Distance between fields

fig, ax = plt.subplots()
cax = fig.axes[0]
for shareness, label,color in zip([0,1], ["Different Arm", "Same Arm"],['red','blue']):
    
    dists = []
    
    for pair in pair_info[arm_shareness==shareness]:
        
        fidA, fidB = pair['fidA'], pair['fidB']
        orien = pair['orien']
        rat, day = pair['rat'], pair['day']
       
        turns = superpop[rat][day]['refturns']
        
        #cell names are the same, I do it twice because i forgot, it' fine 
        
        cellnameA = np.unique(alleydf[(alleydf['FieldID']==fidA)
                                      &(alleydf['Orientation']==orien)]['CellName'])[0]
        
        fieldA = np.unique(alleydf[(alleydf['FieldID']==fidA)
                                      &(alleydf['Orientation']==orien)]['FieldNum'])[0]
        
        
        cellnameB = np.unique(alleydf[(alleydf['FieldID']==fidB)
                                      &(alleydf['Orientation']==orien)]['CellName'])[0]
        
        fieldB = np.unique(alleydf[(alleydf['FieldID']==fidB)
                                      &(alleydf['Orientation']==orien)]['FieldNum'])[0]

        
        comA = superpop[rat][day]['units'][cellnameA].repUnit.PF[fieldA].com
        comB = superpop[rat][day]['units'][cellnameB].repUnit.PF[fieldB].com
        
        dist = np.linalg.norm(np.asarray(comA)-np.asarray(comB))

        dists.append(dist)
            
    
    cax.hist(dists,color=color,label=label,density=True,alpha=0.5)
    
cax.set_xlabel("Euclidian Distance between COMs")
cax.set_ylabel("Normalized Frequency")
    
plt.title(" Distance Between Fields of Repeating Field Pairs")
plt.legend()


#%% Correlate sampling over time 


nptsinterp = 200
fig, ax = plt.subplots()
cax = fig.axes[0]
allcorrs = []
prods = []
sharenesses = []
for shareness, label,color in zip([0,1], ["Different Arm", "Same Arm"],['red','blue']):
    
    corrs = []
    
    for pair in pair_info[arm_shareness==shareness]:
        
        alleyA, alleyB = pair['alleyA'], pair['alleyB']
        orien = pair['orien']
        rat, day = pair['rat'], pair['day']
               
        turns = superpop[rat][day]['refturns']
        
        interp_ts = []
        
        for alley in [alleyA, alleyB]:
            
            turntimes = turns[turns['Alley+']==str(alley)]['Ts exit'].astype(float)
            turndirs = turns[turns['Alley+']==str(alley)]['Allo+'].values.astype(int)
            
            if pair['orien'] == 'V':
                turndirs[turndirs==1] = 1 # N, E are 1. S,W are 0 
                turndirs[turndirs==3] = 0
            elif pair['orien'] == 'H':
                turndirs[turndirs==2] = 1
                turndirs[turndirs==4] = 0
        
            alleyDirSampling = np.column_stack((turntimes, turndirs))        

            ts = np.linspace(alleyDirSampling[0,0], alleyDirSampling[-1,0], num=nptsinterp)
            pcf = pchip(alleyDirSampling[:,0], alleyDirSampling[:,1])
            pcf_interp = pcf(ts)
            
            # if rat == 'R781' and day == 'D3':
            #     line = plt.plot(ts,pcf_interp)
            #     plt.plot(alleyDirSampling[:,0], alleyDirSampling[:,1], color=line[0].get_color())
            
            interp_ts.append(np.column_stack((ts, pcf_interp)))
        
        corr = scipy.stats.pearsonr(interp_ts[0][:,1], interp_ts[1][:,1])[0]
        
        corrs.append(corr)
        allcorrs.append(corr)
        prods.append(pair['prod'])
        sharenesses.append(shareness)
            
    
    cax.hist(corrs,color=color,label=label,density=True,alpha=0.5)
    
cax.set_xlabel("Pearson Correlation between PCHIP-interpolated time series")
cax.set_ylabel("Normalized Frequency")
    
plt.title("Correlation between Sampling of Each Field in Pairs")
plt.legend()


#%% Plot scatter of directionality product against the timeseries correlation of directional sampling
sharenesses = np.asarray(sharenesses)
allcorrs = np.asarray(allcorrs)
prods = np.asarray(prods)

fig, ax = plt.subplots()

colors = ['lightcoral' if i == 1 else 'cornflowerblue' for i in sharenesses]

for i,label,c in zip([0,1],['Different Arm', 'Same Arm'],['cornflowerblue', 'lightcoral']):
    ax.scatter(allcorrs[sharenesses==i], prods[sharenesses==i], 
               c=c,
               edgecolor='black',
               linewidth=2,
               s=200,
               label = label
               )

ax.set_xlabel("Pearson Correlation between PCHIP-interpolated time series",fontsize=25)
ax.set_ylabel("Signed Product of Directional Tuning of Field Pair",fontsize=25)
ax.set_title("Directional Sampling Correlation Versus Directionality Product",fontsize=30)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.tick_params(axis='both', which='major', labelsize=20)
lgnd = plt.legend(prop={'size':25})
lgnd.get_frame().set_linewidth(2)





