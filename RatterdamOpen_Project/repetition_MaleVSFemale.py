"""
Repetition Project
Analyzing male vs female rats for main results
will look at 
1) degree of repetition (operatioanally and OAS)
2) directionality
3) same arm diff arm analysis
"""

import pickle, numpy as np, os, sys
from matplotlib import pyplot as plt
from scipy.stats import binom_test
import repetition_manuscript_defaults as MDef
from matplotlib.ticker import MaxNLocator
import pandas as pd
from collections import Counter

plt.ion()

#Load dict containing all recording day datasets
# structure is rat > day > unit
# unit is a Unit() class 
with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   
    
alleydf = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv")

rat_sex_dict = {'R765':'M',
                'R781':'M',
                'R808':'F',
                'R859':'F',
                'R886':'M'
                }

male_rats = ['R765', 'R781', 'R886']
female_rats = ['R808', 'R859']


#%% Quantification of repetition using operational definition
from scipy.stats import chi2_contingency

rep_mf = {'M':[],
          'F':[]
          }

rat_reps = {rat:{day:[] for day in superpop[rat].keys()} for rat in superpop.keys()}

for rat in superpop.keys():
    ratsex = rat_sex_dict[rat]
    for day in superpop[rat].keys():
        for unit in superpop[rat][day]['units'].values():
            rep_mf[ratsex].append(unit.repeating)


#%% Directionality by Sex
# 
from scipy.stats import mannwhitneyu

pvals_mf = {'M':[],
            'F':[]
            }

meanDiff_mf = {'M':[],
               'F':[]
            }

for rat, rdf in alleydf.groupby('Rat'):
    ratsex = rat_sex_dict[rat]

    for orien,odf in rdf.groupby('Orientation'):
        for fname, fgroup in odf.groupby("FieldID"):
            
            dirs = np.unique(fgroup.CurrDir)
            dirA = fgroup[fgroup.CurrDir==dirs[0]]
            dirB = fgroup[fgroup.CurrDir==dirs[1]]
            
            try:
                pvals_mf[ratsex].append(mannwhitneyu(dirA.Rate, dirB.Rate).pvalue)
                meanDiff_mf[ratsex].append(abs(dirA.Rate.mean()-dirB.Rate.mean()))
            except:
                pass

for s in ['M', 'F']:
    print(s)
    pvals_mf[s] = np.asarray(pvals_mf[s])
    print(sum(pvals_mf[s]<0.05)/len(pvals_mf[s]))
    

chi2_contingency([
        [sum(pvals_mf['M']<0.05),sum(pvals_mf['F']<.05)],
        [len(pvals_mf['M']), len(pvals_mf['F'])]
                ])

#%% Same arm / different arm 
import itertools
from scipy.stats import linregress 


ratsex_target = 'F'

if ratsex_target == 'M':
    multiratdf = alleydf[alleydf.Rat.isin(male_rats)]
elif ratsex_target == 'F':
    multiratdf = alleydf[alleydf.Rat.isin(female_rats)]

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

for cellid, cell in multiratdf.groupby("CellID"):
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

olsAll = linregress(x,y)
olsSame = linregress(x[arm_shareness==1], y[arm_shareness==1])
olsDiff = linregress(x[arm_shareness==0], y[arm_shareness==0])

size=250
fig, ax = plt.subplots()
ax.set_aspect('equal')


#ax.scatter(x,y, c='black',label='All Data')
ax.scatter(x[arm_shareness==0], y[arm_shareness==0], color='cornflowerblue', edgecolor='navy', s=size)
ax.scatter(x[arm_shareness==1], y[arm_shareness==1], color='lightcoral', edgecolor='firebrick',s=size)

 
ax.text(2,-2,"$R^2_{all}$ ="+ f"{round(olsAll.rvalue**2,3)}",fontsize=45)
ax.text(2,-3,"$R^2_{same}$ =" + f"{round(olsSame.rvalue**2,3)}",fontsize=45)
ax.text(2,-4,"$R^2_{diff}$ =" + f"{round(olsDiff.rvalue**2,3)}",fontsize=45)

ax.set_ylabel("Signed Normalized\n Directionality Field A", fontsize=MDef.ylabelsize)
ax.set_xlabel("Signed Normalized\n Directionality Field B", fontsize=MDef.xlabelsize)
ax.set_title("Pairwise Directionality Comparison of Aligned Repeating Fields \nwithin each Repeating Neuron",
             fontsize=MDef.titlesize)

xlim, ylim = ax.get_xlim(), ax.get_ylim()
xfit = np.linspace(xlim[0], xlim[1],100)

for fit, linecolor, errcolor,label in zip([olsAll, olsSame, olsDiff], 
                                    ['black', 'red', 'blue'], 
                                    ['grey', 'lightcoral', 'cornflowerblue'],
                                    ['All Data', 'Same Segment', 'Different Segment']):

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

lgnd = ax.legend(prop={'size':32})
ax.set_title(ratsex_target)
# %%
