"""
Repetition Manuscript Supplementary Figure Code 
WH 2022-05-10
"""
#%%
import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from yaml import MappingNode
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import ratterdam_Defaults as Def
from scipy.stats import sem
import repetition_manuscript_defaults as MDef
import pickle
import itertools

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)   


#%% Process data 
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
samearm_alleys, diffarm_alleys = [], []

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

                if armA == armB:
                    samearm_alleys.append([rat, day, alleyA, alleyB])
                else:
                    diffarm_alleys.append([rat, day, alleyA, alleyB])

                arm_shareness.append(int(armA==armB))

arm_shareness = np.asarray(arm_shareness)
pair_info = np.asarray(pair_info)
x = np.asarray(x)
y = np.asarray(y)           

#%% Fig S4 - Supplementary . Shuffled GLM Emmeans analysis of retrospective/prospective/current direction
df = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Supplementary_Figures\\shuffle_nonlocalDirection\\2022-05-10_shuffle_emmeans.csv")

bins = np.linspace(0,0.3,20)

fig, ax = plt.subplots(1,3,figsize=(10,8))
#proportions are hardcoded from data, see repetition_manuscript_statistics
# and they are current as of 5/10/2022
for i,(dirtype, dirlabel, realprop, cax) in enumerate(zip(["total_current_responsive",
                              "total_previous_responsive",
                              "total_next_responsive"],
                              [
                               "Current Direction",
                               "Retrospective Direction",
                               "Prospective Direction"
                              ],
                              [32/127,
                              27/160,
                              20/151],
                              fig.axes
                                )):
    cax.hist(df[dirtype],bins=bins,facecolor='grey',edgecolor='black',linewidth=2)
    cax.vlines(np.percentile(df[dirtype], 95), 0, cax.get_ylim()[1],color='k')
    cax.vlines(realprop,0, cax.get_ylim()[1],color='r')
    cax.set_title(dirlabel, fontsize=25)
    if i == 0:
        cax.set_ylabel("Frequency", fontsize=25)
    if i == 1:
        cax.set_xlabel("Proportion of Cells Responding to Turn Type", fontsize=25)
    cax.spines['right'].set_visible(False)
    cax.spines['top'].set_visible(False)
    cax.tick_params(axis='both', which='major', labelsize=22)


#%% S5A-Field distances
fig, ax = plt.subplots()
cax = fig.axes[0]
all_dists = {}
for shareness, label,fcolor,ecolor,a in zip([0,1], 
                                    ["Different Corridor", "Same Corridor"],
                                    ['blue', 'red'],
                                    ['navy','firebrick'],
                                    [1,0.6]):
    all_dists[shareness] = []
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
        all_dists[shareness].append(dist)
            
    
    cax.hist(dists,
                color=fcolor,
                edgecolor=ecolor,
                linewidth=2,
                label=label,
                density=True,
                alpha=a
    )
    
cax.set_xlabel("Euclidian Distance between COMs (bins)", fontsize=25)
cax.set_ylabel("Normalized Frequency", fontsize=25)
cax.spines['right'].set_visible(False)
cax.spines['top'].set_visible(False)
cax.tick_params(axis='both', which='major', labelsize=22)
    
plt.title(" Distance Between Fields of Repeating Field Pairs", fontsize=25)
plt.legend(fontsize=20)

# %% S5B - Path length

## Analysis changed 5/31/22 see stats file. basically now we are looking at average path length
# from unique alley pairs. not all paths through alley pairs and do that each time the alley pair 
# was encountered in data. purpose is making samples independent 

samearm_interPairLength = []
diffarm_interPairLength = []

unique_samearm_alleys = [list(y) for y in set([tuple(x) for x in samearm_alleys])]
unique_diffarm_alleys = [list(y) for y in set([tuple(x) for x in diffarm_alleys])]

for entry in unique_samearm_alleys:
    rat, day, alleyA, alleyB = entry
    
    pair_lengthToParner = []

    refturns = superpop[rat][day]['refturns']
    
    # logic is we are using turn df that is unfiltered so indices are each visit # in order
    # so for each visit to Alley A, find min diff of that index to all the indices of alley B
    # i.e. the number of visits (to any alley) between them
    for visit in refturns[refturns['Alley+']==str(alleyA)].index:
       lengthToPartner = np.min(abs(visit-refturns[refturns['Alley+']==str(alleyB)].index))
       pair_lengthToParner.append(lengthToPartner)
    
    samearm_interPairLength.append(np.nanmean(pair_lengthToParner))
       
       
for entry in unique_diffarm_alleys:
    rat, day, alleyA, alleyB = entry

    pair_lengthToParner = []
    
    refturns = superpop[rat][day]['refturns']
    
    # logic is we are using turn df that is unfiltered so indices are each visit # in order
    # so for each visit to Alley A, find min diff of that index to all the indices of alley B
    # i.e. the number of visits (to any alley) between them
    for visit in refturns[refturns['Alley+']==str(alleyA)].index:
       lengthToPartner = np.min(abs(visit-refturns[refturns['Alley+']==str(alleyB)].index))
       pair_lengthToParner.append(lengthToPartner-1)
    diffarm_interPairLength.append(np.nanmean(pair_lengthToParner))
#%%
fig, ax = plt.subplots()

mymax = max(max(samearm_interPairLength), max(diffarm_interPairLength))

ax.hist(samearm_interPairLength,
        density=False, 
        bins=np.linspace(0,mymax,25),
        color='red',
        edgecolor='firebrick',
        linewidth=2,
        label='Same Corridor'
        )

ax.hist(diffarm_interPairLength,
        density=False,
        stacked=False,
        bins=np.linspace(0,mymax,25),
        color='blue',
        edgecolor='navy',
        linewidth=2,
        alpha=0.6,
        label='Different Corridor'
        )


ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_xlim([0,40])
ax.set_ylabel("Normalized, Stacked Frequency", fontsize=MDef.ylabelsize)
ax.set_xlabel("Distance in alley visits between visits to alleys in pair", fontsize=MDef.xlabelsize)
ax.set_title("Relationship between arms alleys are on and the path lengths between alleys", fontsize=MDef.titlesize)
plt.legend(fontsize=20)
# %%
from scipy.stats import mannwhitneyu
print(mannwhitneyu(samearm_interPairLength, diffarm_interPairLength))


