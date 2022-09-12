#%%

import numpy as np, matplotlib.pyplot as plt, pandas as pd
from scipy.stats.stats import ttest_1samp
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
import scipy.stats

from matplotlib import cm
from matplotlib import colors as mpl_colors
import pickle 
import repetition_manuscript_defaults as MDef 
import newAlleyBounds as nab 
from matplotlib import path


plt.ion()
alleydatapath = "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
alleydf = pd.read_csv(alleydatapath)

if 'Code' not in alleydf.columns:
    codes = []
    for r, row in alleydf.iterrows():
       code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
       codes.append(code)
    alleydf = alleydf.assign(Code=codes)

with open("E:\\Ratterdam\\R_data_repetition\\20220405-124315_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)  

# %% Correlate cluster quality with number of fields
# here this is being done for included clusters, i.e. quality >=3
all_quals, all_nfields = [], []
for rat in list(superpop.keys()):
    for day in list(superpop[rat].keys()):
        dataset = superpop[rat][day]
        df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
        clustList, clustQuals = util.getClustList(df)

        print(max(clustQuals))

        for unitName, unit in dataset['units'].items():
            quality = clustQuals[clustList.index(unitName)]
            nfields = len(unit.fields)
            all_quals.append(quality)
            all_nfields.append(nfields)
# %%
all_quals = np.asarray(all_quals)
all_nfields = np.asarray(all_nfields)

qnfs = []
for q in np.unique(all_quals):
    qnf = all_nfields[all_quals==q]
    qnfs.append(qnf)
plt.violinplot(qnfs)
# %% do for all cells regardless of quality
all_quals, all_nfields = [], []
Def.includeAllDetectedFields == True

for rat in list(superpop.keys()):
    for day in list(superpop[rat].keys()):
        dataset = superpop[rat][day]
        df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
        clustList, clustQuals = util.getClustList(df)

        for i, unitName in enumerate(clustList):
            print(rat,day,unitName)
            unit = RepCore.loadRepeatingUnit(rat,day,unitName)
            all_quals.append(clustQuals[i])
            all_nfields.append(len(unit.fields))
#%%
all_quals = np.asarray(all_quals)
all_nfields = np.asarray(all_nfields)
qnfs = []
for q in np.unique(all_quals):
    qnf = all_nfields[all_quals==q]
    qnfs.append(qnf)
plt.violinplot(qnfs)

# %% based on JK suggestion, test avg difference in quality rep vs nonrep

all_quals, rep = [], []
for rat in list(superpop.keys()):
    for day in list(superpop[rat].keys()):
        dataset = superpop[rat][day]
        df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
        clustList, clustQuals = util.getClustList(df)

        for unitName, unit in dataset['units'].items():
            quality = clustQuals[clustList.index(unitName)]
            rep.append(unit.repeating)
            all_quals.append(quality)

all_quals = np.asarray(all_quals)
rep = np.asarray(rep)

# %% Remove off track firing. roll this fx into RepCore when done
# The bounds of alleys or intersections are in [[xmin, xmax],[ymin,ymax]] format 

rat, day, unitName = 'R781', 'D4', 'TT3\\cl-maze1.2'
unit = RepCore.loadRepeatingUnit(rat, day, unitName)
raborders = nab.loadAlleyBounds(rat, day)

ontrackSpikesIdx = np.empty((0), dtype=int)
ontrackPositionIdx = np.empty((0), dtype=int)

for regionName, b in ratborders.alleyInterBounds.items():
    ul, ll, ur, lr = [b[0][0], b[1][1]], \
                    [b[0][0], b[1][0]], \
                    [b[0][1], b[1][1]], \
                    [b[0][1], b[1][0]]

    contour = path.Path([ll, lr, ur, ul])
    posIdx = np.where(contour.contains_points(unit.position[:,1:]))[0]
    spkIdx = np.where(contour.contains_points(unit.spikes[:,1:]))[0]

    ontrackPositionIdx = np.concatenate((ontrackPositionIdx, 
                                        posIdx.astype(int)))

    ontrackSpikesIdx = np.concatenate((ontrackSpikesIdx, 
                                        spkIdx.astype(int)))

ontrackPosition = unit.position[ontrackPositionIdx]
ontrackSpikes = unit.spikes[ontrackSpikesIdx]

# %% testing the above 
rat, day, unitName = 'R781', 'D4', 'TT3\\cl-maze1.2'
unit = RepCore.loadRepeatingUnit(rat, day, unitName)
raborders = nab.loadAlleyBounds(rat, day)



# %%
