"""
Repetition Manuscript Supplementary Figure Code 
WH 2022-05-10
"""
#%%
import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import ratterdam_Defaults as Def
from scipy.stats import sem
import repetition_manuscript_defaults as MDef



#%% Fig S4 - Supplementary 
df = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Supplementary_Figures\\shuffle_nonlocalDirection\\2022-05-10_shuffle_emmeans.csv")

bins = np.linspace(0,0.3,20)

fig, ax = plt.subplots(1,3,figsize=(10,8))
#proportions are hardcoded from data, see repetition_manuscript_statistics
# and they are current as of 5/10/2022
for i,(dirtype, realprop, cax) in enumerate(zip(["total_current_responsive",
                              "total_previous_responsive",
                              "total_next_responsive"],
                              [37/127,
                              27/160,
                              20/151],
                              fig.axes
                                )):
    cax.hist(df[dirtype],bins=bins,facecolor='grey',edgecolor='black',linewidth=2)
    cax.vlines(np.percentile(df[dirtype], 95), 0, cax.get_ylim()[1],color='k')
    cax.vlines(realprop,0, cax.get_ylim()[1],color='r')
    cax.set_title(dirtype.replace("_"," "), fontsize=25)
    if i == 0:
        cax.set_ylabel("Frequency", fontsize=25)
    if i == 1:
        cax.set_xlabel("Proportion of Cells Responding to Turn Type", fontsize=25)
    cax.spines['right'].set_visible(False)
    cax.spines['top'].set_visible(False)
    cax.tick_params(axis='both', which='major', labelsize=22)

