# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 10:43:59 2021

@author: whockei1

This script loads Beltway data from single recording day into a dictionary.
Data is thresholded based on activity (see util.checkInclusion, in brief it
relies on # passes with sufficient activity and an average field envelope that is
certain # bins long) and quality (3+ on 1-5 scale, 5=best).

This data is converted into a "longform" pandas frame where each firing rate
in a spatial bin is associated with a spatial bin, texture, trial, reward status,
and alley. This is saved as csv in E:\\Ratterdam\\R_data for use in (G)LMER 
analysis in R. R code in E:\\UserData\\Documents\\GitHub\\ratterdam\\\Beltway_Project\\

"""


import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt
import pandas as pd
from collections import OrderedDict

rat = 'R859'
expCode = 'BRD5'
datafile = f'E:\\Ratterdam\\{rat}\\{rat}{expCode}\\'

qualThresh = 3
alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
clustlist, clustQuals = util.getClustList(datafile) # clustList and clustQuals share same order. ith entry in each is the name and qual of same cell. 
population = OrderedDict()

# Load Data
for i,clust in enumerate(clustlist):
    if clustQuals[i] >= qualThresh:
        unit = Core.UnitData(clust, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
        unit.loadData_raw()
        validalleys = []
        valid, acorr, alleys = util.checkInclusion(unit, 3, 10,fieldCheck=False) # 2nd arg to util.checkInclusion is how many comps made per alley. This 
                                                            # value (usually 3) is not being saved here and is defined in the relevant R code so ignore it here
        if valid:
            print(clust)
            unit.acorr = acorr
            unit.validAlleys = alleys
            population[clust] = unit



# Create "Long Form" of Data. Single Day
#  Nested Structure: Cell > trial > texture > reward > spatial bin. 
            
#%% Create dataframe 
            
# Reassignment fx, toggles and params
alleyReassignmentSchedule = {1:[2,4,6,8,9],
                            2:[1,4,6,8,9],
                            4:[1,2,6,8,9],
                            6:[1,2,4,8,9],                           
                            8:[1,2,4,6,9],
                            9:[1,2,4,6,8]}
            
def getReassignmentData(unit, alley, reassignmentOrder, i):
    newalley = alleyReassignmentSchedule[Def.beltwayAlleyLookup[alley]][reassignmentOrder]
    newvisit = unit.alleys[alleyLookup_rev[newalley]][i] 
    return newvisit['metadata']['stimulus'], newvisit['metadata']['reward']            
            
doingReassignment = False
reassignmentOrder = 4
alleyLookup_rev = {v:k for k,v in Def.beltwayAlleyLookup.items()}

if doingReassignment == True:
    cells, cellNames, trials, textures, rewards, spatialBins, alleys, firingRate = [], [], [], [], [], [], [], []
    nbins = Def.singleAlleyBins[0]-1
    for cellID, (unitname, unit) in enumerate(population.items()):
        print(unitname)
        for alley in unit.validAlleys:
            if alley in [16,17,1,7,10,11]:
                for i,visit in enumerate(unit.alleys[alley]):
                    reTxt, reR = getReassignmentData(unit, alley, reassignmentOrder, i)
                    alleys.extend([alley]*nbins)
                    cells.extend([cellID]*nbins)
                    cellNames.extend([unitname]*nbins)
                    trials.extend([i]*nbins)
                    textures.extend([reTxt]*nbins)
                    rewards.extend([reR]*nbins)
                    spatialBins.extend(range(nbins))
                    firingRate.extend(visit['ratemap1d'])
else:
    cells, cellNames, trials, textures, rewards, spatialBins, alleys, firingRate = [], [], [], [], [], [], [], []
    nbins = Def.singleAlleyBins[0]-1
    for cellID, (unitname, unit) in enumerate(population.items()):
        print(unitname)
        for alley in unit.validAlleys:
            for i,visit in enumerate(unit.alleys[alley]):
                alleys.extend([alley]*nbins)
                cells.extend([cellID]*nbins)
                cellNames.extend([unitname]*nbins)
                trials.extend([i]*nbins)
                txt = visit['metadata']['stimulus']
                textures.extend([txt]*nbins)
                rewards.extend([visit['metadata']['reward']]*nbins)
                spatialBins.extend(range(nbins))
                firingRate.extend(visit['ratemap1d'])
            
firingRate = np.log(np.asarray(firingRate)+1)
data = {'rate':firingRate, 'cell':cells, 'name':cellNames, 'alley': alleys, 'trial':trials, 'texture':textures, 'reward':rewards, 'spatialBin':spatialBins}
alldata = pd.DataFrame(data=data)
alldata.dropna(inplace=True) 

# Remove spatial bins which do not have enough sampling across textures to
# accurately estimate the mean firing rate. Rationale 5/10/21 is that low # samples
# that co-occur w firing on those trials (eg scanning at the end of an expanded alley bounds)
# gives low/no variance or CIs and trivially a CI overlap check returns a pass.
samplethresh = 5
for unitName in alldata['name'].unique():
    udf = alldata.loc[alldata['name']==unitName]
    for alley in udf['alley'].unique():
        for sb in range(Def.singleAlleyBins[0]-1):
           na = udf.loc[(udf['alley']==alley) & (udf['spatialBin']==sb) & (udf['texture']=='A')].shape[0] 
           nb = udf.loc[(udf['alley']==alley) & (udf['spatialBin']==sb) & (udf['texture']=='B')].shape[0]
           nc = udf.loc[(udf['alley']==alley) & (udf['spatialBin']==sb) & (udf['texture']=='C')].shape[0] 
           if na >= samplethresh and nb >= samplethresh and nc >= samplethresh:
               pass
           else:
               alldata = alldata.drop(alldata[(alldata['name']==unitName)&(alldata['alley']==alley)&(alldata['spatialBin']==sb)].index)

        

#save data
stamp = util.genTimestamp()
filename = f"{stamp}_{rat}{expCode}_{Def.velocity_filter_thresh}vfilt_\
{Def.smoothing_1d_sigma}stepsmooth_{Def.singleAlleyBins[0]-1}bins_\
{Def.includeRewards}R_{qualThresh}qual.csv"
               
alldata.to_csv(f"E:\\Ratterdam\\R_data_beltway\\{filename}", header=True, index=False)