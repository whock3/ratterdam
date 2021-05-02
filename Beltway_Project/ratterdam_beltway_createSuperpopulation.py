# -*- coding: utf-8 -*-
"""
Created on Sun May  2 13:34:47 2021

@author: whockei1

Create nested dict with Unit data for all good days defined by a list
No filters applied other than cluster quality. So they can be done differently as needed

"""

import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt
import pandas as pd
from collections import OrderedDict


Def.velocity_filter_thresh = 0
Def.includeRewards = 2

datasets = [("R781","BRD2"),
            ("R781","BRD3"),
            ("R781","BRD4"),
            ("R808","BRD4"),
            ("R808","BRD6"),
            ("R808","BRD7"),
            ("R859","BRD1"),
            ("R859","BRD3"),
            ("R859","BRD5"),
            ("R886","BRD1"),
            ("R886","BRD2")
            ]

for rat,expCode in datasets:
    print(rat, expCode)
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
            valid, acorr, alleys = util.checkInclusion(unit, 3) # 2nd arg to util.checkInclusion is how many comps made per alley. This 
                                                                # value (usually 3) is not being saved here and is defined in the relevant R code so ignore it here
            if valid:
                print(clust)
                unit.acorr = acorr
                unit.validAlleys = alleys
                population[clust] = unit
