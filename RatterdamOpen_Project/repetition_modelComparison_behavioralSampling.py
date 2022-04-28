# -*- coding: utf-8 -*-
"""
Created on Wed April 27th 2022

@author: whockei1


Bootstrapping GLM Model Comparison
Looking at effects of behavioral sampling on firing. 
Specifically, 
    - total passes through field
    - number unique types of paths (Prev+Curr+Next Directions)
    - Directional Bias
Directional bias is more of a confound of the other two as it can be thought of as a readout of preferred/anti-preferred sampling
which would less interestingly drive change in FR (of directional field)
"""

import pandas as pd, matplotlib.pyplot as plt, numpy as np
import utility_fx as util 
import statsmodels.api as sm
import statsmodels.formula.api as smf


rat_list = ['R765',
            'R765',
            'R781', 
            'R781', 
            'R808', 
            'R808', 
            'R859', 
            'R859', 
            'R886', 
            'R886']

day_list = ['RFD5',
            'DFD4',
            'D3', 
            'D4',
            'D6',
            'D7',
            'D1',
            'D2',
            'D1',
            'D2']

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220404-210901_superPopAlleyBehaviorResponse_1.5vfilt_FieldNormedTrue.csv"
alleydf = pd.read_csv(alleydatapath)

if 'Code' not in alleydf.columns:
    codes = []
    for r, row in alleydf.iterrows():
       code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
       codes.append(code)
    alleydf = alleydf.assign(Code=codes)

for rat,day in zip(rat_list, day_list):
    
    rdf = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]

    #### Time windowing 
    
    included_field_chunks = [] # list of tuples: [(fid,alley)]
    
    rdf.StartTimes = (rdf.StartTimes - rdf.StartTimes.min())/1e6 # s, ref'd to 1st sample
    
    window = 10*60
    offset = 2*60
    wins = []
    begin = 0
    stop = False
    
    while not stop:
        a,b = begin, begin + window
        if b < np.ceil(rdf.StartTimes.max()):
            wins.append((a,b))
            begin += offset
        else:
            stop = True
            
            
    numpaths_time = []
    directionality_time = []
    bias_time = []
    totalpaths_time = []
    
    for win in wins:
        
        numpaths_win = []
        directionality_win = []
        bias_win = []
        totalpaths_win = []

        windf = rdf[(rdf.StartTimes>win[0])&(rdf.StartTimes<=win[1])]

        for fid, field in windf.groupby("FieldID"):
            
            for o, ofield in field.groupby("Orientation"):
    
                dirs = np.unique(ofield.CurrDir)
                if len(dirs)>1:
                    included_field_chunks.append((np.unique(ofield.FieldID)[0], np.unique(ofield.Alleys)[0]))
                    diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                    directionality_win.append(diff)
                    pathcount = 0
                    for cname, code in ofield.groupby("Code"):
                        pathcount +=1
                    numpaths_win.append(pathcount)
                    totalpaths_win.append(ofield.shape[0]) # this is total # filtered passes thru field. not overall # passes thru alley.
                    ndira = ofield[ofield.CurrDir==dirs[0]].shape[0]
                    ndirb = ofield[ofield.CurrDir==dirs[1]].shape[0]
                    bias_win.append(max(ndira,ndirb)/ofield.shape[0])
                
                    
                
        numpaths_time.append(numpaths_win)
        totalpaths_time.append(totalpaths_win)
        directionality_time.append(directionality_win)
        bias_time.append(bias_win)
    
    # So the data above were calculated within time windows and now 
    # we are pooling all that data together. 
    diffs  = np.asarray([i for j in directionality_time for i in j])
    ntrajs = np.asarray([i for j in numpaths_time for i in j])
    totaltrajs = np.asarray([i for j in totalpaths_time for i in j])
    bias = np.asarray([i for j in bias_time for i in j])
    included_field_chunks = np.asarray(included_field_chunks)      
    
    print(f"==============={rat} {day} GLM Model Results=================")


    model_data = pd.DataFrame(data={
        'Rate':diffs+1,
        'Bias':bias,
        'TotalPaths':totaltrajs,
        'UniquePaths':ntrajs
    })

    formula = "Rate ~ TotalPaths + UniquePaths + Bias"
    mod = smf.glm(formula=formula, data=model_data, family=sm.families.Gaussian()).fit()
    print(formula)
    print(mod.summary())

    formula = "Rate ~ Bias"
    mod = smf.glm(formula=formula, data=model_data, family=sm.families.Gaussian()).fit()
    print(formula)
    print(mod.summary())

    formula = "Rate ~ TotalPaths"
    mod = smf.glm(formula=formula, data=model_data, family=sm.families.Gaussian()).fit()
    print(formula)
    print(mod.summary())

    formula = "Rate ~ UniquePaths"
    mod = smf.glm(formula=formula, data=model_data, family=sm.families.Gaussian()).fit()
    print(formula)
    print(mod.summary())

