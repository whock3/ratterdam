# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 17:09:28 2021

@author: whockei1

Script with functions to filter dataframes containing data from all datasets
Filtering to ensure minimum sampling for what tests you want to do
Datasets are as of script creation superpopulation csvs with visits 
to each visit (each alley comprising field gets its own entry.)
"""
import pandas as pd, numpy as np

def __filterUnitDirSampling(uOriented, direction, passThresh):
    """
    DEPRECATED 2/18/22
    
    Take pd df corresponding to one unit
    It has been filtered into V or H passes already
    Here: ID fields and directions within with enough sampling
    and return them as a new df
    
    passThresh - number of passes for each direction to
    be included

    Returns
    -------
    Df with fields and directions with enough sampling
    
    
    DEPRECATED 2/18/22 - code logic is too convoluted and may not be giving proper results
    Do everything in filterAlleyDatasets in more direct way looping over field ids
    instead of cells and checking fields to be included 

    """
    directions = [['N','S'] if direction=='V' else ['E', 'W']][0]
    fielddirs = {}
    for fn, fg in uOriented.groupby('FieldNum'):
        fielddirs[fn] = {}
        for dn, dg in fg.groupby('CurrDir'):
            fielddirs[fn][dn] = dg.shape[0] 
    includedFields = []
    for fnum, dircounts in fielddirs.items():
        if directions[0] in dircounts.keys() and directions[1] in dircounts.keys():
            if dircounts[directions[0]] >= passThresh and dircounts[directions[1]] >= passThresh:
                includedFields.append(fnum)
    if len(includedFields)>=1:
        uValid = uOriented[uOriented['FieldNum'].isin(includedFields)]
    else:
        uValid = None
    return uValid


def filterAlleyDatasets(df,passThresh=2,filterR=True,filterT=True):
    """
    FOr alleys 
    Take pd df for whole dataset (rats + days)
    And return new df with filtered units
    
    filterR - remove alley traversals that were rewarded
    filterT - remove turnarounds in alleys
    passThresh - min number passes in each direction (of current direction) needed
    """
    
    if filterT:
        df = df[df.Traversal==True]
    if filterR:
        df = df[df.Reward==False]
        
    newDf = []
    for fid, field in df.groupby("FieldID"):
        for orien, ofield in field.groupby("Orientation"):
            dirs = np.unique(ofield.CurrDir)
            if len(dirs)>2:
                print(f"{fid} {orien} has too many current direction types")
            elif len(dirs) == 1:
                pass
            else:
                dA = ofield[ofield.CurrDir==dirs[0]].shape[0]
                dB = ofield[ofield.CurrDir==dirs[1]].shape[0]
                if dA >= passThresh and dB >= passThresh:
                    newDf.append(ofield)
    newDf = pd.concat(newDf)
                
    return newDf
 