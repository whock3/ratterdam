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

def filterUnitDirSampling(uOriented, direction, passThresh=3):
    """
    Take pd df corresponding to one unit
    It has been filtered into V or H passes already
    Here: ID fields and directions within with enough sampling
    and return them as a new df
    
    passThresh - number of passes for each direction to
    be included

    Returns
    -------
    Df with fields and directions with enough sampling

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


def filterAlleyDatasets(df,filterFx):
    """
    FOr alleys 
    Take pd df for whole dataset (rats + days)
    Apply filterFx to each cell within it 
    And return new df with filtered units
    """
    newDf = []
    for cellid in df['CellID'].unique():
        u = df[df['CellID']==cellid]
        for d in ['V','H']:
            uOriented = u[u['Orientation']==d]
            filtdf = filterUnitDirSampling(uOriented,d)
            if filtdf is not None:
                newDf.append(filtdf)
    newDf = pd.concat(newDf)
    return newDf
 