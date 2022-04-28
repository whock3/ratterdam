# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:48:57 2022

@author: whockei1

Script to break down directionality by track region and repetition status
3 types of labels:
    - Repeating cells versus single-fielded cells (ignoring, for now, multi-field nonrepeating)
    - Alleys versus intersections
    - Interior vs Perimeter
8 total groups. Will compare them individually or collapsing across different labels

This script combines the functionality of repetition_AlleysVsIntersecxtions.py 
and repetition_InteriorVsPerimeter.py, as well as distinguishing repeating type
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt, copy

def filter_alleyDf(alleydf, toggles):
    #Repeating toggle
    rep = toggles['repeating_type']
    if rep == 'R':
        alleydf = alleydf[alleydf.Repeating==True]
    elif rep == 'SF':
        alleydf = alleydf[alleydf.NumFields==1]
    elif rep == 'All':
        pass
    else:
        print("Toggle error")


    #Region location toggle
    loc = toggles['region_location']
    if loc == 'P':
        alleydf = alleydf[alleydf.Location=='P']
    elif loc == 'I':
        alleydf = alleydf[alleydf.Location=='I']
    elif loc == 'All':
        pass
    else:
        print("Toggle error")
        
    return alleydf

def filter_interDf(interdf, toggles):
    #Repeating toggle
    rep = toggles['repeating_type']
    if rep == 'R':
        interdf = interdf[interdf.Repeating==True]
    elif rep == 'SF':
        interdf = interdf[interdf.NumFields==1]
    elif rep == 'All':
        pass
    else:
        print("Toggle error")


    #Region location toggle
    loc = toggles['region_location']
    if loc == 'P':
        interdf = interdf[~((interdf.Inters=='F')|(interdf.Inters=='G'))]
    elif loc == 'I':
        interdf = interdf[(interdf.Inters=='F')|(interdf.Inters=='G')]
    elif loc == 'All':
        pass
    else:
        print("Toggle error")


    return interdf


def calculate_subgroup_directionality(alleydf, interdf, toggles):
    """
    Take dataframes from alleys and intersections. These have been filtered
    for visits that "go into" the field enough, but have not been filtered
    for sufficient sampling in each direction at a given track compartment
    
    Each dataframe is filtered to only include data from certain track locations
    or cell types
    
    toggles keys:
    'repeating_type':   'R' - repeating cells
                        'SF' - single fielded cells
                        'All' - all cells
    'region_type':      'A' - alleys
                        'I' - intersections
                        'All' - alleys and intersections (i.e. everything)
    'region_location':  'P' - perimeter regions
                        'I' - interior regions
                        'All' - all regions
                        
    Having filtered the data, the directionality of each field is computed.
    For alleys it is abs(mean(dir X) - mean(dir Y)). For intersections it is
    max(mean for each direction) - min(mean for each direction). 
    
    For intersections, current direction is actually egocentric turn, e.g. L,R,B,S
    since turning through an intersection involves 2 headings 
    
    Returns vector of all the fields' mean diffs under the given filters
    """

    interdiffs = []
    alleydiffs = []
    
    f_alleydf = filter_alleyDf(alleydf, toggles)
    f_interdf = filter_interDf(interdf, toggles)
    
    # Field Chunk Directionality for Alleys. If field overlaps multiple
    # intersections, analyze each piece separately 
    for fid, field in f_alleydf.groupby("FieldID"):
        #field can spill over multiple alleys within perimeter or interior alleysets
        oriens = np.unique(field.Orientation)
        for o in oriens:
            try:
                ofield = field[field.Orientation==o]
                dirs = np.unique(ofield.CurrDir)
                diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                alleydiffs.append(diff)
            except:
                pass
        
        
    # Field Chunk Directionality for Intersections. If field overlaps multiple
    # intersections (less likely than for alleys), analyze each piece separately 
    
    # To compute intersection directionality, take the direction with the max firing
    # and subtract the FR under the opposite direction. So it's a fair comparison
    # to alleys. And the direction you use needs to be chosen, here using prevdir
    opposite_dir_lookup = {'N':'S',
                           'E':'W',
                           'W':'E',
                           'S':'N'
                            }
    
    for fid, field in f_interdf.groupby("FieldID"):
    
        # most fields will overlap only 1 intersection (if any), but some are long enough to overlap 2
        inters = np.unique(field.Inters)
        for inter in inters:
            ifield = field[field.Inters==inter]
            dirs = np.unique(ifield.PrevDir)
            if len(dirs)>1:
                meanDirDiffs = {}
                for direction in dirs:
                    meanDirDiffs[direction] = ifield[ifield.PrevDir==direction].Rate.mean()
                maxdir = max(meanDirDiffs, key = lambda x: meanDirDiffs[x])
                maxdirrate = meanDirDiffs[maxdir]
                
                #try to get opposite dir FR, if that doesnt exist just pick randomly
                try:
                    mindirrate = meanDirDiffs[opposite_dir_lookup[maxdir]]
                except: 
                    # this convoluted line randomly selects a direction to be 
                    # compared with the max FR direction to get a directionality value
                    # the reshaping business is bc choice needs 1d and argwhere returns 2d
                    mindir = dirs[np.random.choice(np.argwhere(dirs!=maxdir).reshape(np.argwhere(dirs!=maxdir).shape[0],))]
                    mindirrate = meanDirDiffs[mindir]
                    

                diff = maxdirrate - mindirrate
                interdiffs.append(diff)

        
    loc = toggles['region_type']
    if loc == 'A':
        diffs = alleydiffs
    elif loc == 'I':
        diffs = interdiffs
    elif loc == 'All':
        diffs = alleydiffs + interdiffs
        
    return diffs
    