# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:15:28 2018

@author: whockei1
"""

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket, sys
sys.path.insert(0, 'E:\\UserData\\Documents\\GitHub\\ratterdam\\')

import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_Defaults as Def


######################
# Velocity Filtering #
######################

def velocity_filtering(position, thresh = Def.velocity_filter_thresh, winsz=50):
    gradts, gradx, grady = np.gradient(position[:,0]), np.gradient(position[:,1]), np.gradient(position[:,2])
    gradx = [np.mean(gradx[0+i:winsz+i]) for i in range(len(gradx))]
    grady = [np.mean(grady[0+i:winsz+i]) for i in range(len(grady))]
    gradx = np.asarray([i/Def.ptsCm for i in gradx])
    grady = np.asarray([i/Def.ptsCm for i in grady])
    
    vx = np.asarray([1e6*(a/b) for a,b in zip(gradx,gradts)])
    vy = np.asarray([1e6*(a/b) for a,b in zip(grady,gradts)])
    v =  np.sqrt((vx**2)+(vy**2))  
    
    sv = [np.mean(v[0+i:winsz+i]) for i in range(len(v))]
    sv = np.asarray(sv)
    
    vf_pos = position[sv > thresh]
    
    return vf_pos


###########################
# Check minimum activity  #
###########################
    
def checkMiniumUnitActivity(unit, alleyVisits, threshold=50):
    """
    Returns a boolean as to whether the unit emitted as many or 
    more spikes defined by threshold. If so true, else false. Uses first pos 
    in alley 1 (A16) to last pos in A9 (A11). This is because lap starts themselves
    include start box time while textures are being updated.
    """
    spikeCount = 0
    for lap in range(len(alleyVisits[0])):
        begin, end  = alleyVisits[15][lap][0], alleyVisits[10][lap][1]
        spikesInInterval =  unit.spikes[(unit.spikes[:,0]>=begin)&(unit.spikes[:,0]<end)]
        spikeCount += spikesInInterval.shape[0]
    if spikeCount >= threshold:
        return True
    else:
        return False
    
    
def checkMinimumPassesActivity(unit,alley,fr_thresh=1.0,pass_thresh=9):
    """"
    Filtering function based on N passes that have >=Y FR
    Check each pass to see if its nanmax exceeds the fr_thresh
    If at least pass_thresh # do, then then return TRUE else FALSE
    
    
    Pass thresh default was 3 until 20200606 when it became 9
    """
    pass_thresh_count = 0
    for visit in unit.alleys[alley]:
        maxfr = np.nanmax(visit['ratemap1d'])
        if maxfr >= fr_thresh:
            pass_thresh_count += 1
    if pass_thresh_count >= pass_thresh:
        return True
    else:
        return False
    
def checkOverallAlleyActivity(unit, alley, thresh=1.):
    """
    Filtering function based on overall 1d alley RM
    Check to see that overall 1d RM has a mean value at least thresh Hz
    Alley is 1 indexed (all unit object alley indices are)
    """
    mean = np.nanmean(unit.alleyRMS[alley]['overall'])
    if mean >= thresh:
        return True
    else:
        return False