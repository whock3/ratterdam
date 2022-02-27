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
import matplotlib.path as path
from placeFieldBorders import reorderBorder


def computeSpeed(position, winsz=50):
    """
    Compute 1D speed over time, smoothed
    winsz - window size for rolling average. Units are data samples.
    Returns - 1d vector of speeds at each time point
    """
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
    return sv


def velocity_filtering(position, thresh = Def.velocity_filter_thresh, winsz=50):
    sv = computeSpeed(position, winsz=winsz)
    vf_pos = position[sv > thresh]
    
    return vf_pos

    
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
    
    
def checkMinimumPassesActivity(unit,alley,fr_thresh=2.0,pass_thresh=9):
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

def unitVelocityFilter(ts, position, clust):
    '''
    Remove spikes emitted when the animal's
    velocity was below threshold. 
    
    This class does not explicitly store
    the thresh, or even whether a velocity
    filter was applied. That occurs in the BehavioralData
    class and is reflected implicitly in the position array
    having been filtered. the ts dict is unfiltered and is used
    as a reference.
    '''
    
    allSpikeTs = np.asarray([util.takeClosest(ts, i) for i in clust])
    filtTs = clust[np.isin(allSpikeTs, position[:,0])]
    return filtTs

def filterField(unit, index, rateThresh=0.2, pctThresh=10):
    """
    This was deprecated in early september 2021 in favor of manually curtating
    and redrawing place fields. all fields were processed in this way. see 
    RepCore fx for details. 
    
    Inputs: unit - Unit class object
            index - index of field to filter
            rateThresh - pixel rate that will be used to see % pixels below it
            pctThresh -  % of pixels (0-100) with rate below rateThresh above which field should be discarded
    Returns: True/False, as to whether field passes. True = pass
    Function to filter detected place fields by what percent of the field is 'too quiet'. Count pct
    of pixels with rate below rateThresh. if this pct exceeds pctThresh, return False.
    """
    perim = unit.perimeters[index]
    contour = path.Path(perim)
    spkIn = unit.spikes[contour.contains_points(unit.spikes[:,1:])]
    occIn = unit.position[contour.contains_points(unit.position[:,1:])]
    rm=util.makeRM(spkIn,occIn)
    area = np.sum(~np.isnan(rm.flatten()))
    binsBelowThresh = np.where(rm.flatten()<=rateThresh)[0].shape[0]
    pct = (binsBelowThresh/area)*100
    if pct > pctThresh:
        return False
    else:
        return True
    
def filterFields(unit):
    perims = []
    fields = []
    for i in range(len(unit.perimeters)):
        if filterField(unit, i):
            perims.append(unit.perimeters[i])
            fields.append(unit.fields[i])
    unit.perimeters = []
    unit.fields = []
    unit.perimeters = perims
    unit.fields = fields
    return unit


def normalizeTrajectory(unit, tsStart, tsEnd, contour, alley, border):
    """

    Returns
    -------
    None.

    """
    
    spk_traj = unit.spikes[(unit.spikes[:,0]>tsStart)&(unit.spikes[:,0]<= tsEnd)]
    behav_traj = unit.position[(unit.position[:,0]>tsStart)&(unit.position[:,0]<= tsEnd)]
    
    spk_traj = spk_traj[contour.contains_points(spk_traj[:,1:])]
    behav_traj = behav_traj[contour.contains_points(behav_traj[:,1:])]
    
    spk_sess = unit.spikes[contour.contains_points(unit.spikes[:,1:])]
    behav_sess = unit.position[contour.contains_points(unit.position[:,1:])]
    

    verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
    horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]
    if alley in horizontals:
        xbins = np.linspace(border[0][0],border[0][1],num=15+1) # 0.91 cm bins
        ybins = np.linspace(border[1][0],border[1][1],num=8+1) # 0.91 cm bins 
    elif alley in verticals:
        xbins = np.linspace(border[0][0],border[0][1],num=8+1)
        ybins = np.linspace(border[1][0],border[1][1],num=15+1)
    
    behav_traj_hist = np.histogram2d(behav_traj[:,2],behav_traj[:,1],bins=[ybins,xbins])
    spk_traj_hist = np.histogram2d(spk_traj[:,2],spk_traj[:,1],bins=[ybins,xbins])
    
    behav_sess_hist = np.histogram2d(behav_sess[:,2],behav_sess[:,1],bins=[ybins,xbins])
    spk_sess_hist = np.histogram2d(spk_sess[:,2],spk_sess[:,1],bins=[ybins,xbins])
    
    
    sessrm = (spk_sess_hist[0]/behav_sess_hist[0])*30 # adjust FPS to get Hz 
    trajrm = (spk_traj_hist[0]/behav_traj_hist[0])*30 
    
    return np.nanmean(trajrm/sessrm)
    