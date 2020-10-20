# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 11:38:14 2018

@author: whockei1

Ratterdam Analysis
Directionality Analysis

Divide alley visit data by whether the animal entered from the NE or SW side
depending on orientation of alley. Plot various routines such as linear rate
maps divided into directional groups. Note that instantaneous direction is
not taken into account, i.e. if he enters from one side it goes into that group
even if there are data points in that visit where he turns and faces other direction.

This file will read a script of units to analyze, loop over them, and either
plot them live or save them to a Google drive / ratterdam / figures dir. 

Definition of direction is whether the first occupancy point of a visit is closer
to the SW corner or NE corner of the alley based on Euclidian distance calc
by numpy.linalg.norm(). 
"""

############
# Imports #
###########
import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket, os
from numpy.linalg import norm as npNorm
from scipy.stats import sem
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as GS
from matplotlib.backends.backend_pdf import PdfPages
from importlib import reload
import matplotlib.ticker as plticker


import utility_fx as util
import ratterdam_ParseBehavior as pBehav
import ratterdam_CoreDataStructures as core
import ratterdam_visBasic as Vis
from ratterdam_Defaults import *


#############
# Functions #
#############

def extractCorners(givenAlleyBounds):
    """Alley bounds gives [[x1, x2], [y1, y2]]. Convert that
    to UL, LL, UR, LL (xn, ym) points n,m <- [1,2]
    ul - x1, y2
    ll - x1, y1
    ur = x2, y2
    lr = x2, y1
    
    Returns ul, ll, ur, lr
    """
    b = givenAlleyBounds # for ease of typing
    ul, ll, ur, lr = [b[0][0], b[1][1]], [b[0][0], b[1][0]], [b[0][1], b[1][1]], [b[0][1], b[1][0]]
    return ul, ll, ur, lr

def checkCloserPoint(p1, p2, pt):
    """
    Given two points p1, p2
    where each is [x,y]
    see which pt is closer to
    (also of form [x,y])
    
    Return string "first" or "second"
    meaning its closer to first point arg
    or second point arg. If equal return "error"
    """
    d1 = npNorm(p1 - pt)
    d2 = npNorm(p2 - pt)
    if d1 < d2:
        return "first"
    elif d2 < d1:
        return "second"
    else:
        return None
    
    
def checkVisitEntrySide(visitOccs, bounds):
    """
    visitOccs is [ts,x,y] arr for 1 visit
    bounds is [ul, ll, ur, lr] for alley in question
    Return a label "SW" or "NE"
    """
    begin = visitOccs[0,1:]
    ll, ur = bounds[1], bounds[2]
    closerPt = checkCloserPoint(ll, ur, begin)
    if closerPt is not None:
        if closerPt == "first":
            side = "SW"
        elif closerPt == "second":
            side = "NE"
    else:
        side = None
    return side


def groupVisitsByDir(alley):
    """
    Given a 1-idx alley, consider all
    visits and group them by whether entry
    was from SW or NE side.
    
    Method checks if 1st occ pt is closer
    to LL or UR corner of alley
    
    Returns dict list of visits from SW, NE
    """
    bounds = extractCorners(alleyBounds[alley-1]) # ul, ll, ur, lr format list
    visitDirs = {"SW":[], "NE":[]}
    for i,visit in enumerate(unit.alleys[alley]):
        side = checkVisitEntrySide(visit['occs'], bounds)
        if side is not None:
            visitDirs[side].append(i)
    
    return visitDirs

def visitsByDir_LinearRM(alley):
    """
    Given an alley, separate visits by
    entry on the SW, NE side of it.
    Concat all those linear RM into
    a nXc matrix
    
    Assume unit is in local namespace
    defaults in local namespace
    """
    c = singleAlleyBins[1]-1 #-1 because the nominal val is +1 for histgr reasons
    visitDirs = groupVisitsByDir(alley)
    groupedLinRMs = {"SW":np.empty((0,c)), "NE":np.empty((0,c))}
    for side in ["SW", "NE"]:
        if visitDirs[side] is not []:
            for visitIdx in visitDirs[side]:
                groupedLinRMs[side] = np.vstack((groupedLinRMs[side], unit.alleys[alley][visitIdx]['ratemap1d']))
    return groupedLinRMs


def plot_AvgDirLinRM(ax, alley, groupedLinRMs):
    """
    Given a dict of matrices, one for all
    visits from each side of alley (SW, NE)
    plot them avg +/- sem on same subplot
    provided externally by ax arg
    
    make sure scipy.stats.sem is imp as 'sem'
    
    dir1 will be blue, dir2 red. Don't confuse w/ longrunning txt color codes
    """
    sides = list(groupedLinRMs.keys()) # initially side closer to SW, NE corner but make arb if changed
    dir1, dir2 = groupedLinRMs[sides[0]], groupedLinRMs[sides[1]]
    dir1mean, dir2mean, dir1sem, dir2sem = np.nanmean(dir1, axis=0), np.nanmean(dir2, axis=0), sem(dir1,axis=0, nan_policy='omit'), sem(dir2,axis=0, nan_policy='omit')
    for mean, err, color in zip([dir1mean, dir2mean], [dir1sem, dir2sem], ['b', 'r']):
        ax.plot(mean, color)
        ax.fill_between(range(len(mean)), mean-err, mean+err, color=color, alpha=0.5)
    ax.set_title(f"Alley {alley}, {sides[0]} (b): {dir1.shape[0]}, {sides[1]} (r): {dir2.shape[0]}",fontsize=12)
    ax.set_xticks([])
    
    #make textbox of counts, passes thru dir by txt present
    # counts is nested dir txt -> side
    counts = tabulateTxtByDir(getTxtVisitsByDir(alley))
    annot = "\n".join((
        f"A: SW: {counts['A']['SW']}, NE:{counts['A']['NE']}",
        f"B: SW: {counts['B']['SW']}, NE:{counts['B']['NE']}",
        f"C: SW: {counts['C']['SW']}, NE:{counts['C']['NE']}",
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, annot, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    plt.show()
    

def plotRoutine_AvgDirLinRM():
    """
    """
    fig, ax = plt.subplots(5,4, figsize=(10,10))
    for alley in range(1,18):
        axis = fig.axes[alley-1]
        groupedLinRMs = visitsByDir_LinearRM(alley)
        plot_AvgDirLinRM(axis, alley, groupedLinRMs)
    plt.suptitle(f"{exp} {unit.name}")
    plt.show()
    
def getTxtVisitsByDir(alley):
    """
    Given an alley, group visits
    whether SW/NE entry and crossref
    with txts present to get table
    of txts by dir count
    """
    txtDirTable = {txt:{"SW":[], "NE":[]} for txt in ["A", "B", "C"]}
    visitDirs = groupVisitsByDir(alley)
    for side in ["SW", "NE"]:
        if visitDirs[side] is not []:
            for visitIdx in visitDirs[side]:
                txt = unit.alleys[alley][visitIdx]['metadata']['stimulus'][0] # 0 bc it has dtype as entry 1
                txtDirTable[txt][side].append(visitIdx)
    return txtDirTable

def tabulateTxtByDir(txtDirTable):
    """
    Helper fx to count passes along a dir
    by txt present
    """
    counts = {txt:{"SW":0, "NE":0} for txt in ["A", "B", "C"]}
    for txt in ["A", "B", "C"]:
        for side in ["SW", "NE"]:
            if txtDirTable[txt][side] is not []:
                counts[txt][side] = len(txtDirTable[txt][side])
    return counts


def groupTrials(alley, trialList):
    """
    Given a list of visits and alley
    gather them and vstack
    If empty return None
    """
    if trialList ==[] or trialList == None:
        return None
    else:
        trialMat = np.empty((0, singleAlleyBins[1]-1))
        for trial in trialList:
            rm = unit.alleys[alley][trial]['ratemap1d']
            trialMat = np.vstack((trialMat, rm))
    return trialMat

def plot_VisitsDirLinRM(alley):
    """
    Create a series of heatmap/imshow linear ratemaps for an alley
    Each set corresponds to passes from a given direction over a given txt
    So its a 3 x 2 array (A,B,C) x (NE, SW).
    
    Takes alley as argument. Unit must be in fx namespace.
    """
    table = getTxtVisitsByDir(alley)
    print(table)
    
    fig, ax = plt.subplots(3,2, figsize=(7,7))
    i=0 # crude counter to get current idx for axes
    
    allmats = np.empty((0, singleAlleyBins[1]-1)) # gather all so you can get an overall max.
    
    for txt in ["A", "B", "C"]:
        for side in ["SW", "NE"]:
            mat = groupTrials(alley, table[txt][side])
            if type(mat) is np.ndarray:
                allmats = np.vstack((allmats, mat))
    mymax = util.calcSmartMax(allmats, cutoff=0.90, scale=2.5, bins=100)   
    
    for txt in ["A", "B", "C"]:
        for side in ["SW", "NE"]:
            axis = fig.axes[i]
            mat = groupTrials(alley, table[txt][side])
            if type(mat) is np.ndarray:
                im = axis.imshow(mat, interpolation='None', aspect='auto', cmap = cmap, vmax = mymax)
                axis.set_yticks(range(len(table[txt][side])))
            
                axis.set_yticklabels(table[txt][side])
            i +=1 #even if its empty still increment to move on. 
            
    # annotate with labels
    fig.axes[0].set_title("SW Entry")
    fig.axes[1].set_title("NE Entry")
    fig.axes[0].set_ylabel("A")
    fig.axes[2].set_ylabel("B")
    fig.axes[4].set_ylabel("C")
    
    plt.suptitle(f"{exp} {unit.name} Alley {alley}")
    plot.add_colorbar(fig,im)
    plt.show()
    
    
###################################################################
# Analysis Routine - Directionality for Unit, and Alley w Effect  #
###################################################################


if __name__ == '__main__':

    units = util.readUnitsToLoad()
    cmap = util.makeCustomColormap()
    
    figpath = codeDirBase + "\\KnierimLab\\Ratterdam\\Figures\\Directionality\\"
    if not os.path.isdir(figpath):
        os.mkdir(figpath)
        
    plot = Vis.BasicRateMaps()
    
    
    for u in units:
        
        rat = u['rat']
        exp = u['day']
        clustname = u['cluster']
        
        if 'alley' in u:
            alley = int(u['alley'])
    
        datafile = f"{dataDrive}\\{rat}\\{rat}{exp}\\"
        behav = core.BehavioralData(datafile, exp, velocity_filter_thresh)
        ts, position, alleyTracking, alleyVisits,  txtVisits = behav.loadData()
        unit = core.UnitData(clustname, datafile, exp, alleyBounds, alleyVisits, txtVisits, position, ts)
        unit.loadData_raw()
        
        with PdfPages(f"{figpath}{rat}_{exp}_{unit.name}_vthresh{velocity_filter_thresh}__A{alley}_directionality.pdf") as pdf:
            
            fig,ax = plt.subplots()
            plot.plot_linAlley(fig.axes[0], unit, alley)
            try:
                pdf.savefig()
            except:
                pass
            plt.close()
            
            plot_VisitsDirLinRM(alley)
            
            try:
                pdf.savefig()
            except:
                pass
            plt.close()
            
            plotRoutine_AvgDirLinRM()
            
            try:
                pdf.savefig()
            except:
                pass
            plt.close()
        
        
        
    
