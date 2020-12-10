# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 16:10:36 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import ImageGrid
from string import ascii_uppercase
from collections import namedtuple
from copy import deepcopy
from bisect import bisect_left
from alleyTransitions import alleyTransitions
from newAlleyBounds import R781, R808, R859
from alleyTransitions import crossBorder, crossBorder2
import ratterdam_Defaults as Def
from ratterdam_ParseBehavior import adjustPosCamera
import ratterdam_DataFiltering as Filt
import utility_fx as util
from RDP4_class import RDP4


def directionFilter(pos):
    """
    Returns position filtered by which direction the rat was facing

    """
    directions = np.diff(pos[:, 1:3], axis=0)
    allo = np.arctan2(directions[:, 1], directions[:, 0])
    posN = pos[:-1][np.logical_and((np.pi/4)<=allo, allo<(3/4*np.pi))]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allo, allo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allo, allo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allo, allo<(-3/4*np.pi))]
    
    return [posN, posE, posS, posW] 


def directionFilterS(pos, spikes):
    """
    Returns position and spikes filtered by which direction the rat was facing

    """
    directions = np.diff(pos[:, 1:3], axis=0)
    allo = np.arctan2(directions[:, 1], directions[:, 0])
    posN = pos[:-1][np.logical_and((np.pi/4)<=allo, allo<(3/4*np.pi))]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allo, allo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allo, allo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allo, allo<(-3/4*np.pi))]
    
    #North
    filtTs = Filt.unitVelocityFilter(pos[:,0], posN, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesN = np.column_stack((filtTs,spikexy))
    #East
    filtTs = Filt.unitVelocityFilter(pos[:,0], posE, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesE = np.column_stack((filtTs,spikexy))
    #South
    filtTs = Filt.unitVelocityFilter(pos[:,0], posS, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesS = np.column_stack((filtTs,spikexy))
    #West
    filtTs = Filt.unitVelocityFilter(pos[:,0], posW, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesW = np.column_stack((filtTs,spikexy))    
    return [posN, posE, posS, posW], [spikesN, spikesE, spikesS, spikesW]


def dirFiltWindow(pos, spikes):
    """
    Returns position and spikes filtered by which direction the rat was facing
    Direction based on average point in small windows
    """
    winsz = 10
    directions = np.diff(pos[:, 1:3], axis=0)
    dirx = [np.mean(directions[0+i:winsz+i, 0]) for i in range(len(directions))]
    diry = [np.mean(directions[0+i:winsz+i, 1]) for i in range(len(directions))]
    allo = np.arctan2(diry, dirx)
    posN = pos[:-1][np.logical_and((np.pi/4)<=allo, allo<(3/4*np.pi))]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allo, allo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allo, allo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allo, allo<(-3/4*np.pi))]
    
    #North
    filtTs = Filt.unitVelocityFilter(pos[:,0], posN, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesN = np.column_stack((filtTs,spikexy))
    #East
    filtTs = Filt.unitVelocityFilter(pos[:,0], posE, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesE = np.column_stack((filtTs,spikexy))
    #South
    filtTs = Filt.unitVelocityFilter(pos[:,0], posS, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesS = np.column_stack((filtTs,spikexy))
    #West
    filtTs = Filt.unitVelocityFilter(pos[:,0], posW, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesW = np.column_stack((filtTs,spikexy))    
    return [posN, posE, posS, posW], [spikesN, spikesE, spikesS, spikesW]


def dirFiltRDP(pos, spikes, epsilon=Def.ptsCm):
    """
    Returns position and spikes filtered by which direction the rat was facing
    Trajectory is simplified with RDP and the direction of filtered out points
    is based on the closest remaining points
    """
    posRDP = RDP4(pos, epsilon).ResultList
    filtTs = Filt.unitVelocityFilter(pos[:,0], posRDP, spikes[:,0])
    #spikesRDP = np.column_stack((filtTs, util.getPosFromTs(filtTs, pos)))
    
    directions = np.diff(posRDP[:, 1:3], axis=0)
    allo = np.arctan2(directions[:, 1], directions[:, 0])
    idxRDP = np.arange(len(pos))[np.isin(pos[:,0], posRDP[:,0])] #indices of pts kept after RDP
    idxFilt = np.arange(len(pos))[~np.isin(pos[:,0], posRDP[:,0])] #indices of pts filtered out with RDP
    idxSmaller = [idxRDP[bisect_left(idxRDP, i)-1] for i in idxFilt] #largest index in idxRDP that is smaller than each index in idxFilt
    allIdx = sorted(list(idxRDP) + list(idxSmaller)) 
    allAllo = np.array([allo[np.where(idxRDP == i)[0]][0] for i in allIdx[:-1]])
    
    posN = pos[:-1][np.logical_and((np.pi/4)<=allAllo, allAllo<(3/4*np.pi)),:]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allAllo, allAllo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allAllo, allAllo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allAllo, allAllo<(-3/4*np.pi))]
    
    #North
    filtTs = Filt.unitVelocityFilter(pos[:,0], posN, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesN = np.column_stack((filtTs,spikexy))
    #East
    filtTs = Filt.unitVelocityFilter(pos[:,0], posE, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesE = np.column_stack((filtTs,spikexy))
    #South
    filtTs = Filt.unitVelocityFilter(pos[:,0], posS, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesS = np.column_stack((filtTs,spikexy))
    #West
    filtTs = Filt.unitVelocityFilter(pos[:,0], posW, spikes[:,0])
    spikexy = util.getPosFromTs(filtTs,pos)
    spikesW = np.column_stack((filtTs,spikexy))    
    return [posN, posE, posS, posW], [spikesN, spikesE, spikesS, spikesW]


def graphDirRatemaps(unit, suptitle, epsilon=Def.ptsCm):
    """
    Graphs ratemaps filtered by direction
    """
    #posDir, spikesDir = directionFilterS(unit.position, unit.spikes)
    #posDir, spikesDir = dirFiltWindow(unit.position, unit.spikes)
    posDir, spikesDir = dirFiltRDP(unit.position, unit.spikes, epsilon)
    ns = [util.makeRM(unit.spikes, unit.position)]
    for i in range(len(posDir)):
        ns.append(util.makeRM(spikesDir[i], posDir[i]))
            
    fig = plt.figure(figsize=(8,9))
    vmax = np.nanpercentile(ns, 98)
    if vmax > 1:
        titles = ["North facing", "East facing", "South facing", "West facing"]
        gs = GridSpec(3, 2, figure=fig)
        ax = fig.add_subplot(gs[0,:])
        ax.set_title("Overall")
        ax.set_xlabel("x coordinates (cm)")
        ax.set_ylabel("y coordinates (cm)")
        im = ax.imshow(ns[0], cmap="jet", origin="lower", vmin=0, vmax=vmax)
        cb = fig.colorbar(im, ax=ax)
        cb.set_label("Rate (Hz)")
        
        for i in range(4):
            ax = fig.add_subplot(gs[i//2+1,i%2])
            ax.set_title(titles[i])
            ax.set_xlabel("x coordinates (cm)")
            ax.set_ylabel("y coordinates (cm)")
            im = ax.imshow(ns[i+1], cmap="jet", origin="lower", vmin=0, vmax=vmax)
        
        fig.suptitle(suptitle+f"\nCutoff = 98th percentile, {round(vmax,1)} Hz", y=1.08)
        fig.tight_layout()
    return fig, vmax
    

def bulkGraphDirRatemaps(units, df, ratDayTetrode, timestamp):
    """
    Makes multiple graphs using graphDirRatemaps
    df = path to the tetrode folder
    """
    qualities = util.cellQuality(df)
    
    for i in range(len(units)):
        if not qualities:
            quality = "None"
        else:
            quality = qualities[str(i+1)]
        fig, vmax = graphDirRatemaps(units[i], f"{ratDayTetrode} 1.{i+1}\nvthresh = 3 cm/s  Cell quality = {quality}")
        if vmax > 1:
            fig.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                    + timestamp + " - " + f"{ratDayTetrode} 1.{i+1}" + ".png",
                    bbox_inches="tight")
        fig.close()


def graphDirections(position, suptitle):
    """
    Graphs occupancy heat maps and bar graphs of direction confounds
    """
    pos = directionFilter(position)
    rbins, cbins = Def.wholeAlleyBins
    rows, cols = np.linspace(0, 480, rbins), np.linspace(0,640, cbins)
    hos = []
    for i in range(len(pos)):
        ho = np.histogram2d(pos[i][:,2],pos[i][:,1],bins=[rows, cols])[0]
        hos.append(deepcopy(ho))
        hos[i][np.where(ho==0)] = np.nan
        hos[i] = util.weird_smooth(hos[i],Def.smoothing_2d_sigma)
        hos[i][np.where(ho==0)] = np.nan
    
    fig, axs = plt.subplots(2,2,figsize=(8,6))
    vmax = np.nanpercentile(hos, 95)
    titles = ["North facing", "East facing", "South facing", "West facing"]
    axs = axs.flatten()
    for i in range(len(pos)):
        axs[i].set_title(titles[i])
        axs[i].set_xlabel("x coordinates (cm)")
        axs[i].set_ylabel("y coordinates (cm)")
        im = axs[i].imshow(hos[i], cmap="jet", origin="lower", vmin=0, vmax=vmax)
    fig.suptitle(suptitle, y=1.04)
    cb = fig.colorbar(im)
    cb.set_label("Occupancy")
    fig.tight_layout()
    
    fig, ax = plt.subplots()
    a = [len(pos[i]) for i in range(len(pos))]
    ax.bar(titles, a)
    ax.set_ylabel("Number of points")
    ax.set_title(suptitle)


def graphDirAndOcc(pos, alleyInterBounds, title):
    """
    Makes 10 graphs: 4 heatmaps of occupancy divided by direction
                    4 bar graphs of points in each alley/intersection divided by direction
                    2 bar graphs comparing north vs south and east vs west
    """
    posDir = directionFilter(pos)
    rbins, cbins = [int(round(480/4.72/2)),int(round(640/4.72/2))]
    rows, cols = np.linspace(0, 480, rbins), np.linspace(0,640, cbins)
    hos = []
    for i in range(len(posDir)):
        ho = np.histogram2d(posDir[i][:,2],posDir[i][:,1],bins=[rows, cols])[0]
        hos.append(deepcopy(ho))
        hos[i][np.where(ho==0)] = np.nan
        hos[i] = util.weird_smooth(hos[i],Def.smoothing_2d_sigma)
        hos[i][np.where(ho==0)] = np.nan
    
    inAlleys = np.zeros((29,4))
    letters = [i for i in ascii_uppercase[:12]]
    a = np.hstack((np.arange(17).astype(str), letters))
    for i in range(17):
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        for j in range(4):
            inAlley = np.all((posDir[j][:,1] > aRect.xmin, posDir[j][:,1] <= aRect.xmax,\
                              posDir[j][:,2] > aRect.ymin, posDir[j][:,2] <= aRect.ymax),axis=0)
            inAlleys[i,j] = np.sum(inAlley)
    for i,letter in enumerate(letters):
        alley = alleyInterBounds[letter]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        for j in range(4):
            inAlley = np.all((posDir[j][:,1] > aRect.xmin, posDir[j][:,1] <= aRect.xmax,\
                              posDir[j][:,2] > aRect.ymin, posDir[j][:,2] <= aRect.ymax),axis=0)
            inAlleys[i+17,j] = np.sum(inAlley)
    inAlleys /= sum(inAlleys)
    
    #fig, axs = plt.subplots(4,2,figsize=(12,12),gridspec_kw={"width_ratios":[3,5]})
    vmax = np.nanpercentile(hos, 99)
    titles = ["North facing", "East facing", "South facing", "West facing"]
    #for i in range(4):
    #    axs[i,0].set_title(titles[i])
    #    axs[i,0].set_xlabel("x coordinates (cm)")
    #    axs[i,0].set_ylabel("y coordinates (cm)")
    #    axs[i,0].set_xticks(np.arange(0,80,20))
    #    axs[i,0].set_xticklabels(np.arange(0,80,20)*2)
    #    axs[i,0].set_yticks(np.arange(0,50,10))
    #    axs[i,0].set_yticklabels(np.arange(0,50,10)*2)
    #    im = axs[i,0].imshow(hos[i], cmap="jet", origin="lower", vmin=0, vmax=vmax)
    #cb = fig.colorbar(im, ax=axs[3,0])
    #cb.set_label("Occupancy (points per bin)")
    
    #for i in range(4):
    #    axs[i,1].bar(np.arange(29), inAlleys[:,i])
    #    axs[i,1].set_title(titles[i])
    #    axs[i,1].set_ylabel("Fraction number of points")
    #    axs[i,1].set_xlabel("Alleys")
    #    axs[i,1].set_xticks(np.arange(29))
    #    axs[i,1].set_xticklabels(a)
    #fig.suptitle(title, y=1.04)
    #fig.tight_layout()
    
    
    fig = plt.figure(figsize=(12,15))
    gs = GridSpec(5, 8, figure=fig)
    for i in range(4): #rate maps
        ax = fig.add_subplot(gs[i, :3])
        ax.set_title(titles[i])
        ax.set_xlabel("x coordinates (cm)")
        ax.set_ylabel("y coordinates (cm)")
        ax.set_xticks(np.arange(0,80,20))
        ax.set_xticklabels(np.arange(0,80,20)*2)
        ax.set_yticks(np.arange(0,50,10))
        ax.set_yticklabels(np.arange(0,50,10)*2)
        im = ax.imshow(hos[i], cmap="jet", origin="lower", vmin=0, vmax=vmax)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Occupancy (points per bin)")
    
    horizVert = ["\nH","\nH","\nV","\nV","\nH","\nV","\nH","\nV","\nH","\nV","\nH","\nV","\nH","\nH","\nV","\nH","\nV"] + ["" for _ in range(12)]
    for i in range(4): #bar graphs
        ax = fig.add_subplot(gs[i, 3:])
        ax.bar(np.arange(29), inAlleys[:,i])
        ax.set_title(titles[i])
        ax.set_ylabel("Fraction number of points")
        ax.set_xlabel("Alleys and Intersections")
        ax.set_xticks(np.arange(29))
        ax.set_xticklabels(np.core.defchararray.add(a, np.array(horizVert)))
    fig.suptitle(title, y=1.04)
    fig.tight_layout()
    
    titles2 = ["North - South", "East - West"]
    for i in range(2): #north vs south and east vs west
        ax = fig.add_subplot(gs[4, i*4:i*4+4])
        ax.bar(np.arange(17), inAlleys[:17,i]-inAlleys[:17,i+2])
        ax.set_title(titles2[i])
        ax.set_ylabel("Fraction number of points")
        ax.set_xlabel("Alleys")
        ax.set_xticks(np.arange(17))
        ax.set_xticklabels(np.core.defchararray.add(np.arange(17).astype(str), np.array(horizVert[:17])))
    
    return fig


def occupancy(pos, title=""):
    """
    Makes a bar graph of the amount of time spent in each alley/intersection
    Uses alleyTransitions' definition of which location the rat is in
    """
    posNew, turns = alleyTransitions(pos)
    
    letters = [i for i in ascii_uppercase[:12]]
    times = np.empty(0)
    for i in np.arange(17).astype(str):
        inLoc = np.where(posNew[:,4] == i)
        entries = np.where(posNew[inLoc,3] == str(1))[1]
        exits = np.where(posNew[inLoc,3] == str(2))[1]
        time = sum(posNew[inLoc][exits,0].astype(float) - posNew[inLoc][entries,0].astype(float))
        times = np.hstack((times,time))
    for i in letters:
        inLoc = np.where(posNew[:,4] == i)
        entries = np.where(posNew[inLoc,3] == str(3))[1]
        exits = np.where(posNew[inLoc,3] == str(4))[1]
        time = sum(posNew[inLoc][exits,0].astype(float) - posNew[inLoc][entries,0].astype(float))
        times = np.hstack((times,time))
    
    totalTime = sum(times)
    times /= totalTime
    fig, ax = plt.subplots(figsize=(8,4))
    ax.bar(list(np.arange(17).astype(str))+letters, times)
    ax.set_title(title)
    ax.set_ylabel("Fraction time in the location")
    ax.set_xlabel("Alleys and intersections")
    

def occupancy2(pos, alleyInterBounds, title=""):
    """
    Checks occupancy using contains_points and pos without velocity filtering
    """
    inAlleys = np.empty(0)
    alleySizes = np.empty(0)
    letters = [i for i in ascii_uppercase[:12]]
    a = np.hstack((np.arange(17).astype(str), letters))
    for i in a:
        alley = alleyInterBounds[i]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                               [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
        alley = path.Path(alleyPerim)
        inAlley = np.sum(alley.contains_points(pos[:,1:3]))
        inAlleys = np.hstack((inAlleys, inAlley))
        alleySize = (aRect.xmax-aRect.xmin) * (aRect.ymax-aRect.ymin)
        alleySizes = np.hstack((alleySizes, alleySize))
    
    inAlleysPts = inAlleys/sum(inAlleys)
    inAlleySize = inAlleys/alleySizes
    #print(inAlleys, alleySizes, inAlleySize)
    fig, axs = plt.subplots(2, 1, figsize=(8,8))
    axs[0].bar(a, inAlleysPts)
    axs[0].set_title(title + "\n\nNormalized by total number of points")
    axs[0].set_ylabel("Fraction number of points in location")
    axs[0].set_xlabel("Alleys and intersections")
    
    axs[1].bar(a, inAlleySize)
    axs[1].set_title("Normalized by alley/intersection size")
    axs[1].set_ylabel("Normalized number of points in location")
    axs[1].set_xlabel("Alleys and intersections")
    fig.tight_layout()
    

def findRewards(df, pos, alleyInterBounds, a=b"0x0002"):
    """
    Find which alley the animal was in and which direction it was going in 
    at the time of rewards
    """
    hrd, rec = util.readNEV(df)
    ts = np.array([i[3] for i in rec if a in i[10]])
    ts = ts[np.where(ts <= pos[-1,0])]
    rewardPos = [pos[bisect_left(pos[:,0], i)] for i in ts]
    print("Number of rewards: ", len(rewardPos))
    rewardPos = np.column_stack((rewardPos, ts))
    inAlleys = np.empty(0)
    inAlleysDir = np.zeros((17, 4))
    axTypes = alleyType(alleyInterBounds) #0 = horiz, 1 = vert
    for i in range(17):
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        
        inAlley = rewardPos[np.all((rewardPos[:,1] > aRect.xmin, rewardPos[:,1] < aRect.xmax,\
                                    rewardPos[:,2] > aRect.ymin, rewardPos[:,2] < aRect.ymax),axis=0)]
        inAlleys = np.hstack((inAlleys, len(inAlley)))
        for reward in inAlley:
            for j in np.arange(0.5,10,0.5):
                end = bisect_left(pos[:,0], reward[3]+j*1e6)
                cross = crossBorder(reward, pos[end], aRect)
                if cross == axTypes[i] or cross == axTypes[i]+2: #going N/S in a horiz alley or going E/W in a vert alley
                    if axTypes[i] == 0: #horizontal alley
                        border1 = np.array([[(aRect.xmin+aRect.xmax)/2,aRect.ymax], [aRect.xmax,aRect.ymax], [aRect.xmax,aRect.ymin], [(aRect.xmin+aRect.xmax)/2,aRect.ymin]]) #East half
                        border2 = np.array([[(aRect.xmin+aRect.xmax)/2,aRect.ymax], [aRect.xmin,aRect.ymax], [aRect.xmin,aRect.ymin], [(aRect.xmin+aRect.xmax)/2,aRect.ymin]]) #West half
                        cross = crossBorder2(reward, pos[end], border1, border2)+1
                        
                            
                        
                    elif axTypes[i] == 1: #vertical alley
                        border1 = np.array([[aRect.xmin,(aRect.ymin+aRect.ymax)/2], [aRect.xmin,aRect.ymax], [aRect.xmax,aRect.ymax], [aRect.xmax,(aRect.ymin+aRect.ymax)/2]]) #North half
                        border2 = np.array([[aRect.xmin,(aRect.ymin+aRect.ymax)/2], [aRect.xmin,aRect.ymin], [aRect.xmax,aRect.ymin], [aRect.xmax,(aRect.ymin+aRect.ymax)/2]]) #South half
                        cross = crossBorder2(reward, pos[end], border1, border2)
                        break
                if cross != 999: #999 = no border crossed
                    break
            if cross != 999:
                inAlleysDir[i,cross] += 1
    return inAlleys, inAlleysDir


def alleyType(alleyInterBounds):
    """
    Finds whether an alley is horizontal(0) or vertical(1)
    """
    axTypes = []
    for i in range(17):
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        if (aRect.xmax - aRect.xmin) > (aRect.ymax - aRect.ymin): #horizontal
            axTypes.append(0)
        elif (aRect.xmax - aRect.xmin) < (aRect.ymax - aRect.ymin): #vertical
            axTypes.append(1)
    return np.array(axTypes)


def graphRewards(df, pos, alleyInterBounds, title="", a=b"0x0002"):
    """
    Makes bar graphs of rewards in each alley and rewards in each alley divided by direction
    """
    inAlleys, inAlleysDir = findRewards(df, pos, alleyInterBounds, a)
    
    fig, axs = plt.subplots(2,1)
    x = np.arange(17)
    width = 0.4
    axs[0].bar(x - 0.5*width, inAlleysDir[:,0], width, label="North")
    axs[0].bar(x - 0.5*width, inAlleysDir[:,1], width, label="East")
    axs[0].bar(x + 0.5*width, inAlleysDir[:,2], width, label="South")
    axs[0].bar(x + 0.5*width, inAlleysDir[:,3], width, label="West")
    axs[0].set_ylabel("Number of rewards")
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(x)
    axs[0].legend()
    
    axs[1].bar(x, inAlleys)
    axs[1].set_ylabel("Number of rewards")
    axs[1].set_xlabel("Alleys")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(x)
    fig.suptitle(title)
    
    return fig
    

def directionTime(pos, posDir, alleyInterBounds):
    """
    Finds the direction the animal was facing in each alley over time
    """
    rbins, cbins = [18, 20]
    cols = np.linspace(pos[0,0], pos[-1,0], cbins)
    hists = []
    axLimits = []
    axTypes = alleyType(alleyInterBounds)+1
    #histArray = np.empty((0,19))
    for i in range(17):
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        axType = axTypes[i]
        if axType == 1: #horizontal
            rows = np.linspace(aRect.xmin, aRect.xmax, rbins)
            axLimit = [aRect.xmin, aRect.xmax]
        elif axType == 2: #vertical
            rows = np.linspace(aRect.ymin, aRect.ymax, rbins)
            axLimit = [aRect.ymin, aRect.ymax]
        
        inAlley1 = np.all((posDir[axType][:,1] > aRect.xmin, posDir[axType][:,1] < aRect.xmax,\
                           posDir[axType][:,2] > aRect.ymin, posDir[axType][:,2] < aRect.ymax),axis=0)
        hist1 = np.histogram2d(posDir[axType][inAlley1,axType],posDir[axType][inAlley1,0],bins=[rows,cols])[0]
        inAlley2 = np.all((posDir[(axType+2)%4][:,1] > aRect.xmin, posDir[(axType+2)%4][:,1] < aRect.xmax,\
                           posDir[(axType+2)%4][:,2] > aRect.ymin, posDir[(axType+2)%4][:,2] < aRect.ymax),axis=0)
        hist2 = np.histogram2d(posDir[(axType+2)%4][inAlley2,axType],posDir[(axType+2)%4][inAlley2,0],bins=[rows,cols])[0]
        
        hists.append([axType,deepcopy(hist1),deepcopy(hist2)])
        hists[i][1][np.where(hist1==0)] = np.nan
        hists[i][1] = util.weird_smooth(hists[i][1],1)#Def.smoothing_2d_sigma)
        hists[i][1][np.where(hist1==0)] = np.nan
        hists[i][2][np.where(hist2==0)] = np.nan
        hists[i][2] = util.weird_smooth(hists[i][2],1)#Def.smoothing_2d_sigma)
        hists[i][2][np.where(hist2==0)] = np.nan
        #histArray = np.vstack((histArray,hists[i][1],hists[i][2]))
        axLimits.append(axLimit)
    
    #vmax = np.nanpercentile(histArray,95)
    return hists, axLimits


def graphDirTime(pos, alleyInterBounds, title=""):
    """
    Graphs direction over time (2D histograms and line graphs)
    """
    #posDir = directionFilter(pos)
    posDir, _ = dirFiltWindow(pos, np.empty((0,3)))
    hists, axLimits = directionTime(pos, posDir, alleyInterBounds)
    xlabels = np.round(np.linspace(pos[0,0]/1e6, pos[-1,0]/1e6, 8), 0).astype(int)
    fig1 = plt.figure(figsize=(12,18))
    for i in range(17):
        axs = ImageGrid(fig1, (i%2*3/4,i/8-i%2/8,0.6,0.6), (1,2), axes_pad=0.1, cbar_mode="single")
        
        ylabels = np.round(np.linspace(axLimits[i][0], axLimits[i][1], 8), 0).astype(int)
        vmax = np.nanpercentile(np.hstack((hists[i][1],hists[i][2])),99)
        for j in range(2):
            im = axs[j].imshow(hists[i][j+1], origin="lower",vmax=vmax)
            axs[j].set_xlabel("Time (s)")
            axs[j].set_xticks(np.linspace(-0.5,18.5,8))
            axs[j].set_xticklabels(xlabels[:-1])
            axs[j].set_yticks(np.linspace(-0.5,16.5,8))
            axs[j].set_yticklabels(ylabels)
        cax = axs.cbar_axes[0]
        axs.cbar_axes[0].colorbar(im)
        axis = cax.axis[cax.orientation]
        axis.label.set_text("Points per bin")
        
        if hists[i][0] == 1: #horizontal
            axs[0].set_ylabel("x coordinates (camera coordinates)")
            axs[0].set_title(f"Alley {int(i)}, East")
            axs[1].set_title(f"Alley {int(i)}, West")
        elif hists[i][0] == 2: #vertical
            axs[0].set_ylabel("y coordinates (camera coordinates)")
            axs[0].set_title(f"Alley {int(i)}, South")
            axs[1].set_title(f"Alley {int(i)}, North")
    fig1.suptitle(title,y=2.43, x=0.68)
    
    fig2, axs = plt.subplots(9,2,figsize=(12,15))
    timeBin = (pos[-1,0]-pos[0,0])/1e6/20
    xlabels2 = np.linspace(pos[0,0]/1e6+timeBin/2, pos[-1,0]/1e6-timeBin/2, 19)
    for i in range(17):
        row, col = [8-(i)//2, 1-(i+1)%2]
        hists[i][1][~np.isfinite(hists[i][1])] = 0
        hists[i][2][~np.isfinite(hists[i][2])] = 0
        medians1 = np.nanmedian(hists[i][1],axis=0)
        medians2 = np.nanmedian(hists[i][2],axis=0)
        if hists[i][0] == 1: #horizontal
            axs[row,col].plot(xlabels2, medians1, label="East")
            axs[row,col].plot(xlabels2, medians2, label="West")
        elif hists[i][0] == 2: #vertical
            axs[row,col].plot(xlabels2, medians1, label="South")
            axs[row,col].plot(xlabels2, medians2, label="North")
        axs[row,col].legend()
        axs[row,col].set_title(f"Alley {int(i)}")
        axs[row,col].set_ylabel("Median pts/bin")
        axs[row,col].set_xlabel("Time (s)")
    fig2.suptitle(title, y=1.02)
    fig2.tight_layout()
    return fig1, fig2

def graphConfounds(dfData, dfGraph, title, rat):
    """
    Puts all confound graphs in a multi-page pdf
    rat: string; "R781", "R808", or "R859"
    """
    if rat == "R781":
        alleyInterBounds = R781.alleyInterBounds
        a = b"0x0040"
    elif rat == "R808":
        alleyInterBounds = R808.alleyInterBounds
        a = b"0x0080"
    elif rat == "R859":
        alleyInterBounds = R859.alleyInterBounds
        a = b"0x0002"
    with open(dfData+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = util.read_pos(dfData)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = adjustPosCamera(dfData, pos, ts)
    position = np.column_stack((ts, posx, posy))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    position = position[np.all(position > np.array([0, 0, 0]), axis=1)]
    posFilt = Filt.velocity_filtering(position, 3)
    
    with PdfPages(dfGraph+title+" confounds.pdf") as pdf:
        plot = graphDirAndOcc(position, alleyInterBounds, title+"\nNo velocity filtering")
        pdf.savefig(plot, bbox_inches="tight")
        plot = graphDirAndOcc(posFilt, alleyInterBounds, title+"\nVelocity filtered with threshold = 3 cm/s")
        pdf.savefig(plot, bbox_inches="tight")
        plot = graphRewards(dfData, position, alleyInterBounds, title+"\nNo velocity filtering", a)
        pdf.savefig(plot, bbox_inches="tight")
        plot1, plot2 = graphDirTime(position, alleyInterBounds, title+"\nNo velocity filtering")
        pdf.savefig(plot1, bbox_inches="tight")
        pdf.savefig(plot2, bbox_inches="tight")
        plot1, plot2 = graphDirTime(posFilt, alleyInterBounds, title+"\nVelocity filtered with threshold = 3 cm/s")
        pdf.savefig(plot1, bbox_inches="tight")
        pdf.savefig(plot2, bbox_inches="tight")
        
Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")


def firingDir(pos, spikes, alleyInterBounds, title=""):
    """
    Plots firing rate vs difference between the median occupancy of opposite directions
    """
    posDir = directionFilter(pos)
    hists, axLimits = directionTime(pos, posDir, alleyInterBounds)
    axTypes = alleyType(alleyInterBounds) #0=h, 1=v
    timeBins = np.linspace(pos[0,0], pos[-1,0], 20)
    tBinSize = (pos[-1,0]-pos[0,0])/19
    fig, axs = plt.subplots(6,3,figsize=(10,16))
    axs = axs.flatten()
    titles2 = ["East - West", "North - South"]
    colors = plt.cm.viridis(np.arange(19)/19)
    samples, diffs, ratess = [], [], []
    for i in range(17): #np.where(axTypes == 0)[0]: #horiz alleys
        hists[i][1][~np.isfinite(hists[i][1])] = 0
        hists[i][2][~np.isfinite(hists[i][2])] = 0
        if axTypes[i] == 0:
            diff = np.median(hists[i][1],axis=0) - np.median(hists[i][2],axis=0) #E-W
        elif axTypes[i] == 1:
            diff = np.median(hists[i][2],axis=0) - np.median(hists[i][1],axis=0) #N-S
        
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        oInAlley = np.all((pos[:,1] > aRect.xmin, pos[:,1] < aRect.xmax,\
                           pos[:,2] > aRect.ymin, pos[:,2] < aRect.ymax),axis=0)
        sInAlley = np.all((spikes[:,1] > aRect.xmin, spikes[:,1] < aRect.xmax,\
                           spikes[:,2] > aRect.ymin, spikes[:,2] < aRect.ymax),axis=0)
        rates = []
        for timeBin in timeBins[:-1]:
            o = pos[oInAlley][np.all((pos[oInAlley][:,0] > timeBin, pos[oInAlley][:,0] <= timeBin+tBinSize),axis=0)]
            s = spikes[sInAlley][np.all((spikes[sInAlley][:,0] > timeBin, spikes[sInAlley][:,0] <= timeBin+tBinSize),axis=0)]
            if len(o) == 0:
                rates.append(np.nan)
            else:
                rates.append(len(s)/len(o)*30)
        
        samples.append(np.logical_or(np.median(hists[i][1],axis=0) != 0, np.median(hists[i][2],axis=0) != 0))
        #print(diff,"\n",rates)
        diffs.append(diff)
        ratess.append(rates)
        #print(diff,"\n",rates)
        
    diffs = np.array(diffs)
    ratess = np.array(ratess)
    xmin = np.min(diffs)
    xmax = np.max(diffs)
    ratemax = np.nanpercentile(ratess,99)
    norm = Normalize(vmin=0, vmax=ratemax)
    for i in range(17):
        #x = discrepancy, y = rate, color = time
        #axs[i].scatter(diffs[i][samples[i]], ratess[i][samples[i]], color=colors[samples[i]])
        #x = time, y = discrepancy, color = rate
        axs[i].scatter(np.arange(19)[samples[i]], diffs[i][samples[i]], color=plt.cm.jet(ratess[i][samples[i]]/ratemax))
        #axs[i].set_xlim(xmin, xmax)
        axs[i].set_title(f"Alley {int(i)}")
        axs[i].set_ylabel(f"{titles2[axTypes[i]]} median pts/bin")#Firing rate (Hz)")
        axs[i].set_xlabel("Start to end of session")#f"{titles2[axTypes[i]]} median pts/bin")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="jet"))
    cb.set_label("Firing rate (Hz)")#"0=start 1=end of session")
    fig.suptitle(title,y=1.02)
    fig.tight_layout()


def firingDirGraphs(unit, title, timestamp):
    """
    Generates and saves a graph using firingDir
    """
    firingDir(unit.position, unit.spikes, R859.alleyInterBounds, title+"\nVelocity filtered with threshold = 3 cm/s")
    plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                    + timestamp + " - " + title + ".png", bbox_inches="tight")
    
    
    