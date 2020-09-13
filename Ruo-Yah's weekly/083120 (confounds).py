# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 16:10:36 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
from string import ascii_uppercase
from matplotlib import path
from collections import namedtuple
from copy import deepcopy
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
from alleyTransitions import alleyTransitions
from newAlleyBounds import alleyInterBounds
import ratterdam_Defaults as Def
from utility_fx import weird_smooth, read_pos
from ratterdam_ParseBehavior import adjustPosCamera
import ratterdam_DataFiltering as Filt


def directionFilter(pos):#, spikes):
    """
    Returns position and spikes filtered by which direction the rat was facing

    """
    directions = np.diff(pos[:, 1:3], axis=0)
    allo = np.arctan2(directions[:, 1], directions[:, 0])
    posN = pos[:-1][np.logical_and((np.pi/4)<=allo, allo<(3/4*np.pi))]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allo, allo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allo, allo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allo, allo<(-3/4*np.pi))]
    
    #[spikesN, spikesE, spikesS, spikesW] = [[] for _ in range(4)]
    #for spike in spikes:
    #    if spike[0] in posN[:,0]:
    #        spikesN.append(spike)
    #    elif spike[0] in posE[:,0]:
    #        spikesE.append(spike)
    #    elif spike[0] in posS[:,0]:
    #        spikesS.append(spike)
    #    elif spike[0] in posW[:,0]:
    #        spikesW.append(spike)
    #[spikesN, spikesE, spikesS, spikesW] = [np.array(i) for i in [spikesN, spikesE, spikesS, spikesW]]
    return [posN, posE, posS, posW] #, [spikesN, spikesE, spikesS, spikesW]


def graphDirections(df, suptitle):
    """
    Graphs rate maps (with no firing rate) and bar graphs of direction confounds
    """
    with open(df+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = read_pos(df)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = adjustPosCamera(df, pos, ts)
    position = np.column_stack((ts, posx, posy))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    position = Filt.velocity_filtering(position)
        
    pos = directionFilter(position)
    rbins, cbins = Def.wholeAlleyBins
    rows, cols = np.linspace(0, 480, rbins), np.linspace(0,640, cbins)
    hos = []
    for i in range(len(pos)):
        ho = np.histogram2d(pos[i][:,2],pos[i][:,1],bins=[rows, cols])[0]
        hos.append(deepcopy(ho))
        hos[i][np.where(ho==0)] = np.nan
        hos[i] = weird_smooth(hos[i],Def.smoothing_2d_sigma)
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
    

def occupancy2(pos, title=""):
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
    

def graphOccupancy(df, title, timestamp):
    """
    Load occupancy data without velocity filtering and saves a graph 
    generated by occupancy2
    """
    with open(df+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = read_pos(df)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = adjustPosCamera(df, pos, ts)
    position = np.column_stack((ts, posx, posy))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    
    occupancy2(position, title + " occupancy\nNo velocity filtering")
    plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                + timestamp + " - occupancy " + title + ".png")


Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")