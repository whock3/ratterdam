# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 16:10:36 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import ImageGrid
from string import ascii_uppercase
from collections import namedtuple
from copy import deepcopy
from bisect import bisect_left
from alleyTransitions import alleyTransitions
from newAlleyBounds import alleyInterBounds
from alleyTransitions import crossBorder
import ratterdam_Defaults as Def
from ratterdam_ParseBehavior import adjustPosCamera
import ratterdam_DataFiltering as Filt
import utility_fx as util


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


def graphDirections(position, suptitle):
    """
    Graphs rate maps (with no firing rate) and bar graphs of direction confounds
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


def graphDirAndOcc(pos, title):
    """
    Makes 8 graphs: 4 rate maps (with no firing rate) of occupancy divided by direction
                    4 bar graphs of points in each alley/intersection divided by direction
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
    
    fig, axs = plt.subplots(4,2,figsize=(12,12),gridspec_kw={"width_ratios":[3,5]})
    vmax = np.nanpercentile(hos, 99)
    titles = ["North facing", "East facing", "South facing", "West facing"]
    for i in range(4):
        axs[i,0].set_title(titles[i])
        axs[i,0].set_xlabel("x coordinates (cm)")
        axs[i,0].set_ylabel("y coordinates (cm)")
        axs[i,0].set_xticks(np.arange(0,80,20))
        axs[i,0].set_xticklabels(np.arange(0,80,20)*2)
        axs[i,0].set_yticks(np.arange(0,50,10))
        axs[i,0].set_yticklabels(np.arange(0,50,10)*2)
        im = axs[i,0].imshow(hos[i], cmap="jet", origin="lower", vmin=0, vmax=vmax)
    cb = fig.colorbar(im, ax=axs[3,0])
    cb.set_label("Occupancy (points per bin)")
    
    for i in range(4):
        axs[i,1].bar(np.arange(29), inAlleys[:,i])
        axs[i,1].set_title(titles[i])
        axs[i,1].set_ylabel("Fraction number of points")
        axs[i,1].set_xlabel("Alleys")
        axs[i,1].set_xticks(np.arange(29))
        axs[i,1].set_xticklabels(a)
    fig.suptitle(title, y=1.04)
    fig.tight_layout()
    
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
    

def findRewards(df, pos, a=b"0x0002"):
    """
    Find which alley the animal was in and which direction it was going in 
    at the time of rewards
    """
    hrd, rec = util.readNEV(df)
    ts = [i[3] for i in rec if a in i[10]]
    rewardPos = [pos[bisect_left(pos[:,0], i)-1] for i in ts]
    rewardPos = np.column_stack((rewardPos, ts))
    inAlleys = np.empty(0)
    inAlleysDir = np.zeros((17, 4))
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
                if cross != 999: #999 = no border crossed
                    break
            if cross != 999:
                inAlleysDir[i,cross] += 1
    return inAlleys, inAlleysDir


def graphRewards(df, pos, title="", a=b"0x0002"):
    """
    Makes bar graphs of rewards in each alley and rewards in each alley divided by direction
    """
    inAlleys, inAlleysDir = findRewards(df, pos, a)
    
    fig, axs = plt.subplots(2,1)
    x = np.arange(17)
    width = 0.2
    axs[0].bar(x - 1.5*width, inAlleysDir[:,0], width, label="North")
    axs[0].bar(x - 0.5*width, inAlleysDir[:,1], width, label="East")
    axs[0].bar(x + 0.5*width, inAlleysDir[:,2], width, label="South")
    axs[0].bar(x + 1.5*width, inAlleysDir[:,3], width, label="West")
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
    

def directionTime(pos, posDir):
    """
    Finds the direction the animal was facing in each alley over time
    """
    rbins, cbins = [20, 20]
    rows = np.linspace(pos[0,0], pos[-1,0], rbins)
    hists = []
    axLimits = []
    #histArray = np.empty((0,19))
    for i in range(17):
        alley = alleyInterBounds[str(i)]
        aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
        if (aRect.xmax - aRect.xmin) > (aRect.ymax - aRect.ymin): #horizontal
            cols = np.linspace(aRect.xmin, aRect.xmax, cbins)
            axType = 1
            axLimit = [aRect.xmin, aRect.xmax]
        elif (aRect.xmax - aRect.xmin) < (aRect.ymax - aRect.ymin): #vertical
            cols = np.linspace(aRect.ymin, aRect.ymax, cbins)
            axType = 2
            axLimit = [aRect.ymin, aRect.ymax]
        
        inAlley1 = np.all((posDir[axType][:,1] > aRect.xmin, posDir[axType][:,1] < aRect.xmax,\
                           posDir[axType][:,2] > aRect.ymin, posDir[axType][:,2] < aRect.ymax),axis=0)
        hist1 = np.histogram2d(posDir[axType][inAlley1,0],posDir[axType][inAlley1,axType],bins=[rows,cols])[0]
        inAlley2 = np.all((posDir[(axType+2)%4][:,1] > aRect.xmin, posDir[(axType+2)%4][:,1] < aRect.xmax,\
                           posDir[(axType+2)%4][:,2] > aRect.ymin, posDir[(axType+2)%4][:,2] < aRect.ymax),axis=0)
        hist2 = np.histogram2d(posDir[(axType+2)%4][inAlley2,0],posDir[(axType+2)%4][inAlley2,axType],bins=[rows,cols])[0]
        
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


def graphDirTime(pos):
    """
    Graphs direction over time
    """
    posDir = directionFilter(pos)
    hists, axLimits = directionTime(pos, posDir)
    xlabels = np.round(np.linspace(pos[0,0]/1e6, pos[-1,0]/1e6, 8), 0).astype(int)
    fig = plt.figure(figsize=(12,18))
    for i in range(17):
        axs = ImageGrid(fig, (i%2*3/4,i/8-i%2/8,0.6,0.6), (1,2), axes_pad=0.1, cbar_mode="single")
        
        ylabels = np.round(np.linspace(axLimits[i][0], axLimits[i][1], 8), 0).astype(int)
        vmax = np.nanpercentile(np.hstack((hists[i][1],hists[i][2])),99)
        for j in range(2):
            im = axs[j].imshow(hists[i][j+1], origin="lower",vmax=vmax)
            axs[j].set_xlabel("Time (s)")
            axs[j].set_xticks(np.linspace(-0.5,18.5,8))
            axs[j].set_xticklabels(xlabels[:-1])
            axs[j].set_yticks(np.linspace(-0.5,18.5,8))
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
    return fig

def graphConfounds(dfData, dfGraph, title, a=b"0x0002"):
    """
    Puts all confound graphs in a multi-page pdf
    """
    with open(dfData+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = util.read_pos(dfData)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = adjustPosCamera(dfData, pos, ts)
    position = np.column_stack((ts, posx, posy))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    posFilt = Filt.velocity_filtering(position, 3)
    
    with PdfPages(dfGraph+title+" confounds") as pdf:
        plot = graphDirAndOcc(position, title+"\nNo velocity filtering")
        pdf.savefig(plot, bbox_inches="tight")
        plot = graphDirAndOcc(posFilt, title+"\nVelocity filtered with threshold = 3 cm/s")
        pdf.savefig(plot, bbox_inches="tight")
        plot = graphRewards(dfData, position, title, a)
        pdf.savefig(plot, bbox_inches="tight")
        plot = graphDirTime(posFilt)
        pdf.savefig(plot, bbox_inches="tight")
        
Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")