# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 17:18:08 2020

@author: Ruo-Yah Lai
"""
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import ratterdam_ParseBehavior as Parse
import ratterdam_RepetitionCoreFx as RepCore
import ratterdam_Defaults as Def
import numpy as np
import utility_fx as util
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt
import matplotlib.pyplot as plt
from bisect import bisect_left


def loadRepeatingUnit2(df, clustName, smoothing=2, vthresh=Def.velocity_filter_thresh):
    """
    Same as loadRepeatingUnit but keeps head direction in pos
    
    take a path to a data dir
    load spikes and position into two np arrays
    spikes is (n,1) and pos is (4,n) cols of ts,x,y,dir
    use cameraOrientationInfo.txt to flip axes if needed
    use sessionEpochInfo.txt, specific for open Ratterdam exp
    to get session ts and clip spikes/pos"""
    
    with open(df+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = util.read_pos(df)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = Parse.adjustPosCamera(df, pos, ts)
    headDir = [pos[i][2] for i in ts]
    position = np.column_stack((ts, posx, posy, headDir))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    position = position[np.logical_or(position[:,1]>0, position[:,2]>0)]
    position = Filt.velocity_filtering(position, vthresh)
    clust = np.asarray(util.read_clust(df+clustName))
    clust = Filt.unitVelocityFilter(ts, position, clust)
    clust = clust[(clust >= start) & (clust <= end)]
    spikexy = util.getPosFromTs(clust,position[:,:3])
    spikes = np.column_stack((clust,spikexy))
    
    unit = RepCore.Unit(spikes,position[:,:3], clustName, smoothing)
    unit = Filt.filterFields(unit)
    if smoothing:
        unit.smoothFieldsFx()
        
    unit.position = position
    return unit


def plotVisitPoints(unit, subfield):
    """
    Scatter plot all of the points in visits to a subfield
    """
    points = []
    fig, ax = plt.subplots()
    for visit in unit.visits[subfield]:
        for ts in visit:
            index = np.where(unit.position[:,0]==ts)[0]
            ax.scatter(unit.position[index,1], unit.position[index,2])
            

#from 072220
def thetaSpeed(pos):
    """
    Calculates theta/s and speed between pos points
    """
    posFilt = pos[np.where(pos[:,3] != -99)]
    directions = np.diff(posFilt[:,1:3], axis=0)
    dist = np.linalg.norm(directions, axis=1)[1:] / np.diff(posFilt[:,0])[1:] *1e6 / Def.ptsCm
    
    theta = np.diff(posFilt[1:,3])
    negatives = np.where(theta < -180)
    positives = np.where(theta > 180)
    theta[negatives] = theta[negatives] + 360
    theta[positives] = theta[positives] - 360 #[-180, 180]
    theta = theta / np.diff(posFilt[:,0])[1:] *1e6
    
    return np.column_stack((posFilt[1:-1,0], theta, dist))


def thetaSpeedNearField(unit, subfield):
    """
    Finds ∆θ/s and ∆S/s near a subfield
    position: 4 columns (ts,x,y,headDir)
    Returns [ts, theta/s, speed]
    """
    ranges = np.empty((0,2))
    for i in range(len(unit.visits[subfield])): #the ith visit
        #start and end are indexes of the start and end of visit
        start = bisect_left(unit.position[:,0], unit.visits[subfield][i][0]-400000)
        end = bisect_left(unit.position[:,0], unit.visits[subfield][i][-1]+400000)
        if end-start > 0:
            ranges = np.vstack((ranges, np.reshape([start, end],(1,2))))
            
    postd = np.empty((0, 3))
    for i in range(len(ranges)): #the ith visit
        postd1 = thetaSpeed(unit.position[int(ranges[i,0]):int(ranges[i,1])])
        postd = np.vstack((postd, postd1))
    return postd


def makeRM(spikes, pos, rmax, cmax, smoothing_2d_sigma=2, bins = [50,50]):
    """
    For speed vs change in theta
    """
    #rmax = np.percentile(pos[:,2], clip)
    #cmax = max(abs(np.percentile(pos[:,1], 100-clip)), np.percentile(pos[:,1], clip))
    rows = np.linspace(0, rmax, bins[0])
    cols = np.linspace(-cmax, cmax, bins[1])
    hs,xe,ye = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rows, cols])
    ho = np.histogram2d(pos[:,2],pos[:,1],bins=[rows, cols])[0]
    n = (hs*np.reciprocal(ho))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(ho==0)] = np.nan
    n = util.weird_smooth(n,smoothing_2d_sigma)
    n[np.where(ho==0)] = np.nan
    return n, (int(cmax), int(rmax))


def graphRMSpeedTheta(unit, title, subfieldDim, clip=99, percentile=99, sigma=2):
    """
    unit: from loadRepeatingUnit2 (position includes head direction)
    subfieldDim: list, dimensions (rows,columns) for putting all subfields on the same plot
    clip: percentile for clipping ∆θ/sec and ∆S/sec
    percentile: percentile for setting the max of the colormap
    
    Graphs speed vs change in theta near subfields
    """
    fig, axs = plt.subplots(subfieldDim[0], subfieldDim[1], figsize=(subfieldDim[1]*4,subfieldDim[0]*4))
    thetas = []
    speeds = []
    postds = []
    for i in range(len(unit.perimeters)):
        postd = thetaSpeedNearField(unit, i)
        thetas = np.hstack((thetas, postd[:,1]))
        speeds = np.hstack((speeds, postd[:,2]))
        postds.append(postd)
    rmax = np.percentile(speeds, clip)
    cmax = max(abs(np.percentile(thetas, 100-clip)), np.percentile(thetas, clip))
    
    fig.suptitle(title + f"   vthresh = 3 cm/s\nClipped at {clip}th percentile of ∆θ/sec and ∆S/sec", y=1.05)
    axs = axs.flatten()
    for i in range(len(unit.perimeters)):
        spikets = Filt.unitVelocityFilter(unit.position[:,0], postds[i], unit.spikes[:,0])
        spikexy = util.getPosFromTs(spikets,postds[i])
        spikes = np.column_stack((spikets,spikexy))
        n, minmax = makeRM(spikes, postds[i], rmax, cmax, sigma)
        vmax = np.nanpercentile(n, percentile)
        im = axs[i].imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)#, extent=minmax)
        axs[i].set_xlabel("∆θ/sec (deg/sec)")
        axs[i].set_ylabel("∆S/sec (cm/sec)")
        #ax.set_xlim(minmax[0], minmax[1])
        #ax.set_ylim(minmax[2], minmax[3])
        axs[i].set_xticks([0, 12.5, 25, 37.5, 50])
        axs[i].set_xticklabels([-minmax[0], int(-minmax[0]/2), 0, int(minmax[0]/2), minmax[0]])
        axs[i].set_yticks([0, 12.5, 25, 37.5, 50])
        axs[i].set_yticklabels([0, int(minmax[1]/4), int(minmax[1]/2), int(minmax[1]/4*3), minmax[1]])
        axs[i].set_title(f"Subfield {i}\nCutoff = {percentile}th %tile, {round(vmax,1)} Hz", c=unit.colors[i])
        cb = fig.colorbar(im, ax=axs[i])
        cb.set_label("Rate (Hz)")
    fig.tight_layout()
    
    for i in range(len(unit.perimeters), subfieldDim[0]*subfieldDim[1]):
        axs[i].axis("off")


def graphRM(unit, title="", percentile=98):
    fig, ax = plt.subplots(figsize=(7.5,5))
    n = util.makeRM(unit.spikes, unit.position)
    vmax = np.nanpercentile(n, percentile)
    im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/Def.ptsCm, unit.perimeters[i][:,1]/Def.ptsCm, 
                color=unit.colors[i], label=f"Subfield {i}")
    ax.set_xlim(0, 170)
    ax.legend(loc="lower right")
    ax.set_title(title+f"\nvthresh = 3 cm/s\nCutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()
    
    
def headDirOrientation(df):
    with open(df+"cameraOrientationInfo.txt","r") as cam:
         info = cam.readlines()
    info = [i.rstrip() for i in info]
    info = {info[0][0]:info[0][2:], info[1][0]:info[1][2:]}
    if (info['x'] == 'None' and info['y'] == 'None') or (info['x'] != 'None' and info['y'] != 'None'):
        return "cw"
    elif (info['x'] == 'None' and info['y'] != 'None') or (info['x'] != 'None' and info['y'] == 'None'):
        return "ccw"

def bulkGraphs(file, title, timestamp, df, df2="W:/Ratterdam/R859/ratterdam_tabulations/", dfGraph="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/Speed vs theta change/"):
    """
    file: file name of the file with repetition tabulations
    title: part of the title of graphs (include rat and day)
    df: path to data
    df2: path to repetition tabulations
    dfGraph: path to where to save the graphs
    
    Make speed vs change in head direction rate maps and regular rate maps
    """
    tooManyFields = []
    cells = util.readRepeatingCells(file, df2)
    for cell in cells:
        cellTitle = cell.replace("\\cl-maze", " ")
        unit = loadRepeatingUnit2(df, cell)
        graphRM(unit, title+" "+cellTitle)
        plt.savefig(dfGraph + f"{timestamp} - {title} {cellTitle}.png",
                    bbox_inches="tight")
        plt.close()
        if len(unit.perimeters) <= 2:
            subfieldDim = [2,1]
        elif 2 < len(unit.perimeters) <= 4:
            subfieldDim = [2,2]
        elif 4 < len(unit.perimeters) <= 6:
            subfieldDim = [2,3]
        elif 6 < len(unit.perimeters) <= 9:
            subfieldDim = [3,3]
        elif 9 < len(unit.perimeters) <= 12:
            subfieldDim = [3,4]
        elif 12 < len(unit.perimeters) <= 16:
            subfieldDim = [4,4]
        else:
            tooManyFields.append([cell, len(unit.perimeters)])
            continue
        graphRMSpeedTheta(unit, f"{title} {cellTitle}\n+∆θ = {headDirOrientation(df)}", subfieldDim, 98)
        plt.savefig(dfGraph + f"{timestamp} - {title} {cellTitle}.png",
                    bbox_inches="tight")
        plt.close()
    print(tooManyFields)