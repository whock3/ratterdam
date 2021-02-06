# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:36:49 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import ratterdam_Defaults as Def


def windowSize(pos, title=""):
    """
    Plots histograms: point to point distance
                      point to point distance and speed
    """
    pos = pos / np.array([1,Def.ptsCm,Def.ptsCm])
    dist = np.sqrt((pos[1:,1]-pos[:-1,1])**2 + (pos[1:,2]-pos[:-1,2])**2)
    speed = dist / (pos[1:,0]-pos[:-1,0]) * 1e6
    speed2 = speed[np.logical_and(dist < np.percentile(dist,95),
                                  speed < np.percentile(speed,95))]
    dist2 = dist[np.logical_and(dist < np.percentile(dist,95),
                                speed < np.percentile(speed,95))]
    fig, ax = plt.subplots()
    ax.hist(dist2, 30)#, np.arange(0,8,0.5))
    #ax.set_xticks(np.arange(0,39,3))
    ax.set_xlabel("Point to point distance (cm)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    
    fig, ax  = plt.subplots()
    hist = np.histogram2d(speed2, dist2, bins=[np.arange(0,52.5,2.5),np.arange(0,1.8,0.05)])[0]
    vmax = np.percentile(hist,99)
    im = ax.imshow(hist, origin="lower", vmin=0, vmax=vmax)
    print(max(speed2), max(dist2), np.mean(dist2), np.mean(dist))
    ax.set_xticks(np.arange(-0.5,34.5,4))
    ax.set_xticklabels(np.arange(0,1.8,0.2).round(1))
    ax.set_yticks(np.arange(-0.5,21.5,2))
    ax.set_yticklabels(np.arange(0,55,5))
    ax.set_xlabel("Point to point distance (cm)")
    ax.set_ylabel("Speed (cm/s)")
    ax.set_title(title+f"\n Cutoff = 99th percentile, {int(vmax)}")
    cb = fig.colorbar(im)
    cb.set_label("Count")