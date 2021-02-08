# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 17:05:23 2020

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import ast


def spikePlots(turns, unit, subfield, filename, df, title=""):
    """
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    df: path to the files with which locations a field is in
    """
    with open(df+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    fieldLoc = ast.literal_eval(data[subfield+1][1])
    turns2 = np.empty((0,5)) #turns in the field
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            turns2 = np.vstack((turns2, turn[:5].astype(float))) #allo before turn, ego, allo after turn, ts of exit, ts of entry
    
    fig, axs = plt.subplots(10,10,figsize=(20,20))
    axs = axs.flatten()
    for i in range(100):
        p = unit.position[np.logical_and(turns2[i,3]-0.5e6 < unit.position[:,0],
                                         unit.position[:,0] < turns2[i,4]+0.5e6)]
        s = unit.spikes[np.logical_and(turns2[i,3]-0.5e6 < unit.spikes[:,0],
                                         unit.spikes[:,0] < turns2[i,4]+0.5e6)]
        axs[i].plot(p[:,1], p[:,2], color=(0.5,0.5,0.5))
        axs[i].scatter(s[:,1], s[:,2], color="red", s=4, zorder=3)
        axs[i].scatter(p[0,1], p[0,2], c="green", s=4, zorder=4)
        axs[i].scatter(p[-1,1], p[-1,2], c="blue", s=4, zorder=4)
    fig.tight_layout()
    fig.suptitle(f"{title}\nGreen = beginning   Blue = end", y=1.02)