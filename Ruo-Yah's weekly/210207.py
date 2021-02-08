# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 19:59:22 2021

@author: Ruo-Yah Lai
"""

import numpy as np
import matplotlib.path as path
import matplotlib.pyplot as plt
import scipy
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import utility_fx as util
import ratterdam_Defaults as Def
from confounds import directionFilterS
from ratterdam_RepetitionCoreFx import loadRepeatingUnit


def a(unit, array, title=""):
    """
    Plot the location of local peaks in firing rate
    """
    peaks = np.array(list(zip(array[1], array[0])))*10
    
    fig, ax = plt.subplots(figsize=(7.5,5))
    n = util.makeRM(unit.spikes, unit.position)
    vmax = np.nanpercentile(n, 98)
    im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/Def.ptsCm, unit.perimeters[i][:,1]/Def.ptsCm, 
                color=unit.colors[i], label=f"Subfield {i}")
    for peak in peaks:
        ax.scatter(peak[0]/Def.ptsCm, peak[1]/Def.ptsCm, color="black", zorder=5, s=4)
    ax.set_xlim(0, 170)
    ax.legend(loc="lower right")
    ax.set_title(title+f"\nvthresh = 3 cm/s\nCutoff = 98th percentile, {round(vmax,1)} Hz")
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()
    

def rateDir(unit, title):
    """
    Plots a bar graph of the average rate vs NESW for each subfield
    """
    posDir, spikesDir = directionFilterS(unit.position, unit.spikes)
    avgRate = np.empty((len(unit.perimeters), 4)) #rows=subfields, columns=NESW
    rm = []
    for i in range(len(unit.perimeters)):
        perim = unit.perimeters[i]
        contour = path.Path(perim)
        rm.append([])
        for j in range(4):
            spkIn = spikesDir[j][contour.contains_points(spikesDir[j][:,1:])]
            occIn = posDir[j][contour.contains_points(posDir[j][:,1:])]
            rm[-1].append(util.makeRM(spkIn,occIn))
            avgRate[i,j] = np.nanmean(rm[-1][-1])
    
    if len(unit.perimeters) <= 2:
        subfieldDim = [2,1]
    elif 2 < len(unit.perimeters) <= 4:
        subfieldDim = [2,2]
    elif 4 < len(unit.perimeters) <= 6:
        subfieldDim = [2,3]
    elif 6 < len(unit.perimeters) <= 9:
        subfieldDim = [3,3]
    fig, axs = plt.subplots(subfieldDim[0], subfieldDim[1])
    axs = axs.flatten()
    sig = ""
    for i in range(len(axs)):
        rmNoNan = []
        for j in range(4):
            rmNoNan.append(rm[i][j].flatten()[~np.isnan(rm[i][j].flatten())])
        print(len(rmNoNan[2]))
        f, p = scipy.stats.f_oneway(rmNoNan[0], rmNoNan[1], rmNoNan[2], rmNoNan[3])
        print(f, p)
        if p < 0.05:
            sig = "*"
        axs[i].bar(np.arange(4), avgRate[i])
        axs[i].set_title(f"Subfield {i} {sig}")
        axs[i].set_ylabel("Average firing rate (Hz)")
        axs[i].set_xticks(np.arange(4))
        axs[i].set_xticklabels(["N", "E", "S", "W"])
    fig.suptitle(title, y=1.04)
    fig.tight_layout()
    
    
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
    

def bulkGraphs(file, title, timestamp, df, title2="", df2="W:/Ratterdam/R859/ratterdam_tabulations/", dfGraph="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"):
    """
    file: file name of the file with repetition tabulations
    title: part of the title of graphs (include rat and day)
    df: path to data
    df2: path to repetition tabulations
    dfGraph: path to where to save the graphs
    
    Make speed vs change in head direction rate maps and regular rate maps
    """
    with open(df2+file,"r") as f:
        lines = f.readlines()
        tabulations = [line.split(",") for line in lines]
    tabulations = np.array(tabulations)
    cells = tabulations[1:, 0]
    for cell in cells:
        cellTitle = cell.replace("\\cl-maze", " ")
        unit = loadRepeatingUnit(df, cell)
        graphRM(unit, title+" "+cellTitle+title2)
        plt.savefig(dfGraph + f"{timestamp} - {title} {cellTitle} rateThresh 0.2.png",
                    bbox_inches="tight")
        plt.close()