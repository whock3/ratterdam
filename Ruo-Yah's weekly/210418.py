# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 12:05:03 2021

@author: Ruo-Yah Lai
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as path
from confounds import directionFilterS
import utility_fx as util
import ratterdam_Defaults as Def
import williamDefaults as wmDef
from scipy.stats import binom_test


def graphFieldDir(unit, title):
    """
    Bar graph: 4 bars (for the 4 directions) for each subfield
    """
    visitRates = []
    for i in range(len(unit.perimeters)): #i = subfield
        visitRates.append([[],[],[],[]])
        contour = path.Path(unit.perimeters[i])
        PinC = unit.position[contour.contains_points(unit.position[:,1:])]
        SinC = unit.spikes[contour.contains_points(unit.spikes[:,1:])]
        for visit in range(len(unit.visits[i])):
            
            visitPos = np.logical_and(PinC[:,0] >= unit.visits[i][visit][0],
                                      PinC[:,0] <= unit.visits[i][visit][-1])
            visitSpk = np.logical_and(SinC[:,0] >= unit.visits[i][visit][0],
                                      SinC[:,0] <= unit.visits[i][visit][-1])
            
            posDir, spikesDir = directionFilterS(PinC[visitPos], SinC[visitSpk])
            for j in range(4):
                if len(posDir[j]) == 0:
                    pass
                else:
                    visitRate = len(spikesDir[j]) / len(posDir[j]) * 30
                    visitRates[i][j].append(visitRate)
    
    visitRatesMeans = np.empty((len(unit.perimeters),4))
    for i in range(len(unit.perimeters)):
        for j in range(4):
            visitRatesMeans[i,j] = np.mean(visitRates[i][j])
    
    fig, ax = plt.subplots()
    x = np.arange(len(unit.perimeters))
    width = 0.4
    colors = ["#b2e2e2", "#66c2a4", "#238b45", "#006d2c"]
    ax.bar(x - 0.75*width, visitRatesMeans[:,0], 0.5*width, label="North", color=colors[0])
    ax.bar(x - 0.25*width, visitRatesMeans[:,1], 0.5*width, label="East", color=colors[1])
    ax.bar(x + 0.25*width, visitRatesMeans[:,2], 0.5*width, label="South", color=colors[2])
    ax.bar(x + 0.75*width, visitRatesMeans[:,3], 0.5*width, label="West", color=colors[3])
    ax.set_ylabel("Average rate (Hz)")
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.set_xlabel("Subfields")
    ax.legend()
    fig.suptitle(title)
    
    return fig


def graphRM(unit, title="", percentile=98):
    """
    Rate map with smaller pixel size (1 cm^2)
    """
    fig, ax = plt.subplots(figsize=(7.5,5))
    n = util.makeRM(unit.spikes, unit.position)
    vmax = np.nanpercentile(n, percentile)
    im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/Def.ptsCm, unit.perimeters[i][:,1]/Def.ptsCm, 
                color=unit.colors[i], label=f"Subfield {i}")
    #ax.set_xlim(0, 170)
    #ax.legend(loc="lower right")
    ax.set_title(title+f"\nvthresh = 3 cm/s\nCutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()


def graphRM2(unit, title=""):
    """
    Rate map with the same pixel size as RateMapClass
    """
    fig, ax = plt.subplots(figsize=(7.5,5))
    vmax = np.nanpercentile(unit.repUnit.rateMap2D, 98)
    im = ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                   cmap="jet", vmax=vmax, extent=[wmDef.xedges[0], wmDef.xedges[-1]/Def.ptsCm, 
                                                  wmDef.yedges[0], wmDef.yedges[-1]/Def.ptsCm])
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/Def.ptsCm, unit.perimeters[i][:,1]/Def.ptsCm,
                color=unit.colors[i], label=f"Subfield {i}")
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    ax.axis("equal")
    ax.set_title(title+f"\nvelocity threshold = 3 cm/s\nCutoff = 98th percentile, {round(vmax,1)} Hz")
    fig.tight_layout()
    return fig
    
    
def graphStackedBars():
    labels = ["R781 D2", "R781 D3", "R781 D4", "R808 D6", "R808 D7", "R859 D1",
              "R859 D2", "R859 D3", "R886 D1", "R886 D2", "R886 D3"]
    nonRep = [12,18,9,0,2,7,14,7,6,2,2]
    repNoInter = [1,2,2,6,1,8,3,12,1,2,0]
    repInter = [2,0,0,1,0,11,23,8,1,2,0]
    width = 0.35       # the width of the bars: can also be len(x) sequence
    x = np.arange(len(labels))
    
    fig, ax = plt.subplots()
    colors = ["#fecc5c", "#fd8d3c", "#e31a1c"]
    ax.bar(x, repInter, width, label='Repeating interaction', color=colors[0])
    ax.bar(x, repNoInter, width, bottom=repInter,
           label='Repeating no interaction', color=colors[1])
    ax.bar(x, nonRep, width, bottom=np.array(repInter)+np.array(repNoInter),
           label='Non-repeating', color=colors[2])
    
    ax.set_ylabel('Number of cells')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    for tick in ax.get_xticklabels():
        tick.set_rotation(55)
    #ax.set_title('Scores by group and gender')
    ax.legend()
    
    for i in range(len(labels)):
        if binom_test(repInter[i],repInter[i]+repNoInter[i],0.05) < 0.05:
            plt.annotate("*", (i-0.12,nonRep[i]+repNoInter[i]+repInter[i]))
    plt.tight_layout()
    return fig