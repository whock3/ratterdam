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
    ax.bar(x - 0.75*width, visitRatesMeans[:,0], 0.5*width, label="North")
    ax.bar(x - 0.25*width, visitRatesMeans[:,1], 0.5*width, label="East")
    ax.bar(x + 0.25*width, visitRatesMeans[:,2], 0.5*width, label="South")
    ax.bar(x + 0.75*width, visitRatesMeans[:,3], 0.5*width, label="West")
    ax.set_ylabel("Average rate (Hz)")
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.set_xlabel("Subfields")
    ax.legend()
    fig.suptitle(title)
    
    return fig


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
    #ax.set_xlim(0, 170)
    #ax.legend(loc="lower right")
    #ax.set_title(title+f"\nvthresh = 3 cm/s\nCutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()
    
