# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 16:48:12 2021

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import utility_fx as util
import ratterdam_Defaults as Def
from williamDefaults import binWidth
from newAlleyBounds import R859


def perim(PFmask):
    r, c = PFmask.shape
    expanded = np.zeros((r*2,c*2), dtype=PFmask.dtype)
    expanded[::2, ::2] = PFmask
    print(expanded)
    # filling in using values to the left and right
    for i in range(0,r*2,2):
        for j in range(1,c*2-1,2):
            expanded[i,j] = expanded[i,j-1] and expanded[i,j+1]
    print(expanded)
    # filling in using values above and below
    for i in range(1,r*2-1,2):
        for j in range(c*2):
            expanded[i,j] = expanded[i-1,j] and expanded[i+1,j]
            
    expanded = np.vstack((np.zeros((1,expanded.shape[1])), expanded))
    expanded = np.hstack((np.zeros((expanded.shape[0],1)), expanded))
    return expanded


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
    

def borderComparison(unit, subfield, a, title=""):
    fig, axs = plt.subplots(1,2)
    axs[0].scatter(unit.repUnit.PF[subfield].perimeter[1]*binWidth/a, unit.repUnit.PF[subfield].perimeter[0]*binWidth/a)
    axs[0].axis("equal")
    axs[1].plot(unit.perimeters[subfield][:,0], unit.perimeters[subfield][:,1])
    axs[1].axis("equal")
    fig.suptitle(title, y=1.04)
    fig.tight_layout()
    

