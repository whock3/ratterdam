# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 12:46:45 2020

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import utility_fx as util


def graphRM(unit, title="", percentile=98):
    fig, ax = plt.subplots()
    n = util.makeRM(unit.spikes, unit.position)
    vmax = np.nanpercentile(n, percentile)
    im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    colors = ["b","g","r","k","c","m","y","orange","b","g","r","k","c","m","y","orange"]
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/4.72, unit.perimeters[i][:,1]/4.72, color=colors[i], 
                label=f"Subfield {i}")
    ax.legend()
    ax.set_title(title+f"\nCutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()