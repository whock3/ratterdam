# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:41:47 2020

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from organized import makeRM2
from RDP4_class import RDP4
from turn3 import turn3
import numpy as np

def multipageTurnRM(pos, spikes, filename):
    """
    pos and spikes from shiftPosP and shiftPosS
    """
    vmax = 23.8
    titles = ["","Left", "Back", "Right", "North", "East", "South", "West"]
    rows, cols = [0,5,5,10,4,5,10,6], [0,8,8,12,6,6,10,10]
    figs = []
    
    for i in range(1,8):
        plt.rc('text', usetex=False)
        fig, axs = plt.subplots(rows[i], cols[i], figsize=(cols[i]*1.5,rows[i]*1.5), dpi=750)
        plt.suptitle(titles[i], y=0.99)
        for j, ax in enumerate(axs.flatten()):
            if j < len(pos[i]):
                n = makeRM2(spikes[i][j], pos[i][j])
                im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax, extent=(-7,7,-7,7))
                ax.axis("equal")
        cb = fig.colorbar(im)
        cb.set_label("Rate (Hz)")
        fig.tight_layout()
        figs.append(fig)
        plt.close(fig)
        
    with PdfPages(filename) as pdf:
        for fig in figs:
            pdf.savefig(fig)


def graph(pos, title=""): 
    """
    makes 10*10 graphs of 10 second segments of turn3
    pos: in cm before RDP
    """
    fig, axs = plt.subplots(10, 10, figsize=(30,30), dpi=500)
    b = 1000
    for ax in axs.flatten():
        RList=RDP4(pos[b:b+300],1).ResultList
        idx,_,_,_ = turn3(RList, np.pi/12, np.pi/4)
        ax.set_xlim(0, 120)
        ax.set_ylim(0, 100)
        
        ax.plot(RList[:,1], RList[:,2], zorder=1)
        ax.scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first", zorder=2)
        ax.scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last", zorder=2)
        ax.scatter(RList[idx,1], RList[idx,2], s = 16, color = "r", label='turning points', zorder=2)
        
        ax.set_title(title)
        ax.axis("equal")
        b += 300
    fig.tight_layout()