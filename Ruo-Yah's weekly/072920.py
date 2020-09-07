# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 21:14:19 2020

@author: Ruo-Yah Lai
"""
from williamDefaults import alleyInterBounds
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import datetime
from string import ascii_uppercase


#modified from 060220
def rect(bounds, color):
    bottom_left = [bounds[0][0], bounds[1][0]]
    top_right = [bounds[0][1], bounds[1][1]]
    return patches.Rectangle(tuple(bottom_left), (top_right[0]-bottom_left[0]), 
                             (top_right[1]-bottom_left[1]), color=color, alpha=0.4)

def graph(pos):
    """
    Graphs the alleys and intersections
    """
    fix,ax = plt.subplots(1)
    for i in range(17):
        ax.add_patch(rect(alleyInterBounds[str(i)], "r"))
        ax.annotate(i, (alleyInterBounds[str(i)][0][0], alleyInterBounds[str(i)][1][0]),
                    xytext=(alleyInterBounds[str(i)][0][0], alleyInterBounds[str(i)][1][0]))
    for i in ascii_uppercase[:12]:
        ax.add_patch(rect(alleyInterBounds[i], "g"))
        ax.annotate(i, (alleyInterBounds[i][0][0], alleyInterBounds[i][1][0]),
                    xytext=(alleyInterBounds[i][0][0], alleyInterBounds[i][1][0]))
    ax.plot(pos[:,1], pos[:,2])
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 450)
    ax.axis("equal")