# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:45:04 2021

@author: Ruo-Yah Lai
"""
from RDP4_class import RDP4
from turn3 import turn3
import ratterdam_Defaults as Def
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from newAlleyBounds import R859
from string import ascii_uppercase


#Rpos = RDP4(unit.position, Def.ptsCm).ResultList
def graphTurns(Rpos):
    """
    Graphs turns as labeled by turn3
    """
    a = np.array([1, Def.ptsCm, Def.ptsCm])
    Rpos = Rpos/a[None,:]
    
    fig, ax = plt.subplots()
    ax.scatter(Rpos[0,1], Rpos[0,2], marker = "+", color = "r", label = "first", zorder=3)
    ax.scatter(Rpos[-1,1], Rpos[-1,2], marker = "x", color = "r", label = "last", zorder=3)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    ax.plot(Rpos[:,1], Rpos[:,2])
    ax.axis('equal')
    
    idx3,_,_,_ = turn3(Rpos, np.pi/12, np.pi/4, 5, 15)
    
    ax.scatter(Rpos[idx3,1], Rpos[idx3,2], s = 16, color = "r", label='turns', zorder=3)
    plt.legend()
    

def rect(bounds, color):
    bottom_left = [bounds[0][0], bounds[1][0]]
    top_right = [bounds[0][1], bounds[1][1]]
    return patches.Rectangle(tuple(bottom_left), (top_right[0]-bottom_left[0]), 
                             (top_right[1]-bottom_left[1]), color=color, alpha=0.3)


def graph(alleyInterBounds=R859.alleyInterBounds):
    """
    Graphs the alleys and intersections
    """
    fix,ax = plt.subplots(1)
    for i in range(17):
        ax.add_patch(rect(alleyInterBounds[str(i)], "r"))
        x = (alleyInterBounds[str(i)][0][0] + alleyInterBounds[str(i)][0][1])/2-5
        y = (alleyInterBounds[str(i)][1][0] + alleyInterBounds[str(i)][1][1])/2-5
        ax.annotate(i, (x, y), xytext=(x, y))
    for i in ascii_uppercase[:12]:
        ax.add_patch(rect(alleyInterBounds[i], "g"))
        x = (alleyInterBounds[str(i)][0][0] + alleyInterBounds[str(i)][0][1])/2-5
        y = (alleyInterBounds[str(i)][1][0] + alleyInterBounds[str(i)][1][1])/2-5
        ax.annotate(i, (x, y), xytext=(x, y))
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 450)
    ax.axis("equal")
