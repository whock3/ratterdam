# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 21:14:22 2020

@author: Ruo-Yah Lai
"""
import numpy as np
from williamDefaults import binWidth
from copy import deepcopy


def reorderBorder(border):
    """
    Reorder the points that make up the border of subfields
    """
    border = list(zip(border[1], border[0])) #[x,y]
    border2 = deepcopy(border)
    firstPoint = border[0]
    point = border[0]
    newBorder = np.array(point)
    while True:
        if (point[0]+1, point[1]) in border: #to the right
            border.remove(point)
            point = (point[0]+1, point[1])
            newBorder = np.vstack((newBorder, point))
        elif (point[0], point[1]+1) in border: #up
            border.remove(point)
            point = (point[0], point[1]+1)
            newBorder = np.vstack((newBorder, point))
        elif (point[0]-1, point[1]) in border: #left
            border.remove(point)
            point = (point[0]-1, point[1])
            newBorder = np.vstack((newBorder, point))
        elif (point[0], point[1]-1) in border: #down
            border.remove(point)
            point = (point[0], point[1]-1)
            newBorder = np.vstack((newBorder, point))
        
        elif (point[0], point[1]-1) in border2: #down
            border.remove(point)
            border2.remove(point)
            point = (point[0], point[1]-1)
            newBorder = np.vstack((newBorder, point))
            border.append(point)
        elif (point[0]+1, point[1]) in border2: #right
            border.remove(point)
            border2.remove(point)
            point = (point[0]+1, point[1])
            newBorder = np.vstack((newBorder, point))
            border.append(point)
        elif (point[0], point[1]+1) in border2: #up
            border.remove(point)
            border2.remove(point)
            point = (point[0], point[1]+1)
            newBorder = np.vstack((newBorder, point))
            border.append(point)
        elif (point[0]-1, point[1]) in border2: #left
            border.remove(point)
            border2.remove(point)
            point = (point[0]-1, point[1])
            newBorder = np.vstack((newBorder, point))
            border.append(point)
        
        if point == firstPoint:
            break
        
    return newBorder*binWidth+binWidth/2


def reorderBorders(unit):
    borders = []
    for i in range(len(unit.repUnit.PF)):
        borders.append(reorderBorder(unit.repUnit.PF[i].perimeter))
    return borders

"""
newBorder = reorderBorder(unit1.repUnit.PF[0].perimeter)
fig, axs = plt.subplots(1,3)
axs[0].scatter(unit1.repUnit.PF[0].perimeter[1], unit1.repUnit.PF[0].perimeter[0])
axs[0].axis("equal")
axs[1].plot(newBorder[:,0], newBorder[:,1])
axs[1].axis("equal")
axs[2].plot(unit1.perimeters[0][:,0], unit1.perimeters[0][:,1])
axs[2].axis("equal")
fig.suptitle("R859 D3 T6 1.1 subfield 0", y=1.04)
fig.tight_layout()
"""