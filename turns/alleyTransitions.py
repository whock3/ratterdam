# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:13:25 2020

@author: Ruo-Yah Lai
"""

from williamDefaults import alleyInterBounds, alleyLines
from string import ascii_uppercase
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy


def getVisits(data, maxgap=0.5e6, minTime=0.5e6):
    """
    Leaves the field for <= maxgap and is in the field for > minTime
    """
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    for group in groups:
        if abs(group[0] - group[-1]) < minTime:
            groups.remove(group)
    return groups


#for alleys then intersections:
    #find point that are inside
    #getVisits
    #all points between 1st and last pt of a visit is in that location
    #put into a new array with [ts,x,y,location,whether it's the last pt of the visit]
    #remove these points from the overall list
#sort the array by ts
def alleyTransitions(pos):
    letters = [i for i in ascii_uppercase[:12]]
    a = np.hstack((np.arange(17).astype(int), letters))
    posList = deepcopy(pos)
    posNew = np.empty((0,5))
    for i in a:
        alley = alleyInterBounds[i]
        inAlley = np.all((alley[0][0] < posList[:,1], posList[:,1] < alley[0][1],
                          alley[1][0] < posList[:,2], posList[:,2] < alley[1][1]), axis=0)
        
        if len(posList[inAlley]) == 0:
            continue
        
        visits = getVisits(posList[inAlley,0], 0.5e6)
        #print(i, len(visits[-1]))
        for visit in visits:
            idxPos = np.all((visit[0]<=posList[:,0], posList[:,0]<=visit[-1]), axis=0) #bool

            idxStart = np.where(idxPos)[0][0]
            idxEnd = np.where(idxPos)[0][-1]+1
            posnew = np.column_stack((posList[idxPos], np.full(idxEnd-idxStart,i), np.zeros(idxEnd-idxStart)))
            posnew[-1,4] = 1 #last point of this visit
            posNew = np.vstack((posNew, posnew))
            
            posList = list(posList)
            del(posList[idxStart:idxEnd])
            posList = np.array(posList)
        
        if len(posList) == 0:
            break
    posNew = posNew[posNew[:,0].argsort()]
    
    
    fig, ax = plt.subplots()
    ax.plot(posNew[:,1].astype(float), posNew[:,2].astype(float))
    ax.scatter(float(posNew[0,1]), float(posNew[0,2]), marker = "+", color = "r", label = "first", zorder=2)
    ax.scatter(float(posNew[-1,1]), float(posNew[-1,2]), marker = "x", color = "r", label = "last", zorder=2)
    
    for i in range(14):
        plt.plot(alleyLines[i,:,0], alleyLines[i,:,1], color="k", alpha=0.5)
    
    for i in np.where(posNew[:,4].astype(float) == 1)[0]:
        #print(posNew[i])
        ax.annotate(posNew[i,3], (float(posNew[i,1])+2, float(posNew[i,2])+2))
    ax.axis("equal")
    ax.legend()
    return posNew