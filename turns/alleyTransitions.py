# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:13:25 2020

@author: Ruo-Yah Lai
"""

from williamDefaults import alleyInterBounds, alleyLines
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from collections import namedtuple
import csv


def getVisits(data, maxgap=0.5e6):
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
    return groups


def ccw(a, b, c):
    """
    From https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])

def intersect(A,B,C,D):
    """
    From https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def crossBorder(pos1, pos2, alley):
    if intersect((alley.xmin,alley.ymax), (alley.xmax,alley.ymax), pos1[1:3], pos2[1:3]):
        border = 0 #North
    elif intersect((alley.xmax,alley.ymin), (alley.xmax,alley.ymax), pos1[1:3], pos2[1:3]):
        border = 1 #East
    elif intersect((alley.xmin,alley.ymin), (alley.xmax,alley.ymin), pos1[1:3], pos2[1:3]):
        border = 2 #South
    elif intersect((alley.xmin,alley.ymin), (alley.xmin,alley.ymax), pos1[1:3], pos2[1:3]):
        border = 3 #West
    else:
        border = 999 #not a turn
    return border

def turn(a, b):
    """
    First border crossed is a. Second border crossed is b.
    """
    if abs(a-b) == 2:
        egoturn = "s"
    elif a-b == 1 or a-b == -3:
        egoturn = "r"
    elif a-b == -1 or a-b == 3:
        egoturn = "l"
    elif a == b:
        egoturn = "b"
    return egoturn


def egoAllo(ego, allo):
    """
    Returns the allocentric direction that the animal was facing before turning
    ego: 1=S, 2=R, 3=B, 4=L
    allo(after turning): 1=N, 2=E, 3=S, 4=W
    """
    turn = (allo - (ego-1)) % 4
    if turn == 0:
        return 4
    else:
        return turn
    
    
#for alleys then intersections:
    #find point that are inside
    #getVisits
    #all points between 1st and last pt of a visit is in that location
    #put into a new array with [ts,x,y,location,whether it's the last pt of the visit]
    #remove these points from the overall list
#sort the array by ts
def alleyTransitions(pos, minTime=0.5e6):
    posList = deepcopy(pos)
    posNew = np.empty((0,5)) #ts, x, y, alley, 1st or last pt of visit
    alleys = []
    for i in np.arange(17):
        alley = alleyInterBounds[str(i)]
        alley = Rectangle(alley[0][0], alley[0][1], alley[1][0], alley[1][1])
        alleys.append(alley)
        inAlley = np.all((alley.xmin < posList[:,1], posList[:,1] < alley.xmax,
                          alley.ymin < posList[:,2], posList[:,2] < alley.ymax), axis=0)
        
        if len(posList[inAlley]) == 0:
            continue
        
        visits = getVisits(posList[inAlley, 0], 0.5e6)
        for visit in visits:            
            idxPos = np.all((visit[0]<=posList[:,0], posList[:,0]<=visit[-1]), axis=0) #bool
                        
            visitTime = abs(visit[0] - visit[-1])
            xDist = (max(posList[idxPos,1])-min(posList[idxPos,1])) / (alley.xmax-alley.xmin)
            yDist = (max(posList[idxPos,2])-min(posList[idxPos,2])) / (alley.ymax-alley.ymin)
            
            if visitTime > minTime or xDist > 0.7 or yDist > 0.7:            
                idxStart = np.where(idxPos)[0][0]
                idxEnd = np.where(idxPos)[0][-1]+1
                posnew = np.column_stack((posList[idxPos], np.full(idxEnd-idxStart,i), np.zeros(idxEnd-idxStart)))
                posnew[0,4] = 1 #1st point of this visit
                posnew[-1,4] = 2 #last point of this visit
                posNew = np.vstack((posNew, posnew))
            
                posList = list(posList)
                del(posList[idxStart:idxEnd])
                posList = np.array(posList)
        
        if len(posList) == 0:
            break
    posNew = posNew[posNew[:,0].argsort()]
    
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Alley transitions egocentric.csv","r") as file:
        data_iter = csv.reader(file)
        trsnEgo = [data for data in data_iter]
        trsnEgo = np.array(trsnEgo)
        trsnEgo = trsnEgo.astype(int)
    
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Alley transitions allocentric.csv","r") as file:
        data_iter = csv.reader(file)
        trsnAllo = [data for data in data_iter]
        trsnAllo = np.array(trsnAllo)
        trsnAllo = trsnAllo.astype(int)
    
    turns = np.empty((0,3)) #allo before turn, ego, allo after turn
    entries = np.where(posNew[:,4] == 1)[0]
    for Exit in np.where(posNew[:,4] == 2)[0]: #for each alley exit
        if Exit+1 == len(posNew):
            break
        entry = min(entries[entries > Exit])
        
        ego = trsnEgo[int(posNew[Exit,3]), int(posNew[entry,3])]
        if ego == 3: #back around turn
            allo2 = crossBorder(posNew[Exit], posNew[Exit+1], alleys[int(posNew[Exit,3])])+1 #0 index to 1 index
            if allo2 == 1000:
                break
        else:
            allo2 = trsnAllo[int(posNew[Exit,3]), int(posNew[entry,3])]
        allo1 = egoAllo(ego, allo2)
        turns = np.vstack((turns, (allo1, ego, allo2)))
    
    fig, ax = plt.subplots()
    ax.plot(posNew[:,1], posNew[:,2])
    ax.scatter(posNew[0,1], posNew[0,2], marker = "+", color = "r", label = "first", zorder=2)
    ax.scatter(posNew[-1,1], posNew[-1,2], marker = "x", color = "r", label = "last", zorder=2)
    
    for i in range(14):
        plt.plot(alleyLines[i,:,0], alleyLines[i,:,1], color="k", alpha=0.5)
    
    exits = np.where(posNew[:,4] == 2)[0]
    for i,j in zip(exits,turns):
        ax.annotate(f"{int(posNew[i,3])} {j}", (posNew[i,1]+2, posNew[i,2]+2))
    ax.axis("equal")
    ax.legend()
    return posNew, turns


Rectangle = namedtuple("Rectangle", "xmin xmax ymin ymax")