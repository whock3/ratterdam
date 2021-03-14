# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:13:25 2020

@author: Ruo-Yah Lai
"""

from newAlleyBounds import R781, R808, R859, alleyInterType
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from collections import namedtuple
import csv
from string import ascii_uppercase
import ast


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

def crossBorder2(pos1, pos2, border1, border2):
    """
    border1 and border2: 3 segments that make up the N/S or E/W halves of the alley's border
                        arrays of shape 4 (points) x 2 (x and y)
    """
    for i in range(3):
        if intersect((border1[i,0],border1[i,1]), (border1[i+1,0],border1[i+1,1]), pos1[1:3], pos2[1:3]):
            border = 0
            break
        elif intersect((border2[i,0],border2[i,1]), (border2[i+1,0],border2[i+1,1]), pos1[1:3], pos2[1:3]):
            border = 2
            break
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

def alleyTransitions(pos, rat, graph=False, minTime=0.5e6):
    """
    Finds turns based on crossing alley bounds 
    rat : named tuple R781, R808, or R859 from newAlleyBounds


    Returns
    -------
    posNew : [ts, x, y, 1st or last pt of visit, alley/intersection]
    turns : [allo before turn, ego, allo after turn, ts of exit, ts of entry,
             alley exited, intersection, alley entered]

    """
    letters = [i for i in ascii_uppercase[:12]]
    a = np.hstack((np.arange(17).astype(str), letters))
    posList = deepcopy(pos)
    posNew = np.empty((0,5)) #ts, x, y, 1st or last pt of visit, alley/intersection
    alleys = []
    for i in a:
        alley = rat.alleyInterBounds[str(i)]
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
                posnew = np.column_stack((posList[idxPos], np.zeros(idxEnd-idxStart), np.full(idxEnd-idxStart,i)))
                if i.isdigit(): #if i is an alley
                    posnew[0,3] = 1 #1st point of this visit to alley
                    posnew[-1,3] = 2 #last point of this visit to alley
                else:
                    posnew[0,3] = 3 #1st point of this visit to intersection
                    posnew[-1,3] = 4 #last point of this visit to intersection
                posNew = np.vstack((posNew, posnew))
            
                posList = list(posList)
                del(posList[idxStart:idxEnd])
                posList = np.array(posList)
        
        if len(posList) == 0:
            break
    posNew = posNew[posNew[:,0].argsort()]
    posNewF = posNew[:,:4].astype(float) #ts, x, y, 1st/last pt of visit
    
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
    
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Alley transitions intersections.csv","r") as file:
        data_iter = csv.reader(file)
        inters = [data for data in data_iter]
        inters = np.array(inters)
    
    turns = np.empty((0,8)) #allo before turn, ego, allo after turn, ts of exit, ts of entry, alley exited, intersection, alley entered
    entries = np.where(posNewF[:,3] == 1)[0]
    
    for Exit in np.where(posNewF[:,3] == 2)[0]: #for each alley exit
        
        if len(entries[entries > Exit]) == 0:
            break
        entry = min(entries[entries > Exit])
        
        ego = trsnEgo[int(posNew[Exit,4]), int(posNew[entry,4])]
        if ego == 3: #back around turn
            allo1 = crossBorder(posNewF[Exit], posNewF[Exit+1], alleys[int(posNew[Exit,4])])+1 #0 index to 1 index
            if allo1 == 1000:
                continue
            allo2 = egoAllo(ego, allo1)
            if entry-Exit == 1: #no intersection data
                intersection = 999
            else:
                intersection = posNew[Exit+1,4]
        else:
            allo2 = trsnAllo[int(posNew[Exit,4]), int(posNew[entry,4])]
            intersection = inters[int(posNew[Exit,4]), int(posNew[entry,4])]
            allo1 = egoAllo(ego, allo2)
                
        turns = np.vstack((turns, (allo1, ego, allo2, posNew[Exit,0], posNew[entry,0], posNew[Exit,4], intersection, posNew[entry,4])))
    
    turns = turns[turns[:,3].astype(float).argsort()]
    
    if graph:
        fig, ax = plt.subplots()
        ax.plot(posNewF[:,1], posNewF[:,2])
        ax.scatter(posNewF[0,1], posNewF[0,2], marker = "+", color = "r", label = "first", zorder=2)
        ax.scatter(posNewF[-1,1], posNewF[-1,2], marker = "x", color = "r", label = "last", zorder=2)
    
        for i in range(14):
            plt.plot(rat.alleyLines[i,:,0], rat.alleyLines[i,:,1], color="k", alpha=0.5)
    
        #exits = np.where(turns[:,5].astype(int) == 2)[0]
        for turn in turns:
            a = np.where(posNewF[:,0] == float(turn[3]))[0]
            ax.annotate(f"{int(turn[5])} {turn[:3].astype(int)}", (posNewF[a,1], posNewF[a,2]))
            ax.axis("equal")
            ax.legend()
    return posNew, turns

def turnsInFieldIO(turns, unit, fieldLocs, subfield):
    """
    Finds the turns associated with a subfield, separated by into/inside/out of/through
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    
    Returns: [allo before turn, ego, allo after turn, ts of exit, ts of entry,
             alley exited, intersection, alley entered, into/inside/out of/through]
    """
    fieldLoc = ast.literal_eval(fieldLocs[subfield+1][1])
    turns2 = np.empty((0,9))
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if subfield == 3:
            print(allLoc)
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            #if neither the alley exited nor the alley entered is in the field, it's a through turn
            if len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])) and len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn,3))))
            #elif alley exited isn't in the field, it's an into turn
            elif len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])):
                turns2 = np.vstack((turns2, np.hstack((turn,0))))
            #elif alley entered isn't in the field, it's an out of turn
            elif len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn,2))))
            #else (alley exited and alley entered are both in the field), it's an inside turn
            else:
                turns2 = np.vstack((turns2, np.hstack((turn,1))))
    return turns2


def closestTurnToVisit(unit, rat, filename, df):
    """
    Finds the turns immediately before or after the 1st ts of each visit to a field
    """
    turns = alleyTransitions(unit.position, rat)
    #with open(df+filename, "r") as csv_file:
    #    data_iter = csv.reader(csv_file)
    #    fieldLocs = [data for data in data_iter]
    prevTurns = []
    nextTurns = []
    for subfield in range(len(unit.visits)):
        prevTurns.append([])
        nextTurns.append([])
        #turns2 = turnsInFieldIO(turns, unit, fieldLocs, subfield)
        turnsF = turns[:,:5].astype(float)
        
        for visit in unit.visits[subfield]:
            turnsBefore = turns[turnsF[:,4] < visit[0]] #ts of entry < 1st ts of visit
            turnsAfter = turns[turnsF[:,3] > visit[0]] #ts of exit > 1st ts of visit
            if len(turnsBefore) > 0: 
                prevTurns[-1].append(turnsBefore[-1])
            else:
                prevTurns[-1].append(np.empty((0,9)))
            if len(turnsAfter) > 0:
                nextTurns[-1].append(turnsAfter[0])
            else:
                nextTurns[-1].append(np.empty((0,9)))
                
            
    return prevTurns, nextTurns
            
        
def alleyTransitions2(pos, rat, graph=False, minTime=0.5e6):
    """
    Finds turns based on crossing alley bounds. Different ts.
    rat : named tuple R781, R808, or R859 from newAlleyBounds

    Returns
    -------
    posNew : [ts, x, y, 1st or last pt of visit, alley/intersection]
    turns : [allo before turn, ego, allo after turn, ts of entering 1st alley, 
             ts of exiting 2nd alley, alley exited, intersection, alley entered]

    """
    letters = [i for i in ascii_uppercase[:12]]
    a = np.hstack((np.arange(17).astype(str), letters))
    posList = deepcopy(pos)
    posNew = np.empty((0,5)) #ts, x, y, 1st or last pt of visit, alley/intersection
    alleys = []
    for i in a:
        alley = rat.alleyInterBounds[str(i)]
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
                posnew = np.column_stack((posList[idxPos], np.zeros(idxEnd-idxStart), np.full(idxEnd-idxStart,i)))
                if i.isdigit(): #if i is an alley
                    posnew[0,3] = 1 #1st point of this visit to alley
                    posnew[-1,3] = 2 #last point of this visit to alley
                else:
                    posnew[0,3] = 3 #1st point of this visit to intersection
                    posnew[-1,3] = 4 #last point of this visit to intersection
                posNew = np.vstack((posNew, posnew))
            
                posList = list(posList)
                del(posList[idxStart:idxEnd])
                posList = np.array(posList)
        
        if len(posList) == 0:
            break
    posNew = posNew[posNew[:,0].argsort()]
    posNewF = posNew[:,:4].astype(float) #ts, x, y, 1st/last pt of visit
    
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
    
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Alley transitions intersections.csv","r") as file:
        data_iter = csv.reader(file)
        inters = [data for data in data_iter]
        inters = np.array(inters)
    
    turns = np.empty((0,8)) #allo before turn, ego, allo after turn, ts of exit, ts of entry, alley exited, intersection, alley entered
    entries = np.where(posNewF[:,3] == 1)[0]
    exits = np.where(posNewF[:,3] == 2)[0]
    
    for exit1 in exits: #for each exit from alley1
        if len(entries[entries > exit1]) == 0:
            break
        
        entry1 = max(entries[entries < exit1]) #entry to alley1
        entry2 = min(entries[entries > exit1]) #entry to alley2
        exit2 = min(exits[exits > entry2]) #exit from alley2
        
        ego = trsnEgo[int(posNew[entry1,4]), int(posNew[exit2,4])]
        if ego == 3: #back around turn
            allo1 = crossBorder(posNewF[exit1], posNewF[exit1+1], alleys[int(posNew[exit1,4])])+1 #0 index to 1 index
            if allo1 == 1000:
                continue
            allo2 = egoAllo(ego, allo1)
            if entry2-exit1 == 1: #no intersection data
                intersection = 999
            else:
                intersection = posNew[exit1+1,4]
        else:
            allo2 = trsnAllo[int(posNew[entry1,4]), int(posNew[exit2,4])]
            intersection = inters[int(posNew[entry1,4]), int(posNew[exit2,4])]
            allo1 = egoAllo(ego, allo2)
                
        turns = np.vstack((turns, (allo1, ego, allo2, posNew[entry1,0], posNew[exit2,0], posNew[entry1,4], intersection, posNew[exit2,4])))
        
    if graph:
        fig, ax = plt.subplots()
        ax.plot(posNewF[:,1], posNewF[:,2])
        ax.scatter(posNewF[0,1], posNewF[0,2], marker = "+", color = "r", label = "first", zorder=2)
        ax.scatter(posNewF[-1,1], posNewF[-1,2], marker = "x", color = "r", label = "last", zorder=2)
    
        for i in range(14):
            plt.plot(rat.alleyLines[i,:,0], rat.alleyLines[i,:,1], color="k", alpha=0.5)
    
        #exits = np.where(turns[:,5].astype(int) == 2)[0]
        for turn in turns:
            a = np.where(posNewF[:,0] == float(turn[3]))[0]
            ax.annotate(f"{int(turn[5])} {turn[:3].astype(int)}", (posNewF[a,1], posNewF[a,2]))
            ax.axis("equal")
            ax.legend()
            
    return posNew, turns


#def turnsBasedOnFieldVisits(unit, posNew):
#    posNewF = posNew[:,:2].astype(float)
#    for subfield in range(len(unit.perimeters)):
#        
#        for i, visit in enumerate(unit.visits[subfield]):
#            #find which location 1st and last pts of visit are in
#            startBool = np.logical_and(posNewF[:,0]<=visit[0], visit[0]<=posNewF[:,1]) 
#            endBool = np.logical_and(posNewF[:,0]<=visit[-1], visit[-1]<=posNewF[:,1])
#            start = np.where(startBool)[0] #index in posNew
#            end = np.where(endBool)[0] #index in posNew
#            
#            if end-start != 0: #1st and last pt of visit are in different locations
#                #if one of the locations is an alley and the other is an intersection
#                if posNew[start,2].isdigit(): #start of visit is in alley
#                    alley1 = posNew[start,2]
#                else: #start of visit is in intersection
#                    intersection = posNew[start,2]
#                    if posNew[start-1,2].isdigit(): #the location visited before the visit to the field is an alley
#                        alley1 = posNew[start-1,2]
#                    else:
#                        print("Can't find the 1st alley", )
#                if posNew[start,2]
#            
#            else: #1st and last pt of visit are in the same visit to a location
#                if ~posNew[start,2].isdigit(): #if in intersection
#                    start = start-1
#                    end = end+1
#                    
#                else: #if in alley
#                    startPos = unit.position[np.where(unit.position[:,0] == visit[0])[0]]
#                    endPos = unit.position[np.where(unit.position[:,0] == visit[-1])[0]]
#                    alleyOrient = alleyInterType[posNew[start,2]][1]
#                    if alleyOrient == "horizontal":
#                        maxLongAxis = max(unit.perimeters[subfield][:,0]) #max in x-axis
#                        minLongAxis = min(unit.perimeters[subfield][:,0]) #min in x-axis
#                        
#                        #find which end of the field the start and end of the visit are closer to
#                        startSide = np.argmin(np.abs(np.array([maxLongAxis,minLongAxis])-startPos[1])) #0=max x side, 1=min x side
#                        endSide = np.argmin(np.abs(np.array([maxLongAxis,minLongAxis])-endPos[1]))
#                        
#                    elif alleyOrient == "vertical":
#                        maxLongAxis = max(unit.perimeters[subfield][:,1]) #max in y-axis
#                        minLongAxis = min(unit.perimeters[subfield][:,1]) #min in y-axis
#                    
#                        #find which end of the field the start and end of the visit are closer to
#                        startSide = np.argmin(np.abs(np.array([maxLongAxis,minLongAxis])-startPos[2])) #0=max y side, 1=min y side
#                        endSide = np.argmin(np.abs(np.array([maxLongAxis,minLongAxis])-endPos[2]))



Rectangle = namedtuple("Rectangle", "xmin xmax ymin ymax")


into0 = [[False, False, False, True],
         [False, True, False, False],
         [False, False, False, False],
         [False, False, False, True]]
into1 = [[False, False, False, True],
         [False, True, False, True],
         [False, False, False, False],
         [False, True, False, True]]
into2= [[True, False, False, False],
        [True, False, False, False],
        [False, False, False, False],
        [False, False, True, False]]
into3 = [[True, False, False, False],
         [True, False, True, False],
         [False, False, False, False],
         [True, False, True, False],]
into4 = [[False, True, False, True], 
         [False, True, False, False],
         [False, False, False, False],
         [False, False, False, True]]
into5 = [[True, False, False, False],
         [True, False, True, False],
         [False, False, False, False],
         [True, False, True, False],]
into6 = [[False, True, False, False], 
         [False, True, False, False],
         [False, False, False, False],
         [False, False, False, True]]
into7 = [[True, False, False, False],
         [False, False, True, False],
         [False, False, False, False],
         [True, False, False, False]]
into8 = [[False, True, False, False],
         [False, True, False, True],
         [False, False, False, False],
         [False, True, False, True]]
into9 = [[False, False, True, False],
         [False, False, True, False],
         [False, False, False, False],
         [True, False, False, False]]
into10 = [[False, True, False, False],
          [False, False, False, True],
          [False, False, False, False],
          [False, True, False, False]]
into11 = [[False, False, True, False],
         [True, False, True, False],
         [False, False, False, False],
         [True, False, True, False],]
into12 = [[False, True, False, True],
         [False, True, False, True],
         [False, False, False, False],
         [False, True, False, True]]
into13 = [[False, True, False, True],
          [False, False, False, True],
          [False, False, False, False],
          [False, True, False, False]]
into14 = [[False, False, True, False],
         [True, False, True, False],
         [False, False, False, False],
         [True, False, True, False],]
into15 = [[False, False, False, True],
          [False, False, False, True],
          [False, False, False, False],
          [False, True, False, False]]
into16 = [[False, False, True, False],
        [True, False, False, False],
        [False, False, False, False],
        [False, False, True, False]]


insideH = [[False, False, False, False],
           [False, False, False, False],
           [False, True, False, True],
           [False, False, False, False]]
insideV = [[False, False, False, False],
           [False, False, False, False],
           [True, False, True, False],
           [False, False, False, False],]


outof0 = [[False, True, False, False], 
         [False, False, True, False],
         [False, False, False, False],
         [False, False, True, False]]
outof1 = [[False, True, False, False],
          [True, False, True, False],
          [False, False, False, False],
          [True, False, True, False]]
outof2 = [[False, False, True, False],
          [False, True, False, False],
          [False, False, False, False],
          [False, True, False, False]]
outof3 = [[False, False, True, False],
          [False, True, False, True],
          [False, False, False, False],
          [False, True, False, True]]
outof4 = [[False, True, False, True], 
         [False, False, True, False],
         [False, False, False, False],
         [False, False, True, False]]
outof5 = [[False, False, True, False],
          [False, True, False, True],
          [False, False, False, False],
          [False, True, False, True]]
outof6 = [[False, False, False, True], 
         [False, False, True, False],
         [False, False, False, False],
         [False, False, True, False]]
outof7 = [[False, False, True, False],
          [False, False, False, True],
          [False, False, False, False],
          [False, False, False, True]]
outof8 = [[False, False, False, True],
          [True, False, True, False],
          [False, False, False, False],
          [True, False, True, False]]
outof9 = [[True, False, False, False],
          [False, False, False, True],
          [False, False, False, False],
          [False, False, False, True]]
outof10 = [[False, False, False, True],
           [True, False, False, False],
           [False, False, False, False],
           [True, False, False, False]]
outof11 = [[True, False, False, False],
          [False, True, False, True],
          [False, False, False, False],
          [False, True, False, True]]
outof12 = [[False, True, False, True],
          [True, False, True, False],
          [False, False, False, False],
          [True, False, True, False]]
outof13 = [[False, True, False, True],
           [True, False, False, False],
           [False, False, False, False],
           [True, False, False, False]]
outof14 = [[True, False, False, False],
          [False, True, False, True],
          [False, False, False, False],
          [False, True, False, True]]
outof15 = [[False, True, False, False],
           [True, False, False, False],
           [False, False, False, False],
           [True, False, False, False]]
outof16 = [[True, False, False, False],
          [False, True, False, False],
          [False, False, False, False],
          [False, True, False, False]]


possibleTurnsA = np.array([[into0,into1,into2,into3,into4,into5,into6,into7,into8,into9,into10,into11,into12,into13,into14,into15,into16],
                  [insideH,insideH,insideV,insideV,insideH,insideV,insideH,insideV,insideH,insideV,insideH,insideV,insideH,insideH,insideV,insideH,insideV],
                  [outof0,outof1,outof2,outof3,outof4,outof5,outof6,outof7,outof8,outof9,outof10,outof11,outof12,outof13,outof14,outof15,outof16],
                  [[[False,False,False,False] for _ in range(4)] for _ in range(17)]])


thruA = [[False, False, False, False],
         [False, True, False, False],
         [False, True, True, False],
         [False, False, True, False]]
thruB = [[False, True, False, True],
         [False, True, True, False],
         [False, True, True, True],
         [False, False, True, True]]
thruC = [[False, True, False, True],
         [False, True, True, False],
         [False, True, True, True],
         [False, False, True, True]]
thruD = [[False, False, False, False],
         [False, False, True, False],
         [False, False, True, True],
         [False, False, False, True]]
thruE = [[True, False, True, False],
         [True, True, False, False],
         [True, True, True, False],
         [False, True, True, False]]
thruF = [[True, True, True, True],
         [True, True, True, True],
         [True, True, True, True],
         [True, True, True, True]]
thruG = [[True, True, True, True],
         [True, True, True, True],
         [True, True, True, True],
         [True, True, True, True]]
thruH = [[True, False, True, False],
         [False, False, True, True],
         [True, False, True, True],
         [True, False, False, True]]
thruI = [[False, False, False, False],
         [True, False, False, False],
         [True, True, False, False],
         [False, True, False, False]]
thruJ = [[False, True, False, True],
         [True, False, False, True],
         [True, True, False, True],
         [True, True, False, False]]
thruK = [[False, True, False, True],
         [True, False, False, True],
         [True, True, False, True],
         [True, True, False, False]]
thruL = [[False, False, False, False],
         [False, False, False, True],
         [True, False, False, True],
         [True, False, False, False]]

possibleTurnsI = np.array([[[[False,False,False,False] for _ in range(4)] for _ in range(12)],
                  [[[False,False,False,False] for _ in range(4)] for _ in range(12)],
                  [[[False,False,False,False] for _ in range(4)] for _ in range(12)],
                  [thruA,thruB,thruC,thruD,thruE,thruF,thruG,thruH,thruI,thruJ,thruK,thruL]])
