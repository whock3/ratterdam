# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:13:25 2020

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.patches as patches
import numpy as np
from copy import deepcopy
from collections import namedtuple
import csv
from string import ascii_uppercase
import ast
from bisect import bisect_left


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
    
    
#for alleys then intersections:
    #find point that are inside
    #getVisits
    #all points between 1st and last pt of a visit is in that location
    #put into a new array with [ts,x,y,location,whether it's the last pt of the visit]
    #remove these points from the overall list
#sort the array by ts
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
    
    with open("E:\\Ratterdam\\RYL\\Alley_transitions_egocentric.csv","r") as file:
        data_iter = csv.reader(file)
        trsnEgo = [data for data in data_iter]
        trsnEgo = np.array(trsnEgo)
        trsnEgo = trsnEgo.astype(int)
    
    with open("E:\\Ratterdam\\RYL\\Alley_transitions_allocentric.csv","r") as file:
        data_iter = csv.reader(file)
        trsnAllo = [data for data in data_iter]
        trsnAllo = np.array(trsnAllo)
        trsnAllo = trsnAllo.astype(int)
    
    with open("E:\\Ratterdam\\RYL\\Alley_transitions_intersections.csv","r") as file:
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


def turnsInField(turns, spikes, filename, subfield, title=""):
    """
    Graphs the turns associated with a subfield
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    """
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/"+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    fieldLoc = ast.literal_eval(data[subfield+1][1])
    
    turns2 = np.empty((0,6)) #allo before turn, ego, allo after turn, ts of exit, ts of entry, into/inside/out of
    rates = np.empty(0)
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            #if alley exited isn't in the field, it's an into turn
            if len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],0))))
            #elif alley entered isn't in the field, it's an out of turn
            elif len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],2))))
            #else (alley exited and alley entered are both in the field), it's an inside turn
            else:
                turns2 = np.vstack((turns2, np.hstack((turn[:5],1))))
    
            start = bisect_left(spikes[:,0], float(turn[3]))
            end = bisect_left(spikes[:,0], float(turn[4]))
            rate = (end-start) / (float(turn[4])-float(turn[3])) * 1e6
            rates = np.hstack((rates,rate))
            
            
    ebins = np.arange(1,5)
    abins = np.arange(1,5)
    meanRates = np.zeros((4,4)) #number of ebins x number of abins
    Ns = np.empty((4,4))
    for ebin in ebins:
        for abin in abins:
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin)
            meanRates[ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[ebin-1,abin-1] = np.sum(inBin)
    
    fig, ax = plt.subplots()
    vmax = np.nanpercentile(meanRates,95)
    ax.set_title(title + "\n" + f"95th percentile of firing rate = {round(vmax,1)}")
    plt.xticks(np.arange(1.5,5.5,1), ["N", "E", "S", "W"])
    plt.yticks(np.arange(1.5,5.5,1), ["S", "R", "B", "L"])
    ax.set_xlim(1, 5)
    ax.set_ylim(1, 5)
    ax.set_xlabel("Allocentric turn direction")
    ax.set_ylabel("Egocentric turn direction")
    
    meanRates2 = deepcopy(meanRates)
    meanRates2[np.where(meanRates != meanRates)] = 0
    norm = Normalize(vmax=vmax)
    m = cm.ScalarMappable(norm=norm, cmap="jet")
    colors = m.to_rgba(meanRates2)
    #alphas = Ns/np.percentile(Ns,90) *0.8 + 0.2
    #alphas[np.where(alphas > 1)] = 1
    #colors[:,:,3] = alphas
    colors[np.where(meanRates != meanRates)] = 0
    for i in range(colors.shape[0]):
        for j in range(colors.shape[1]):
            patch = patches.Rectangle((abins[j],ebins[i]), 1, 1, color=colors[i,j,:])
            ax.annotate(int(Ns[i,j]), (j+1.5,i+1.5), color=(0.5,0.5,0.5))
            ax.add_patch(patch)
    
    axcb = fig.colorbar(m)
    axcb.set_label("Rate (Hz)")


def turnsInFieldIO(turns, spikes, filename, subfield, title=""):
    """
    Graphs the turns associated with a subfield, separated by into/inside/out of
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    """
    with open("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/"+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    fieldLoc = ast.literal_eval(data[subfield+1][1])
    
    turns2 = np.empty((0,6)) #allo before turn, ego, allo after turn, ts of exit, ts of entry, into/inside/out of
    rates = np.empty(0)
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            #if alley exited isn't in the field, it's an into turn
            if len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],0))))
            #elif alley entered isn't in the field, it's an out of turn
            elif len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],2))))
            #else (alley exited and alley entered are both in the field), it's an inside turn
            else:
                turns2 = np.vstack((turns2, np.hstack((turn[:5],1))))
            if turn[1] == "3" and turn[2] == "3":
                print(turn)
            start = bisect_left(spikes[:,0], float(turn[3]))
            end = bisect_left(spikes[:,0], float(turn[4]))
            rate = (end-start) / (float(turn[4])-float(turn[3])) * 1e6
            rates = np.hstack((rates,rate))
            
            
    ebins = np.arange(1,5)
    abins = np.arange(1,5)
    meanRates = np.zeros((3,4,4)) #3 x number of ebins x number of abins
    Ns = np.empty((3,4,4))
    for ebin in ebins:
        for abin in abins:
            #into
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 0)
            meanRates[0,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[0,ebin-1,abin-1] = np.sum(inBin)
            
            #inside
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 1)
            meanRates[1,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[1,ebin-1,abin-1] = np.sum(inBin)
            
            #out of
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 2)
            meanRates[2,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[2,ebin-1,abin-1] = np.sum(inBin)
    
    fig, axs = plt.subplots(1,3,figsize=(10,4))
    vmax = np.nanpercentile(meanRates,95)
    axs[1].set_title(title + "\n" + f"95th percentile of firing rate = {round(vmax,1)}"
                     + "\n\nInside the field")
    norm = Normalize(vmax=vmax)
    m = cm.ScalarMappable(norm=norm, cmap="jet")
    axs[0].set_title("Into the field")
    axs[2].set_title("Out of the field")
    for i in range(3):
        axs[i].set_xticks(np.arange(1.5, 5.5, 1))
        axs[i].set_xticklabels(["N", "E", "S", "W"])
        axs[i].set_yticks(np.arange(1.5, 5.5, 1))
        axs[i].set_yticklabels(["S", "R", "B", "L"])
        axs[i].set_xlim(1, 5)
        axs[i].set_ylim(1, 5)
        axs[i].set_xlabel("Allocentric turn direction")
        axs[i].set_ylabel("Egocentric turn direction")
    
        meanRates2 = deepcopy(meanRates[i])
        meanRates2[np.where(meanRates[i] != meanRates[i])] = 0
        colors = m.to_rgba(meanRates2)
        #alphas = Ns/np.percentile(Ns,90) *0.8 + 0.2
        #alphas[np.where(alphas > 1)] = 1
        #colors[:,:,3] = alphas
        colors[np.where(meanRates[i] != meanRates[i])] = 0
        for j in range(colors.shape[0]):
            for k in range(colors.shape[1]):
                patch = patches.Rectangle((abins[k],ebins[j]), 1, 1, color=colors[j,k,:])
                axs[i].annotate(int(Ns[i,j,k]), (k+1.5,j+1.5), color=(0.5,0.5,0.5))
                axs[i].add_patch(patch)
    
    axcb = fig.colorbar(m)
    axcb.set_label("Rate (Hz)")
    fig.tight_layout()

def bulkGraphs(turns,spikes1,spikes2,spikes4,spikes5,spikes6,spikes7,spikes8,timestamp):
    """
    Makes multiple graphs using turnsInFieldIO
    timestamp: current timestamp
    """
    A = [[spikes1,"unit1",0], [spikes2,"unit2",0], [spikes2,"unit2",3], 
         [spikes4,"unit4",0], [spikes5,"unit5",0], [spikes5,"unit5",1],
         [spikes5,"unit5",2], [spikes5,"unit5",3], [spikes6,"unit6",0], 
         [spikes6,"unit6",1], [spikes6,"unit6",2], [spikes6,"unit6",3], 
         [spikes7,"unit7",0], [spikes7,"unit7",1], [spikes8,"unit8",1],
         [spikes8,"unit8",3]]
    titles = ["1.1 subfield 0", "1.2 subfield 0", "1.2 subfield 3",
              "1.4 subfield 0", "1.5 subfield 0", "1.5 subfield 1", 
              "1.5 subfield 2", "1.5 subfield 3", "1.6 subfield 0", 
              "1.6 subfield 1", "1.6 subfield 2", "1.6 subfield 3",
              "1.7 subfield 0", "1.7 subfield 1", "1.8 subfield 1",
              "1.8 subfield 3"]
    for i,a in enumerate(A):
        turnsInFieldIO(turns, a[0], a[1]+" locations", a[2], "R859 D3 T6 "+titles[i])
        plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                    + timestamp + " - " + titles[i] + ".png")

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
