# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:05:21 2020

@author: Ruo-Yah Lai
"""
import matplotlib.pyplot as plt
from bisect import bisect
import numpy as np
import matplotlib.patches as patches
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
from RDP4_class import RDP4
from turn3 import turn3
from newAlleyBounds import cityBlocks, alleyLines



def graph(tath, pos): 
    """
    Graphs each turn on a subplot
    """
    cm = np.array([1, 4.72, 4.72])
    pos = pos/cm[None,:]
    fig, axs = plt.subplots(4, 4, figsize=(8,8))
    axs = axs.flatten()

    for i,turn in enumerate(tath[(tath[:,3] == 0) | (tath[:,3] == 3)]):
        a = bisect(pos[:,0], turn[0]-5e5)
        a1 = bisect(pos[:,0], turn[0])
        a2 = bisect(pos[:,0], turn[4])
        
        axs[i].plot(pos[a:a2+1,1], pos[a:a2+1,2], zorder=1)
        axs[i].scatter(pos[a,1], pos[a,2], marker = "+", color = "r", label = "first", zorder=2)
        axs[i].scatter(pos[a2,1], pos[a2,2], marker = "x", color = "r", label = "last", zorder=2)
        axs[i].scatter(pos[a1,1], pos[a1,2], s = 16, color = "r", label='turning points', zorder=2)
        
        axs[i].set_title(turn[0])
        axs[i].axis("equal")
    fig.tight_layout()
    
    
def graph2(tath, pos, border):
    """
    Graphs all turns together with the field border
    """
    cm = np.array([1, 4.72, 4.72])
    pos = pos/cm[None,:]
    fig, ax = plt.subplots()
    ax.plot(border[:,0], border[:,1])
    for i,turn in enumerate(tath[(tath[:,3] == 0) | (tath[:,3] == 3)]):
    #for i,turn in enumerate(tath[tath[:,3] == 1]):
        a = bisect(pos[:,0], turn[0]-5e5)
        a1 = bisect(pos[:,0], turn[0])
        a2 = bisect(pos[:,0], turn[4])
        
        plt.plot(pos[a:a2+1,1], pos[a:a2+1,2], zorder=1)
        ax.scatter(pos[a,1], pos[a,2], marker = "+", color = "b", label = "first", zorder=2)
        ax.scatter(pos[a2,1], pos[a2,2], marker = "x", color = "g", label = "last", zorder=2)
        ax.scatter(pos[a1,1], pos[a1,2], s = 16, color = "r", label='turning points', zorder=2)
        
        ax.axis("equal")
        #ax.set_title(turn[0])


def graphLive(pos, border=np.empty((0,2)), start=0):
    """
    Live plot of the rat running
    """
    cm = np.array([1, 4.72, 4.72])
    pos = pos/cm[None,:]
    border = border/4.72
    
    x = np.empty(0)
    y = np.empty(0)
    i = start
    while True:
        for j in range(5):
            plt.cla()
            plt.axis([0,130,0,100])
            plt.plot(border[:,0], border[:,1])
            for k in range(len(alleyLines)):
                plt.plot(alleyLines[k][:,0]/4.72, alleyLines[k][:,1]/4.72)
    
            #for k in range(6):
            #    plt.plot(cityBlocks[k,:,0]/4.72, cityBlocks[k,:,1]/4.72)
                            
            x = np.hstack((x, pos[i*150+j*30:i*150+j*30+30, 1]))
            y = np.hstack((y, pos[i*150+j*30:i*150+j*30+30, 2]))
            plt.scatter(x, y, s=2)
            plt.title(pos[i*150+j*30+30,0])
        
            plt.pause(0.3)
        
        a = input("Press Enter to continue, s to stop, or b to go back once: ")
        if a == "s":
            break
        elif a == "b":
            i -= 1
            x = np.delete(x, np.s_[-150:])
            y = np.delete(y, np.s_[-150:])
        else:
            if i > start:
                x = np.delete(x, np.s_[0:150])
                y = np.delete(y, np.s_[0:150])
        i += 1
    print("Last i =", i)


def graphTurns(Rpos, border, tath, nextTurn, title=""):
    """
    Graphs turns associated with a field and labels them as 
    into, inside, or out from the field
    Rpos: position after RDP
    tath: from turnsInField3, only the ones associated with a certain visit
    """
    cm = np.array([1, 4.72, 4.72])
    Rpos = Rpos/cm[None,:]
    border = border/4.72
    
    fig, ax = plt.subplots()
    for i in range(6):
        plt.plot(cityBlocks[i,:,0]/4.72, cityBlocks[i,:,1]/4.72)
    ax.plot(border[:,0], border[:,1])
    first = np.where(Rpos[:,0] == tath[0,0])[0]
    last = np.where(Rpos[:,0] == tath[-1,0])[0]
    ax.plot(Rpos[int(first-5):int(last+5),1], Rpos[int(first-5):int(last+5),2], zorder=1)
    
    ax.scatter(Rpos[first-5,1], Rpos[first-5,2], marker = "+", color = "r", label = "first", zorder=2)
    ax.scatter(Rpos[last+4,1], Rpos[last+4,2], marker = "x", color = "r", label = "last", zorder=2)
    ax.axis("equal")
    ax.set_title(title)
    
    for j,i in enumerate(tath):
        idxFirst = np.where(Rpos[:,0] == i[0])
        idxLast = np.where(Rpos[:,0] == i[4])
        ax.scatter(Rpos[idxFirst,1], Rpos[idxFirst,2], color = "r", label="first of turns", zorder=2)
        ax.scatter(Rpos[idxLast,1], Rpos[idxLast, 2], s = 16, color = "b", label="last of turns", zorder=2)
        
        if i[0] == tath[0,0]:
            ax.legend()
        
        if i[3] == 1:
            ax.annotate("inside", (Rpos[idxFirst,1]+0.5, Rpos[idxFirst,2]+0.5))
        elif (i[3] == 0) | (i[3] == 3):
            ax.annotate("into", (Rpos[idxFirst,1]+0.5, Rpos[idxFirst,2]+0.5))
        else:
            ax.annotate("out from", (Rpos[idxFirst,1]+0.5, Rpos[idxFirst,2]+0.5))
    
    idxNextTurn = np.where(Rpos[:,0] == nextTurn[0])
    ax.scatter(Rpos[idxNextTurn,1], Rpos[idxNextTurn,2], color="g", zorder=2)

def graphVisit(visit, pos, border, title=""):
    """
    Graphs a visit to a field
    """
    cm = np.array([1, 4.72, 4.72])
    pos = pos/cm[None,:]
    border = border/4.72
    
    fig, ax = plt.subplots()
    ax.set_xlim(60, 130)
    ax.plot(border[:,0], border[:,1])
    start = np.where(pos[:,0] == visit[0])[0]
    end = np.where(pos[:,0] == visit[-1])[0]
    ax.plot(pos[int(start):int(end+1),1], pos[int(start):int(end+1),2])
    ax.set_title(title)
      

def turnsInField5(visits, position):
    """
    Finds the turns associated with the subfield
    If the ith turn is outside and the (i+1)th turn is inside, 
        the ith turn is into
    If the ith turn is inside and the (i+1)th turn is outside, 
        the ith turn is out from
    visits: ts of visits to the subfield of interest
    position: before RDP
    Returns: ts of 1st pt, ego, allo, into/inside/out from visit, ts of last pt
    """
    Rpos = RDP4(position, 4.72).ResultList
    idx3, theta_sum2, idx2, idx3_2 = turn3(Rpos, np.pi/12, np.pi/4, 5*4.72, 15*4.72)
    directions = np.diff(Rpos[:, 1:3], axis=0)
    allo = np.arctan2(directions[idx2, 1], directions[idx2, 0]) % (2*np.pi)
    ts, ts2 = Rpos[idx3,0], Rpos[idx2,0]
    #ts ts of 1st pt of each turn
    #ts2: ts of last pt of each turn
    j = 0
    js = np.array([])
    a = np.array([]) #into, inside, out from the visit
    b = np.array([]) #index of the visit
    c = np.empty((0,3)) #last or next turn
    for i in range(len(visits)): #the ith visit
        while j < len(ts):
            if ts[j] > visits[i][-1]:
                break
            if visits[i][0] <= ts[j] <= visits[i][-1]:
                if ts[j+1] > visits[i][-1]:# and \
                    #30*4.72 > np.sqrt((Rpos[idx3[j], 1]-Rpos[idx3[j+1], 1])**2 
                    #                  + (Rpos[idx3[j], 2]-Rpos[idx3[j+1], 2])**2):
                    js = np.hstack((js,j))
                    a = np.hstack((a,2)) #out from
                    b = np.hstack((b,i))
                    c = np.vstack((c,Rpos[idx3[j+1]]))
                elif ts2[j] <= visits[i][-1]:
                    js = np.hstack((js,j))
                    a = np.hstack((a,1)) #inside
                    b = np.hstack((b,i))
                    c = np.vstack((c,Rpos[idx3[j+1]]))
                else:
                    js = np.hstack((js,j))
                    a = np.hstack((a,2)) #out from
                    b = np.hstack((b,i))
                    c = np.vstack((c,Rpos[idx3[j+1]]))
            else:
                if visits[i][0] <= ts[j+1] <= visits[i][-1]:# and \
                    #30*4.72 > np.sqrt((Rpos[idx3[j], 1]-Rpos[idx3[j+1], 1])**2 
                    #                  + (Rpos[idx3[j], 2]-Rpos[idx3[j+1], 2])**2):
                    js = np.hstack((js,j))
                    a = np.hstack((a,0)) #into
                    b = np.hstack((b,i))
                    c = np.vstack((c,Rpos[idx3[j+1]]))
            j += 1

    #rotate allo so that 0 = N, np.pi/2 = E
    allo = -(allo - np.pi/2) % (np.pi*2)
    timeAndThetas = np.column_stack((ts[js.astype(int)], theta_sum2[js.astype(int)], allo[js.astype(int)], a, ts2[js.astype(int)], b))
    return timeAndThetas, c

    

def rect(bounds, color):
    bottom_left = [bounds[0][0]/4.72, bounds[1][0]/4.72]
    top_right = [bounds[0][1]/4.72, bounds[1][1]/4.72]
    return patches.Rectangle(tuple(bottom_left), (top_right[0]-bottom_left[0]), 
                             (top_right[1]-bottom_left[1]), color=color, alpha=0.4)
