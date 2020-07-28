# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 15:20:27 2020

@author: Ruo-Yah Lai
"""
import numpy as np
from RDP4_class import RDP4
from turn3 import turn3
from copy import deepcopy
import matplotlib.pyplot as plt

def turnsInField3(visits, position, subfield):
    """
    Finds the turns where at least 1 point is inside the subfield
    position: before RDP
    Returns: ts of 1st pt, before/during/after visit
    """
    RList = RDP4(position, 4.72).ResultList
    idx3, theta_sum2, idx2, idx3_2 = turn3(RList, np.pi/12, np.pi/4, 5*4.72, 15*4.72)
    directions = np.diff(RList[:, 1:3], axis=0)
    allo = np.arctan2(directions[idx2, 1], directions[idx2, 0]) % (2*np.pi)
    ts, ts2 = RList[idx3, 0], RList[idx2,0]
    #ts ts of 1st pt of each turn
    #ts2: ts of last pt of each turn
    j = 0
    js = np.array([])
    a = np.array([]) #before, during, or after the visit
    for i in range(len(visits[subfield])): #the ith visit
        while j < len(ts):
            if ts[j]-400000 > visits[subfield][i][-1]:
                break
            if visits[subfield][i][0] <= ts[j] <= visits[subfield][i][-1] and ts2[j] <= visits[subfield][i][-1]:
                js = np.hstack((js,j))
                a = np.hstack((a,1))
            elif visits[subfield][i][0] <= ts[j] <= visits[subfield][i][-1] and ts2[j] > visits[subfield][i][-1]:
                js = np.hstack((js,j))
                a = np.hstack((a,2))
            elif visits[subfield][i][0] <= ts2[j] <= visits[subfield][i][-1]:
                js = np.hstack((js,j))
                a = np.hstack((a,0))
            elif ts[j] < visits[subfield][i][0] and visits[subfield][i][-1] <= ts2[j]:
                js = np.hstack((js,j))
                a = np.hstack((a,3))
            elif visits[subfield][i][0] <= ts[j]-400000 <= visits[subfield][i][-1]:
                js = np.hstack((js,j))
                a = np.hstack((a,4))
            elif ts[j]-400000 < visits[subfield][i][0] and visits[subfield][i][-1] < ts2[j]:
                js = np.hstack((js,j))
                a = np.hstack((a,5))
            j += 1

    ts, theta_sum2, allo = ts.reshape((-1,1)), theta_sum2.reshape((-1,1)), allo.reshape((-1,1))
    timeAndThetas = np.hstack((ts, theta_sum2, allo))
    return timeAndThetas[js.astype(int)], a

def classifyTurns2(b1, b2, position, min_angle1=np.pi/12, min_angle2=np.pi/4, turn_range_start=5*4.72, vCutoff=15*4.72, vIncrease=1/7, epsilon=4.72):
    """
    9 ego and 12 allo turn directions
    """
    RList = RDP4(position[b1:b2], epsilon).ResultList
    idx3, theta_sum2, idx2, idx3_2 = turn3(RList, min_angle1, min_angle2, turn_range_start, vCutoff, vIncrease)
    
    ego_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if theta_sum2[i] < (5/12*np.pi):
            ego_turns[i] = 0 #FL
        elif (5/12*np.pi) <= theta_sum2[i] < (7/12*np.pi):
            ego_turns[i] = 1 #L
        elif (7/12*np.pi) <= theta_sum2[i] < (9/12*np.pi):
            ego_turns[i] = 2 #BL
        elif (9/12*np.pi) <= theta_sum2[i] < (11/12*np.pi):
            ego_turns[i] = 3 #LB
        elif (11/12*np.pi) <= theta_sum2[i] < (13/12*np.pi):
            ego_turns[i] = 4 #B
        elif (13/12*np.pi) <= theta_sum2[i] < (15/12*np.pi):
            ego_turns[i] = 5 #RB
        elif (15/12*np.pi) <= theta_sum2[i] < (17/12*np.pi):
            ego_turns[i] = 6 #BR
        elif (17/12*np.pi) <= theta_sum2[i] < (19/12*np.pi):
            ego_turns[i] = 7 #R
        elif (19/12*np.pi) <= theta_sum2[i] < (21/12*np.pi):
            ego_turns[i] = 8 #FR
    ego_turns = np.reshape(ego_turns, (len(idx3), 1))
            
    directions = np.diff(RList[:, 1:3], axis=0)
    allo = np.arctan2(directions[idx2, 1], directions[idx2, 0]) % (2*np.pi) #in RList frame, length of idx3
    allo_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if (1/12*np.pi) <= allo[i] < (3/12*np.pi):
            allo_turns[i] = 0 #ENE
        elif (3/12*np.pi) <= allo[i] < (5/12*np.pi):
            allo_turns[i] = 1 #NNE
        elif (5/12*np.pi) <= allo[i] < (7/12*np.pi):
            allo_turns[i] = 2 #N
        elif (7/12*np.pi) <= allo[i] < (9/12*np.pi):
            allo_turns[i] = 3 #NNW
        elif (9/12*np.pi) <= allo[i] < (11/12*np.pi):
            allo_turns[i] = 4 #WNW
        elif (11/12*np.pi) <= allo[i] < (13/12*np.pi):
            allo_turns[i] = 5 #W
        elif (13/12*np.pi) <= allo[i] < (15/12*np.pi):
            allo_turns[i] = 6 #WSW
        elif (15/12*np.pi) <= allo[i] < (17/12*np.pi):
            allo_turns[i] = 7 #SSW
        elif (17/12*np.pi) <= allo[i] < (19/12*np.pi):
            allo_turns[i] = 8 #S
        elif (19/12*np.pi) <= allo[i] < (21/12*np.pi):
            allo_turns[i] = 9 #SSE
        elif (21/12*np.pi) <= allo[i] < (23/12*np.pi):
            allo_turns[i] = 10 #ESE
        else:
            allo_turns[i] = 11 #E
    allo_turns = np.reshape(allo_turns, (len(idx3), 1))
    
    times = np.reshape(RList[idx3, 0], (len(idx3), 1))
    return np.hstack((times, ego_turns, allo_turns)), RList[idx3_2,0], RList[idx2,0]


def linearize(pos, spikes, timeAndThetas, turnBinSize = 1/4*np.pi, RMrange=7*4.72):
    """
    Linearize 2D turns
    
    pos: [ts,x,y] before RDP
    timeAndThetas: [ts, theta, theta] from turnsInField3
    RMrange: +/- n unit distance of the 1st point of turns to make rate map
    
    Returns: 
        Ps: list of 22 lists each with arrays where each array is 1 turn
    pos (before RDP) within RMrange of the 1st point of turns, 
    centered on the 1st point of turns,
    normalized to be at most [-1,1] in the x axis
        Ss: same as Ps but with spikes
    """
    Ps = [[] for _ in range(int(1 + (3/2*np.pi)/turnBinSize + 2*np.pi/turnBinSize))]
    Ss = [[] for _ in range(int(1 + (3/2*np.pi)/turnBinSize + 2*np.pi/turnBinSize))]
    for i in range(len(timeAndThetas)): #adapted from turn3
        k = np.where(pos == timeAndThetas[i,0])[0] #1st point of the turn, in pos frame
        j = deepcopy(k)
        dist = 0
        while dist < RMrange: #checking points after the turn
            idx_end = deepcopy(j)
            j += 1
            if j > len(pos)-1:
                break
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
        
        j = k - 1 #one point before the turn
        if j >= 0:
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
            idx_start = j+1
        else:
            dist = RMrange+10
            idx_start = 0
        while dist < RMrange: #checking points before the turn
            idx_start = deepcopy(j)
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
        
        P = np.array([pos[k,0], 0])
        lastDist = 0
        S = np.empty((0,2))
        spike_idx = np.where(spikes == pos[j,0])[0]
        for l in spike_idx:
            S = np.vstack((S,np.array([spikes[l,0], 0])))
        for j in range(k[0]+1, idx_end[0]+1):
            dist = np.sqrt((pos[j, 1]-pos[j-1, 1])**2 + (pos[j, 2]-pos[j-1, 2])**2)
            P = np.vstack((P,np.array([pos[j,0], lastDist+dist])))
            spike_idx = np.where(spikes == pos[j,0])[0]
            for l in spike_idx:
                S = np.vstack((S,np.array([spikes[l,0], lastDist+dist])))
            lastDist = lastDist + dist
        
        lastDist = 0
        for j in range(k[0]-1, idx_start[0]-1, -1):
            dist = np.sqrt((pos[j, 1]-pos[j+1, 1])**2 + (pos[j, 2]-pos[j+1, 2])**2)
            P = np.vstack((np.array([pos[j,0], -lastDist-dist]),P))
            spike_idx = np.where(spikes == pos[j,0])[0]
            for l in reversed(spike_idx):
                S = np.vstack((np.array([spikes[l,0], -lastDist-dist]),S))
            lastDist = lastDist + dist
        
        maxDist = max(np.max(P[:,1]), abs(np.min(P[:,1])))
        P = np.divide(P, np.array([1,maxDist]))
        S = np.divide(S, np.array([1,maxDist]))
        
        Ps[0].append(P)
        Ss[0].append(S)
        egoBins = np.arange(1/4*np.pi, 7/4*np.pi, turnBinSize)
        for j, egoBin in enumerate(egoBins):
            if egoBin <= timeAndThetas[i,1] < (egoBin + turnBinSize):
                Ps[j+1].append(P)
                Ss[j+1].append(S)
        alloBins = np.arange(0, 2*np.pi, turnBinSize)
        for j, alloBin in enumerate(alloBins):
            if alloBin <= timeAndThetas[i,2] < (alloBin + turnBinSize):
                Ps[j+1+int((3/2*np.pi)/turnBinSize)].append(P)
                Ss[j+1+int((3/2*np.pi)/turnBinSize)].append(S)
        
        """
        if timeAndTurns[i,1] == 0:
            Ps[1].append(P)
            Ss[1].append(S)
            #shifted_pos2[1] = np.vstack((shifted_pos2[1],srP))
        elif timeAndTurns[i,1] == 1:
            Ps[2].append(P)
            Ss[2].append(S)
            #shifted_pos2[2] = np.vstack((shifted_pos2[2],srP))
        elif timeAndTurns[i,1] == 2:
            Ps[3].append(P)
            Ss[3].append(S)
            #shifted_pos2[3] = np.vstack((shifted_pos2[3],srP))
        elif timeAndTurns[i,1] == 3:
            Ps[4].append(P)
            Ss[4].append(S)
        elif timeAndTurns[i,1] == 4:
            Ps[5].append(P)
            Ss[5].append(S)
        elif timeAndTurns[i,1] == 5:
            Ps[6].append(P)
            Ss[6].append(S)
        elif timeAndTurns[i,1] == 6:
            Ps[7].append(P)
            Ss[7].append(S)
        elif timeAndTurns[i,1] == 7:
            Ps[8].append(P)
            Ss[8].append(S)
        elif timeAndTurns[i,1] == 8:
            Ps[9].append(P)
            Ss[9].append(S)
        
        #allo turns
        if timeAndTurns[i,2] == 0:
            Ps[10].append(P)
            Ss[10].append(S)
            #shifted_pos2[4] = np.vstack((shifted_pos2[4],shiftedP))
        elif timeAndTurns[i,2] == 1:
            Ps[11].append(P)
            Ss[11].append(S)
        elif timeAndTurns[i,2] == 2:
            Ps[12].append(P)
            Ss[12].append(S)
        elif timeAndTurns[i,2] == 3:
            Ps[13].append(P)
            Ss[13].append(S)
        elif timeAndTurns[i,2] == 4:
            Ps[14].append(P)
            Ss[14].append(S)
        elif timeAndTurns[i,2] == 5:
            Ps[15].append(P)
            Ss[15].append(S)
        elif timeAndTurns[i,2] == 6:
            Ps[16].append(P)
            Ss[16].append(S)
        elif timeAndTurns[i,2] == 7:
            Ps[17].append(P)
            Ss[17].append(S)
        elif timeAndTurns[i,2] == 8:
            Ps[18].append(P)
            Ss[18].append(S)
        elif timeAndTurns[i,2] == 9:
            Ps[19].append(P)
            Ss[19].append(S)
        elif timeAndTurns[i,2] == 10:
            Ps[20].append(P)
            Ss[20].append(S)
        elif timeAndTurns[i,2] == 11:
            Ps[21].append(P)
            Ss[21].append(S)
        """
    return Ps, Ss


def makeRM1D(pos, spikes, nBins):
    hs = np.histogram(spikes[:,1], bins=nBins, range=(-1,1))[0]
    hp = np.histogram(pos[:,1], bins=nBins, range=(-1,1))[0]
    n = (hs*np.reciprocal(hp.astype(float)))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(hp==0)] = np.nan
    return n


def graphRM1D(pos, spikes, suptitle, turnBinSize=1/4*np.pi, percentile=99):
    """ 
    graphRM2 adapted for 1D avg RM for turns
    pos and spikes from linearize
    """
    if turnBinSize == 1/8*np.pi:
        nBins = 10
    else:
        nBins = 15
    n = np.empty((0,nBins))
    stds = []
    titles2 = []
    for i in range(int(1 + (3/2*np.pi)/turnBinSize + 2*np.pi/turnBinSize)):
        n1 = []
        for j in range(len(pos[i])):
            n1.append(makeRM1D(pos[i][j], spikes[i][j], nBins))
        n2 = np.nanmean(n1, axis=0)
        stds.append(np.nanstd(n1, axis=0))
        #n2 = weird_smooth(n2,smoothing_2d_sigma)
        #n2[np.where(ho==0)] = np.nan
        if len(n1) == 0:
            n2 = np.zeros((1,nBins))
        n = np.vstack((n,n2))
        titles2.append(f"n = {len(n1)}")
    
    if turnBinSize == 1/8*np.pi:
        fig, axs = plt.subplots(5,6,figsize=(12,10))
        turnBinD = round(turnBinSize/np.pi*180,1)
    elif turnBinSize == 1/6*np.pi:
        fig, axs = plt.subplots(4,6,figsize=(12,8))
        turnBinD = int(turnBinSize/np.pi*180)
    elif turnBinSize == 1/4*np.pi:
        fig, axs = plt.subplots(3,5,figsize=(10,6))
        turnBinD = int(turnBinSize/np.pi*180)
    egoTitles = ["Ego: " + str(Bin) + "°-" + str(Bin+turnBinD) + "°" 
                 for Bin in np.arange(45, 315, turnBinD)]
    alloTitles = ["Allo: " + str(Bin) + "°-" + str(Bin+turnBinD) + "°" 
                 for Bin in np.arange(0, 360, turnBinD)]
    titles = ["All turns"] + egoTitles + alloTitles
    #titles = ["All turns", "Front-left", "Left", "Back-left", "Left-back", 
    #          "Back", "Right-back", "Back-right", "Right", "Front-right",
    #          "ENE", "NNE", "N", "NNW", "WNW", "W", "WSW", "SSW", "S", "SSE",
    #          "ESE", "E"]
    ymax = np.nanpercentile(n, percentile)
       
    x = np.linspace(-1+2/nBins/2, 1-2/nBins/2, nBins)
    for j in range(2):
        if j == 0:
            errorbar = False
        else:
            errorbar = True
            if turnBinSize == 1/8*np.pi:
                fig, axs = plt.subplots(5,6,figsize=(12,10))
            elif turnBinSize == 1/6*np.pi:
                fig, axs = plt.subplots(4,6,figsize=(12,8))
            elif turnBinSize == 1/4*np.pi:
                fig, axs = plt.subplots(3,5,figsize=(10,6))
        fig.suptitle(suptitle, y=1.02)
        for i,ax in enumerate(axs.flatten()):
            if i == int(1 + (3/2*np.pi)/turnBinSize + 2*np.pi/turnBinSize):
                break
            ax.set_title(titles[i] + "\n" + titles2[i])
            ax.set_ylim(0, ymax)
            ax.set_xlim(-1, 1)
            ax.set_ylabel("Rate (Hz)")
            if len(n[i]) > 0:
                if errorbar:
                    ax.errorbar(x, n[i], yerr=stds[i])
                else:
                    ax.plot(x, n[i])
        fig.tight_layout()