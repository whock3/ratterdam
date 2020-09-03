# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 10:23:45 2020

@author: Ruo-Yah Lai
"""

import numpy as np
import matplotlib.pyplot as plt
from organized import classifyTurns
from bisect import bisect
from copy import deepcopy


def barRateNearTurns(spikes, position, b1, b2, suptitle):
    """
    Bar graphs of firing rates within 2s of each type of ego and allo turns
    spikes : after adjusting to ts of pos
    b1, b2 = start and end of segment that goes into RDP
    """
    TimeAndTurns = classifyTurns(b1,b2,position)
    starts_ends = np.empty((0,1))
    for i in range(len(TimeAndTurns)):
        start = bisect(spikes[:,0], (TimeAndTurns[i,0] - 1e6))
        end = bisect(spikes[:,0], (TimeAndTurns[i,0] + 1e6))
        starts_ends = np.vstack((starts_ends, np.array([end-start])))
        #starts_ends = indices in spikes for 2s before and after turns
    TimeAndTurns = np.hstack((TimeAndTurns, starts_ends/2)) #ts, ego, allo, # of spikes within 2s
    
    #egocentric
    avgSpikesL = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,1] == 0)[0]][:,3])
    stSpikesL = np.std(TimeAndTurns[np.where(TimeAndTurns[:,1] == 0)[0]][:,3])
    avgSpikesB = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,1] == 1)[0]][:,3])
    stSpikesB = np.std(TimeAndTurns[np.where(TimeAndTurns[:,1] == 1)[0]][:,3])
    avgSpikesR = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,1] == 2)[0]][:,3])
    stSpikesR = np.std(TimeAndTurns[np.where(TimeAndTurns[:,1] == 2)[0]][:,3])
    
    #allocentric
    avgSpikesN = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,2] == 0)[0]][:,3])
    stSpikesN = np.std(TimeAndTurns[np.where(TimeAndTurns[:,2] == 0)[0]][:,3])
    avgSpikesE = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,2] == 1)[0]][:,3])
    stSpikesE = np.std(TimeAndTurns[np.where(TimeAndTurns[:,2] == 1)[0]][:,3])
    avgSpikesS = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,2] == 2)[0]][:,3])
    stSpikesS = np.std(TimeAndTurns[np.where(TimeAndTurns[:,2] == 2)[0]][:,3])
    avgSpikesW = np.mean(TimeAndTurns[np.where(TimeAndTurns[:,2] == 3)[0]][:,3])
    stSpikesW = np.std(TimeAndTurns[np.where(TimeAndTurns[:,2] == 3)[0]][:,3])
    
    egoLabels = ["Left", "Back", "Right"]
    alloLabels = ["North", "East", "South", "West"]
    egoAvgs = [avgSpikesL, avgSpikesB, avgSpikesR]
    alloAvgs = [avgSpikesN, avgSpikesE, avgSpikesS, avgSpikesW]
    egoSt = [stSpikesL, stSpikesB, stSpikesR]
    alloSt = [stSpikesN, stSpikesE, stSpikesS, stSpikesW]
    
    fig,ax = plt.subplots(1,2)
    ax[0].bar(egoLabels, egoAvgs, yerr=egoSt)
    ax[1].bar(alloLabels, alloAvgs, yerr=alloSt)
    ax[0].set_ylabel("Average spike rate (Hz)")
    fig.suptitle(suptitle)
    ax[0].set_title("Egocentric")
    ax[1].set_title("Allocentric")


#2D histogram of vel vs angle
"""
RList = RDP4(pos_cm, 1).ResultList
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))
directions = np.diff(RList[:,1:3], axis=0)
theta = angle(directions)  #point to point angle, not turn angle
hist,xe,ye = np.histogram2d(v[1:], theta, bins=30)
fig,ax = plt.subplots()
ax.imshow(hist, extent=[0,np.max(v),0,3.14])
ax.set_aspect(np.max(v)/3.14)
ax.set_xlabel("Instantaneous velocity (cm/s)")
ax.set_ylabel("Point to point angle (radians)")
"""


def polarRateAngle(idx, theta, Rpos, pos, spikes, title="", binsize=np.pi/36, RMrange=7*4.72):
    """
    Makes a polar plot of firing rate based on turn angle

    """
    
    rates = np.array([])
    for i in idx: #adapted from turn3
        j = np.where(pos == Rpos[i,0])[0] #1st point of the turn, in pos frame
        dist = 0
        #idx2: which pts are within range of 1st pt of turns, in pos frame
        while dist < RMrange: #checking points after the turn
            idx2_end = deepcopy(j)
            j += 1
            if j > len(pos)-1:
                break
            dist = np.sqrt((Rpos[i, 1]-pos[j, 1])**2 + (Rpos[i, 2]-pos[j, 2])**2)
        
        j = np.where(pos == Rpos[i,0])[0] - 1 #one point before the turn
        if j >= 0:
            dist = np.sqrt((Rpos[i, 1]-pos[j, 1])**2 + (Rpos[i, 2]-pos[j, 2])**2)
            idx2_start = j+1
        else:
            dist = RMrange+10
            idx2_start = 0
        while dist < RMrange: #checking points before the turn
            idx2_start = deepcopy(j)
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((Rpos[i, 1]-pos[j, 1])**2 + (Rpos[i, 2]-pos[j, 2])**2)
        
        start = bisect(spikes[:,0], pos[idx2_start,0])
        end = bisect(spikes[:,0], pos[idx2_end,0])
        rate = (end-start) / (pos[idx2_end,0]-pos[idx2_start,0]) * 1e6
        rates = np.hstack((rates,rate))
    
    #rates = np.array([])
    #for i in range(len(idx)): #spike rate during a turn
    #    start = bisect(spikes[:,0], Rpos[idx[i],0])
    #    end = bisect(spikes[:,0], Rpos[idx2[i],0])
    #    rate = (end-start) / (Rpos[idx2[i],0]-Rpos[idx[1],0]) * 1e6
    #    rates = np.hstack((rates, rate))
    
    tbins = np.arange(0, 2*np.pi, binsize)
    bRates = np.array([])
    for i in tbins:
        bRate = np.mean(rates[np.logical_and(theta>=i, theta<(i+binsize))])
        bRates = np.hstack((bRates,bRate))
        #print(round(i/np.pi*180), round(bRate*1e10,1))
    
    ax = plt.subplot(111, projection="polar")
    ax.plot(tbins, bRates)
    ax.set_title(title)
    #ax.set_rticks([2e-3, 4e-3, 6e-3])
    ax.grid(True)
    