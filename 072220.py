# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from bisect import bisect_left
from matplotlib.collections import LineCollection
from RateMap import weird_smooth
from organized import velocity_filtering, makeRM2
from path import Path
import csv
from RateMap import getPosFromTs


#adapted from 061820
def polarRateAngle(pos, spikes, timeAndThetas, elim, alim, eticks=[], aticks=[], binsize=np.pi/15, RMrange=7*4.72):
    """
    Makes a polar plot of firing rate based on turn angle
    Colored by number of trajectories
    """
    rates = np.array([])
    for i in range(len(timeAndThetas)):
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
        
        start = bisect_left(spikes[:,0], pos[idx_start,0])
        end = bisect_left(spikes[:,0], pos[idx_end,0])
        rate = (end-start) / (pos[idx_end,0]-pos[idx_start,0]) * 1e6
        rates = np.hstack((rates,rate))
        
    ebins = np.arange(np.pi/4, 7/4*np.pi, binsize)
    eRates = np.array([])
    eNs = np.empty(0)
    for i in ebins:
        eRate = np.mean(rates[np.logical_and(timeAndThetas[:,1]>=i, 
                                             timeAndThetas[:,1]<(i+binsize))])
        eRates = np.hstack((eRates,eRate))
        eNs = np.hstack((eNs, np.sum(np.logical_and(timeAndThetas[:,1]>=i, 
                                                    timeAndThetas[:,1]<(i+binsize)))))
        if eNs[-1] == 0:
            print(f"0 trajectory in bin {round(i,2)}")
    
    abins = np.arange(0, 2*np.pi, binsize)
    aRates = np.array([])
    aNs = np.empty(0)
    for i in abins:
        aRate = np.mean(rates[np.logical_and(timeAndThetas[:,2]>=i, 
                                             timeAndThetas[:,2]<(i+binsize))])
        aRates = np.hstack((aRates,aRate))
        aNs = np.hstack((aNs, np.sum(np.logical_and(timeAndThetas[:,2]>=i, 
                                                    timeAndThetas[:,2]<(i+binsize)))))
        if aNs[-1] == 0:
            print(f"0 trajectory in bin {round(i,2)}")
    
    titles = ["Egocentric", "Allocentric"]
    lims = [elim, alim]
    ticks = [eticks, aticks]
    bins = [ebins+binsize/2, abins+binsize/2]
    Rates = [eRates, aRates]
    Ns = [eNs, aNs]
    for i in range(2):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        ax.set_title(titles[i], y=1.14)
        ax.set_ylim(0, lims[i])
        ax.set_rticks(ticks[i])
        ax.plot(bins[i], Rates[i], color="k")
        
        segments = np.zeros((len(bins[i]),1,2))
        segments = np.concatenate([segments[0:-1], segments[1:]], axis=1)
        lc = LineCollection(segments)
        lc.set_array(Ns[i])
        
        #circle = np.full(len(cbins[i]), lims[i]*0.9)
        #circle = np.hstack((cbins[i].reshape((-1,1)), circle.reshape((-1,1))))
        #circle = circle.reshape((-1,1,2))
        #segments = np.concatenate([circle[:-1], circle[1:]], axis=1)
        #lc = LineCollection(segments)
        #lc.set_array(Ns[i])
        #ax.add_collection(lc)
        
        colors = plt.cm.viridis(Ns[i] / max(Ns[i]))
        ax.bar(bins[i], lims[i], width=binsize, color=colors, alpha=0.4)
        axcb = fig.colorbar(lc)
        axcb.set_label("Number of trajectories")
        fig.tight_layout()
        

def thetaDist(pos):
    """
    Calculates theta/s and dist/s between pos points
    """
    directions = np.diff(pos[:,1:3], axis=0)
    dist = np.linalg.norm(directions, axis=1)[1:] / np.diff(pos[:,0])[1:] *1e6 /4.72
    
    theta = np.diff(pos[1:,3])
    negatives = np.where(theta < -180)
    positives = np.where(theta > 180)
    theta[negatives] = theta[negatives] + 360
    theta[positives] = theta[positives] - 360 #[-180, 180]
    theta = theta / np.diff(pos[:,0])[1:] *1e6
    
    ts = np.reshape(pos[1:-1,0], (-1,1))
    theta = theta.reshape((-1,1))
    dist = dist.reshape((-1,1))
    return np.hstack((ts, theta, dist))


#makeRM from RateMap modified
def makeRM3(spikes, pos, bins = [50,50], smoothing_2d_sigma=2):
    """
    For speed vs change in theta
    """
    rmax = np.percentile(pos[:,2], 99)
    cmax = max(abs(np.percentile(pos[:,1], 1)), np.percentile(pos[:,1], 99))
    rows = np.linspace(0, rmax, bins[0])
    cols = np.linspace(-cmax, cmax, bins[1])
    hs,xe,ye = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rows, cols])
    ho = np.histogram2d(pos[:,2],pos[:,1],bins=[rows, cols])[0]
    n = (hs*np.reciprocal(ho))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(ho==0)] = np.nan
    n = weird_smooth(n,smoothing_2d_sigma)
    n[np.where(ho==0)] = np.nan
    return n, (int(cmax), int(rmax))


#graphRM from organized modified
def graphRM3(position, spikes, title, percentile=99, mycm="jet", sigma=2):
    """
    Graphs speed vs change in theta
    """
    n, minmax = makeRM3(spikes, position, smoothing_2d_sigma=sigma)
    vmax = np.nanpercentile(n, percentile)
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title + "\n" + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    im = ax.imshow(n, cmap=mycm, origin="lower", vmin=0, vmax=vmax)#, extent=minmax)
    ax.set_xlabel("∆θ/sec (deg/sec)")
    ax.set_ylabel("∆S/sec (cm/sec)")
    #ax.set_xlim(minmax[0], minmax[1])
    #ax.set_ylim(minmax[2], minmax[3])
    plt.xticks([0, 12.5, 25, 37.5, 50],
               [-minmax[0], int(-minmax[0]/2), 0, int(minmax[0]/2), minmax[0]])
    plt.yticks([0, 12.5, 25, 37.5, 50],
               [0, int(minmax[1]/4), int(minmax[1]/2), int(minmax[1]/4*3), minmax[1]])
    #ax.legend()
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")
    

#read_pos from organized modified
def read_pos2(path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/", file="pos.p.ascii", to_cm = False):
    """
    includes head direction
    """
    if to_cm:
        a = np.array([1, 4.72, 4.72, 1])
        ptsCm = 1
    else:
        a = np.array([1,1,1,1])
        ptsCm = 4.72
    with open(path + file, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        pos = [data for data in data_iter]
    with open(path+"sessionEpochInfo.txt", "r") as file:
        epochs = file.readlines()
    pos = np.array(pos[25:], dtype = "float64")
    start = bisect_left(pos[:, 0], float(epochs[0][:10]))
    end = bisect_left(pos[:, 0], float(epochs[0][11:21]))
    pos = pos[start:end]/a[None,:]
    pos = pos[np.all(pos > np.array([0, 0, 0, 0]), axis=1)]
    return velocity_filtering(pos, 3*ptsCm)


#adapted from turnsInField from organized
def turnsInField4(visits, position, subfield):
    """
    Finds ∆θ and ∆S near a subfield
    position: 4 columns, before RDP
    Returns [ts, theta, dist]
    """
    ranges = np.empty((0,2))
    for i in range(len(visits[subfield])): #the ith visit
        #start and end are indexes of the start and end of visit
        start = bisect_left(position[:,0], visits[subfield][i][0]-400000)
        end = bisect_left(position[:,0], visits[subfield][i][-1]+400000)
        if end-start > 0:
            ranges = np.vstack((ranges, np.reshape([start, end],(1,2))))
            
    postd = np.empty((0, 3))
    for i in range(len(ranges)): #the ith visit
        postd1 = thetaDist(position[int(ranges[i,0]):int(ranges[i,1])])
        postd = np.vstack((postd, postd1))
    return postd


#shiftPos from organized modified for conjunctive turns, no rotation
def shiftPos2(pos, spikes, timeAndTurns, RMrange=7*4.72):
    """
    12 combinations of 3 ego and 4 allo turn directions
    pos: [ts,x,y] before RDP
    timeAndTurns: [ts, ego turn direction, allo turn direction] from classifyTurns
    
    Returns: 
        shifted_pos: list of 12 lists each with arrays where each array is 1 turn
    pos (before RDP) within RMrange of the 1st point of turns, 
    shifted based on the 1st point of turns
        shifted_pos2: list of 12 arrays
        shifted_s: list of 12 lists of spikes
    """
    shifted_pos = [[] for _ in range(12)]
    shifted_pos2 = [np.empty((0,3)) for _ in range(12)]
    shifted_s = [[] for _ in range(12)]
    for i in range(len(timeAndTurns)): #adapted from turn3
        k = np.where(pos == timeAndTurns[i,0])[0] #1st point of the turn, in pos frame
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
        
        shiftedP = np.subtract(pos[int(idx_start):int(idx_end+1)], np.array([0,pos[k,1],pos[k,2]]))       
        if idx_end < len(pos)-1:
            spike_idx = np.where(np.logical_and(spikes>pos[int(idx_start),0], spikes<pos[int(idx_end+1),0]))[0]
        else:
            spike_idx = np.where(np.logical_and(spikes>pos[int(idx_start),0], spikes<pos[int(idx_end),0]))[0]
        shiftedS = np.subtract(spikes[spike_idx], np.array([0,pos[k,1],pos[k,2]]))
        
        
        if timeAndTurns[i,1] == 0:
            if timeAndTurns[i,2] == 0:
                shifted_pos[0].append(shiftedP)
                shifted_pos2[0] = np.vstack((shifted_pos2[0],shiftedP))
                shifted_s[0].append(shiftedS)
            elif timeAndTurns[i,2] == 1:
                shifted_pos[1].append(shiftedP)
                shifted_pos2[1] = np.vstack((shifted_pos2[1],shiftedP))
                shifted_s[1].append(shiftedS)
            elif timeAndTurns[i,2] == 2:
                shifted_pos[2].append(shiftedP)
                shifted_pos2[2] = np.vstack((shifted_pos2[2],shiftedP))
                shifted_s[2].append(shiftedS)
            elif timeAndTurns[i,2] == 3:
                shifted_pos[3].append(shiftedP)
                shifted_pos2[3] = np.vstack((shifted_pos2[3],shiftedP))
                shifted_s[3].append(shiftedS)
        elif timeAndTurns[i,1] == 1:
            if timeAndTurns[i,2] == 0:
                shifted_pos[4].append(shiftedP)
                shifted_pos2[4] = np.vstack((shifted_pos2[4],shiftedP))
                shifted_s[4].append(shiftedS)
            elif timeAndTurns[i,2] == 1:
                shifted_pos[5].append(shiftedP)
                shifted_pos2[5] = np.vstack((shifted_pos2[5],shiftedP))
                shifted_s[5].append(shiftedS)
            elif timeAndTurns[i,2] == 2:
                shifted_pos[6].append(shiftedP)
                shifted_pos2[6] = np.vstack((shifted_pos2[6],shiftedP))
                shifted_s[6].append(shiftedS)
            elif timeAndTurns[i,2] == 3:
                shifted_pos[7].append(shiftedP)
                shifted_pos2[7] = np.vstack((shifted_pos2[7],shiftedP))
                shifted_s[7].append(shiftedS)
        elif timeAndTurns[i,1] == 2:
            if timeAndTurns[i,2] == 0:
                shifted_pos[8].append(shiftedP)
                shifted_pos2[8] = np.vstack((shifted_pos2[8],shiftedP))
                shifted_s[8].append(shiftedS)
            elif timeAndTurns[i,2] == 1:
                shifted_pos[9].append(shiftedP)
                shifted_pos2[9] = np.vstack((shifted_pos2[9],shiftedP))
                shifted_s[9].append(shiftedS)
            elif timeAndTurns[i,2] == 2:
                shifted_pos[10].append(shiftedP)
                shifted_pos2[10] = np.vstack((shifted_pos2[10],shiftedP))
                shifted_s[10].append(shiftedS)
            elif timeAndTurns[i,2] == 3:
                shifted_pos[11].append(shiftedP)
                shifted_pos2[11] = np.vstack((shifted_pos2[11],shiftedP))
                shifted_s[11].append(shiftedS)
        
    return shifted_pos, shifted_pos2, shifted_s


#graphRM2 from organized modified for 12 graphs
def graphRM2_1(position, pos2, spikes, suptitle, percentile=95, mycm="jet", smoothing_2d_sigma=1):
    """ 
    Conjuntive graphs of avg RM for turns
    position, pos2, and spikes from shiftPos2
    """
    rows, cols = np.linspace(-7*4.72, 7*4.72, 15), np.linspace(-7*4.72, 7*4.72, 15)
    n = []
    titles2 = []
    for i in range(12):
        n1 = []
        ho = np.histogram2d(pos2[i][:,2], pos2[i][:,1], bins=[rows,cols])[0]
        for j in range(len(position[i])):
            n1.append(makeRM2(spikes[i][j], position[i][j])[0])
        n2 = np.nanmean(n1, axis=0)
        if len(n1) == 0:
            n2 = np.empty((14,14))
            n2[:] = np.nan
        n2 = weird_smooth(n2,smoothing_2d_sigma)
        n2[np.where(ho==0)] = np.nan
        n.append(n2)
        titles2.append(f"n = {len(n1)}")
    
    fig, axs = plt.subplots(3,4, figsize=(8,7))
    vmax = np.nanpercentile(n, percentile)
    fig.suptitle(suptitle + "\n" + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz", y=1.05)
    titles = ["Left & North", "Left & East", "Left & South", "Left & West",
              "Back & North", "Back & East", "Back & South", "Back & West", 
              "Right & North", "Right & East", "Right & South", "Right & West"]

    ims = []
    for i in range(3):
        for j in range(4):
            axs[i][j].set_title(titles[i*4+j] + "\n" + titles2[i*4+j])
            ims.append(axs[i][j].imshow(n[i*4+j], cmap=mycm, origin="lower", vmin=0, vmax=vmax, extent=(-7,7,-7,7)))
            axs[i][j].axis("equal")
    cb = fig.colorbar(ims[0])
    cb.set_label("Rate (Hz)")
    fig.tight_layout()
    

#trajectories from organized modified for 12 graphs
def trajectories_1(position, suptitle=""):
    """
    Graph trajectories after shifting them
    For conjuntive graphs of avg RM for turns
    """
    fig,axs = plt.subplots(3,4, figsize=(8,7))
    titles = ["Left & North", "Left & East", "Left & South", "Left & West",
              "Back & North", "Back & East", "Back & South", "Back & West", 
              "Right & North", "Right & East", "Right & South", "Right & West"]
    fig.suptitle(suptitle, y=1.02)
        
    for i in range(3):
        for j in range(4):
            axs[i][j].set_title(titles[i*4+j] + "\n" + f"n = {len(position[i*4+j])}")
            for k in range(len(position[i*4+j])):
                axs[i][j].plot(position[i*4+j][k][:,1]/4.72, position[i*4+j][k][:,2]/4.72)
                axs[i][j].axis("equal")
    fig.tight_layout()


def bulkGraphs(pos,unit1,unit2,unit4,unit5,unit6,spikes1,spikes2,spikes4,spikes5,spikes6,timestamp):
    """
    Makes multiple graphs of speed vs change in theta
    pos: 4 columns, from read_pos2
    timestamp: current timestamp
    """
    A = [[unit1,0,spikes1], [unit2,0,spikes2], [unit2,3,spikes2], 
         [unit4,0,spikes4], [unit5,0,spikes5], [unit5,1,spikes5],
         [unit5,2,spikes5], [unit5,3,spikes5], [unit6,0,spikes6], 
         [unit6,1,spikes6], [unit6,2,spikes6], [unit6,3,spikes6]]
    titles = ["1.1 subfield 0", "1.2 subfield 0", "1.2 subfield 3",
              "1.4 subfield 0", "1.5 subfield 0", "1.5 subfield 1",
              "1.5 subfield 2", "1.5 subfield 3", "1.6 subfield 0",
              "1.6 subfield 1", "1.6 subfield 2", "1.6 subfield 3"]
    for i,a in enumerate(A):
        postd = turnsInField4(a[0].visits, pos, a[1])
        spikes = getPosFromTs(a[2][:,0], postd)
        graphRM3(postd, spikes, "R859 D3 T6 "+titles[i])
        plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                    + timestamp + " - " + titles[i] + ".png")