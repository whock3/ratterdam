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


#adapted from 061820
def polarRateAngle(pos, spikes, timeAndThetas, elim, alim, eticks, aticks, binsize=np.pi/20, RMrange=7*4.72):
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
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    ax.set_title("Egocentric", y=1.14)
    ax.set_ylim(0,elim)
    ax.set_rticks(eticks)
    ego = np.hstack((ebins.reshape((-1,1)), eRates.reshape((-1,1))))
    ego = ego.reshape((-1,1,2))
    eSegments = np.concatenate([ego[0:-1], ego[1:]], axis=1)
    lc = LineCollection(eSegments)
    lc.set_array(eNs)
    ax.add_collection(lc)
    axcb = fig.colorbar(lc)
    axcb.set_label("Number of trajectories")
    fig.tight_layout()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    ax.set_title("Allocentric", y=1.14)
    ax.set_ylim(0,alim)
    ax.set_rticks(aticks)
    allo = np.hstack((abins.reshape((-1,1)), aRates.reshape((-1,1))))
    allo = allo.reshape((-1,1,2))
    aSegments = np.concatenate([allo[0:-1], allo[1:]], axis=1)
    lc = LineCollection(aSegments)
    lc.set_array(aNs)
    ax.add_collection(lc)
    axcb = fig.colorbar(lc)
    axcb.set_label("Number of trajectories")
    fig.tight_layout()


def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)


def thetaDist(pos):
    """
    Calculates theta/s and dist/s between pos points
    """
    directions = np.diff(pos[:,1:3], axis=0)
    theta2 = angle2(directions) #[0, 2pi]
    theta = deepcopy(theta2)
    for i in range(len(theta2)):
        if theta[i] > np.pi:
            theta[i] = theta2[i] - 2*np.pi #[-pi, pi]
    theta = theta / np.diff(pos[:,0])[1:] *1e6 /np.pi*180
    dist = np.linalg.norm(directions, axis=1)[1:] / np.diff(pos[:,0])[1:] *1e6 /4.72
    
    theta1 = np.diff(pos[1:,3]) / np.diff(pos[:,0])[1:] *1e6
    
    ts = np.reshape(pos[1:-1,0], (-1,1))
    theta = theta.reshape((-1,1))
    dist = dist.reshape((-1,1))
    theta1 = theta1.reshape((-1,1))
    return np.hstack((ts, theta1, dist))


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
#includes head direction
def read_pos2(path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/", pos="pos.p.ascii", to_cm = False):
    data_folder = Path(path)
    file_to_open = data_folder / pos
    if to_cm:
        a = np.array([1, 4.72, 4.72, 1])
        ptsCm = 1
    else:
        a = np.array([1,1,1,1])
        ptsCm = 4.72
    with open(file_to_open, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    data_np1 = np.array(data[50574:-17], dtype = "float64")
    data_np2 = data_np1[:,0:]/a[None,:]
    data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0, 0]), axis=1)]
    return velocity_filtering(data_np3, 3*ptsCm)


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