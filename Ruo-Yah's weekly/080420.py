# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 11:58:49 2020

@author: Ruo-Yah Lai
"""
import os
os.environ["PROJ_LIB"] = "C:\\Users\\Ruo-Yah Lai\\anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm
from copy import deepcopy
from bisect import bisect_left
from placeFieldBorders import reorderBorders
from organized import visits
from RateMap1D import turnsInField3


"""
#testing path and contains_points
points = np.random.rand(50,2) * 250 + np.array([250, 100])
points = np.vstack((points, [360, 220], [358, 216]))

field = path.Path(borders1[2])
inField = field.contains_points(points)

fig, ax = plt.subplots(1,1)
ax.plot(borders1[2][:,0], borders1[2][:,1])
ax.scatter(points[inField,0], points[inField,1], c="r", s=2, marker="o", label="Inside the border")
ax.scatter(points[~inField,0], points[~inField,1], c="g", s=2, marker="o", label="Outside the border")
ax.axis("equal")
fig.legend()
"""


#adapted from 072220
def polarRateAngle2(pos, spikes, timeAndThetas, title="", binsize=np.pi/4, RMrange=7*4.72):
    """
    Makes 2 conjunctive polar plots of firing rate based on allocentric
    and egocentric turn angles
    A hammer projection and a flat 2D plot where alpha is adjusted by the
    number of trajectories
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
    abins = np.arange(0, 2*np.pi, binsize)
    meanRates = np.zeros((int(1.5*np.pi/binsize), int(2*np.pi/binsize))) #number of ebins x number of abins
    Ns = np.empty((int(1.5*np.pi/binsize), int(2*np.pi/binsize)))
    for i, ebin in enumerate(ebins):
        for j, abin in enumerate(abins):
            inBin = (timeAndThetas[:,1]>=ebin) & (timeAndThetas[:,1]<(ebin+binsize)) \
                    & (timeAndThetas[:,2]>=abin) & (timeAndThetas[:,2]<(abin+binsize))
            meanRates[i,j] = np.nanmean(rates[inBin])
            Ns[i,j] = np.sum(inBin)
            

    abins = np.hstack((abins, abins[-1]+binsize))*180/np.pi
    ebins = np.hstack((ebins, ebins[-1]+binsize))*180/np.pi
    fig, axs = plt.subplots(2, 1, figsize=(8,8))
    vmax = np.nanpercentile(meanRates,95)
    axs[0].set_title(title + "\n" + f"95th percentile of firing rate = {round(vmax,1)}", y=1.1)
    
    #hammer projection
    m = Basemap(projection="hammer", lon_0 = 0, ax=axs[0])
    m.drawparallels(np.arange(-90, 90+binsize*180/np.pi*2/3, binsize*180/np.pi*2/3))
    m.drawmeridians(np.arange(0, 360+binsize*180/np.pi, binsize*180/np.pi))
    x, y = np.meshgrid(abins-180, ebins*2/3-120)
    cb = m.pcolormesh(x, y, meanRates, latlon=True, cmap="viridis", vmax=vmax)
    
    #flat, alpha of colors is number of trajectories
    plt.xticks(np.arange(0,405,45), ["0\nN", "45", "90\nE", "135", "180\nS",
                                     "225", "270\nW", "315", "360\nN"])
    plt.yticks(np.arange(45,360,45), ["45", "L 90", "135", "B 180", "225", "R 270", "315"])
    axs[1].set_xlim(0, 360)
    axs[1].set_ylim(45, 315)
    axs[1].set_xlabel("Allocentric turn direction (degrees)")
    axs[1].set_ylabel("Egocentric turn direction (degrees)")
    axs[1].set_title(f"90th percentile number of trajectories = {int(np.percentile(Ns,90))}")
    alphas = Ns/np.percentile(Ns,90) *0.8 + 0.2
    alphas[np.where(alphas > 1)] = 1
    meanRates2 = deepcopy(meanRates)
    meanRates2[np.where(meanRates != meanRates)] = 0
    colors = cm.viridis(meanRates2/vmax)
    colors[:,:,3] = alphas
    colors[np.where(meanRates != meanRates)] = 0
    for i in range(colors.shape[0]):
        for j in range(colors.shape[1]):
            patch = patches.Rectangle((abins[j],ebins[i]), binsize*180/np.pi, binsize*180/np.pi, color=colors[i,j,:])
            axs[1].add_patch(patch)
    
    axcb = fig.colorbar(cb)
    axcb.set_label("Rate (Hz)")
    
    
def polarRateAngle2_1(pos, spikes, timeAndThetas, title="", binsize=np.pi/4, RMrange=7*4.72):
    """
    Makes 2 conjunctive polar plots of firing rate based on allocentric
    and egocentric turn angles
    Separates turns going into and coming out of the field
    A hammer projection and a flat 2D plot where alpha is adjusted by the
    number of trajectories
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
    abins = np.arange(0, 2*np.pi, binsize)
    meanRates = np.zeros((3, int(1.5*np.pi/binsize), int(2*np.pi/binsize))) #3 x number of ebins x number of abins
    Ns = np.empty((3, int(1.5*np.pi/binsize), int(2*np.pi/binsize)))
    for i, ebin in enumerate(ebins):
        for j, abin in enumerate(abins):
            #into the field
            inBin = (timeAndThetas[:,1]>=ebin) & (timeAndThetas[:,1]<(ebin+binsize)) \
                    & (timeAndThetas[:,2]>=abin) & (timeAndThetas[:,2]<(abin+binsize)) \
                    & ((timeAndThetas[:,3] == 0) | (timeAndThetas[:,3] == 3))
            meanRates[0,i,j] = np.nanmean(rates[inBin])
            Ns[0,i,j] = np.sum(inBin)
            
            #inside the field
            inBin = (timeAndThetas[:,1]>=ebin) & (timeAndThetas[:,1]<(ebin+binsize)) \
                    & (timeAndThetas[:,2]>=abin) & (timeAndThetas[:,2]<(abin+binsize)) \
                    & (timeAndThetas[:,3] == 1)
            meanRates[1,i,j] = np.nanmean(rates[inBin])
            Ns[1,i,j] = np.sum(inBin)
            
            #out from the field
            inBin = (timeAndThetas[:,1]>=ebin) & (timeAndThetas[:,1]<(ebin+binsize)) \
                    & (timeAndThetas[:,2]>=abin) & (timeAndThetas[:,2]<(abin+binsize)) \
                    & ((timeAndThetas[:,3] == 2) | (timeAndThetas[:,3] == 4) | (timeAndThetas[:,3] == 5))
            meanRates[2,i,j] = np.nanmean(rates[inBin])
            Ns[2,i,j] = np.sum(inBin)
            

    abins = np.hstack((abins, abins[-1]+binsize))*180/np.pi
    ebins = np.hstack((ebins, ebins[-1]+binsize))*180/np.pi
    fig, axs = plt.subplots(2, 3, figsize=(16,8))
    vmax = np.nanpercentile(meanRates,95)
    axs[0,1].set_title(title + "\n" + f"95th percentile of firing rate = {round(vmax,1)}"
                       + "\n\nInside the field")
    
    #hammer projection
    for i in range(3):
        m = Basemap(projection="hammer", lon_0 = 0, ax=axs[0,i])
        m.drawparallels(np.arange(-90, 90+binsize*180/np.pi*2/3, binsize*180/np.pi*2/3))
        m.drawmeridians(np.arange(0, 360+binsize*180/np.pi, binsize*180/np.pi))
        x, y = np.meshgrid(abins-180, ebins*2/3-120)
        cb = m.pcolormesh(x, y, meanRates[i], latlon=True, cmap="viridis", vmax=vmax)
    axs[0,0].set_title("Into the field")
    #axs[0,1].set_title("Inside the field")
    axs[0,2].set_title("Out from the field")
    
    #flat, alpha of colors is number of trajectories
    for i in range(3):
        axs[1,i].set_xticks(np.arange(0,405,45))
        axs[1,i].set_xticklabels(["0\nN", "45", "90\nE", "135", "180\nS",
                                  "225", "270\nW", "315", "360\nN"])
        axs[1,i].set_yticks(np.arange(45,360,45))
        axs[1,i].set_yticklabels(["45", "L 90", "135", "B 180", "225", "R 270", "315"])
        axs[1,i].set_xlim(0, 360)
        axs[1,i].set_ylim(45, 315)
        axs[1,i].set_xlabel("Allocentric turn direction (degrees)")
        axs[1,i].set_ylabel("Egocentric turn direction (degrees)")
        percentile = np.percentile(Ns[i][np.where(Ns[i] != 0)], 90)
        axs[1,i].set_title(f"90th percentile \n # of trajectories = {int(percentile)}")
        alphas = Ns[i]/percentile *0.8 + 0.2
        alphas[np.where(alphas > 1)] = 1
        meanRates2 = deepcopy(meanRates[i])
        meanRates2[np.where(meanRates2 != meanRates2)] = 0
        colors = cm.viridis(meanRates2/vmax)
        colors[:,:,3] = alphas
        colors[np.where(meanRates[i] != meanRates[i])] = 0
        for j in range(colors.shape[0]):
            for k in range(colors.shape[1]):
                patch = patches.Rectangle((abins[k],ebins[j]), binsize*180/np.pi, binsize*180/np.pi, color=colors[j,k,:])
                axs[1,i].add_patch(patch)
    
    axcb = fig.colorbar(cb)
    axcb.set_label("Rate (Hz)")
    fig.tight_layout()


def bulkGraphs(pos,unit1,unit2,unit5,unit6,unit7,unit8,timestamp):
    """
    Makes multiple graphs of conjunctive polar plots
    timestamp: current timestamp
    """
    A = [[unit1,0], [unit2,0], [unit2,3], [unit5,1],
         [unit5,2], [unit5,3], [unit6,0], [unit6,1], 
         [unit6,3], [unit7,0], [unit7,1], [unit8,3]]
    titles = ["1.1 subfield 0", "1.2 subfield 0", "1.2 subfield 3",
              "1.5 subfield 1", "1.5 subfield 2", "1.5 subfield 3", 
              "1.6 subfield 0", "1.6 subfield 1", "1.6 subfield 3",
              "1.7 subfield 0", "1.7 subfield 1","1.8 subfield 3"]
    for i,a in enumerate(A):
        borders = reorderBorders(a[0])
        visits_ = visits(borders, pos)
        tath = turnsInField3(visits_, pos, a[1])
        polarRateAngle2_1(pos, a[0].spikes, tath, "R859 D3 T6 "+titles[i])
        plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                    + timestamp + " - " + titles[i] + ".png")
        



"""
m = Basemap(projection="hammer", lon_0 = 0)
m.drawparallels(np.arange(-90,120,30))
m.drawmeridians(np.arange(0,420,60))
x, y = np.arange(-180,210,30), np.arange(-90,120,30)
x, y = np.meshgrid(x, y)
#z = np.arange(36).reshape(6,6)
z = np.random.rand(6,12)
#x, y = m(x,y)
m.pcolormesh(x, y, z, vmax=z.max(), latlon=True)
"""