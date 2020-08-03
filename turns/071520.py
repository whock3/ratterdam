# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 13:56:30 2020

@author: Ruo-Yah Lai
"""

import numpy as np
from pathlib import Path
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from organized import makeRM2, weird_smooth
import matplotlib.cm as cm


def noVelFilter(path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/", pos="pos.p.ascii", to_cm = True):
    """
    read_pos without velocity filtering
    """
    data_folder = Path(path)
    file_to_open = data_folder / pos
    if to_cm == True:
        a = np.array([1, 4.72, 4.72])
    else:
        a = np.array([1,1,1])
    with open(file_to_open, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    data_np1 = np.array(data[50574:-17], dtype = "float64")
    data_np2 = data_np1[:,0:3]/a[None,:]
    data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
    return data_np3

"""
#average time difference between the 1st and last points of a turn
RList = RDP4(pos, 4.72).ResultList
idx3,_,idx3_2,_ = turn3(RList, np.pi/12, np.pi/4)
a = RList[idx3, 0]-RList[idx3_2, 0]
b = deepcopy(a)
b[np.where(a == 0)] = np.nan
"""


def graphRM3D(position, pos2, spikes, suptitle, percentile=99, mycm="jet", smoothing_2d_sigma=1):
    """ 
    Avg RM for turns, showing number of trajectories
    position, pos2, and spikes from shiftPosP and shiftPosS
    """
    rows, cols = np.linspace(-7*4.72, 7*4.72, 15), np.linspace(-7*4.72, 7*4.72, 15)
    n = []
    ho = []
    titles2 = []
    X,Y = np.mgrid[-6.5:6.5:14j, -6.5:6.5:14j]
    for i in range(8):
        n1 = []
        ho1 = np.histogram2d(pos2[i][:,2], pos2[i][:,1], bins=[rows,cols])[0]
        for j in range(len(position[i])):
            n1.append(makeRM2(spikes[i][j], position[i][j]))
        n2 = np.nanmean(n1, axis=0)
        n2 = weird_smooth(n2,smoothing_2d_sigma)
        n2[np.where(ho1==0)] = np.nan
        n.append(n2)
        ho.append(ho1)
        titles2.append(f"n = {len(n1)}")
    
    fig = plt.figure(figsize=(16,8))
    zmax = np.nanpercentile(n[0], percentile)
    fig.suptitle(suptitle + "\n" + f"Cutoff = {percentile}th percentile, {round(zmax,1)} Hz", y=1.08)
    titles = ["All turns", "Left", "Back", "Right", "North", "East", "South", "West"]
    
    #normalize ho
    hoMax = np.empty(8) 
    for i in range(8):
        hoMax[i] = np.nanpercentile(ho[i], percentile)
        ho[i] = ho[i]/hoMax[i]
    
    for i in range(8):
        ax = fig.add_subplot(2, 4, i+1, projection="3d")
        ax.set_title(titles[i] + "\n" + titles2[i] + "\nCutoff = " + f"{round(hoMax[i],1)}")
        ax.plot_surface(X,Y,n[i],facecolors=cm.jet(ho[i]), vmin=0, vmax=hoMax[i])
        ax.set_zlim(0, zmax)\
    #sm = plt.cm.ScalarMappable(cmap="jet")#, norm=[0, hoMax])
    #cb = fig.colorbar(sm)
    #cb.set_label("Rate (Hz)")
    fig.tight_layout()

def sigmoid(x):
    return 100 / (1 + np.exp(-0.1*(x - 50)))

def graphRM2Dt(position, pos2, spikes, suptitle="", percentile=95, mycm="jet", smoothing_2d_sigma=1):
    """ 
    Avg RM for turns, divided by the number of trajectories in each bin
    position, pos2, and spikes from shiftPosP and shiftPosS
    """
    rows, cols = np.linspace(-7*4.72, 7*4.72, 15), np.linspace(-7*4.72, 7*4.72, 15)
    n = []
    titles2 = []
    for i in range(8):
        n1s = []
        ho1s = []
        ho = np.histogram2d(pos2[i][:,2], pos2[i][:,1], bins=[rows,cols])[0]
        for j in range(len(position[i])):
            n1, ho1 = makeRM2(spikes[i][j], position[i][j])
            n1s.append(n1)
            ho1[np.where(ho1 > 1)] = 1
            ho1[np.where(ho1 != ho1)] = 0
            ho1s.append(ho1)
            
        n2 = np.nanmean(n1s, axis=0)
        n2 = weird_smooth(n2,smoothing_2d_sigma)
        n2 = n2/sigmoid(np.sum(ho1s, axis=0))
        n2[np.where(ho==0)] = np.nan
        n.append(n2)
        titles2.append(f"n = {len(n1s)}")
    
    fig, axs = plt.subplots(2,4)
    vmax = np.nanpercentile(n[0], percentile)
    fig.suptitle(suptitle + "\nNormalized by the number of trajectories per bin\n" 
                 + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz", y=1.12)
    titles = ["All turns", "Left", "Back", "Right", "North", "East", "South", "West"]

    ims = []
    for i in range(2):
        for j in range(4):
            axs[i][j].set_title(titles[i*4+j] + "\n" + titles2[i*4+j])
            ims.append(axs[i][j].imshow(n[i*4+j], cmap=mycm, origin="lower", vmin=0, vmax=vmax, extent=(-7,7,-7,7)))
            axs[i][j].axis("equal")
    cb = fig.colorbar(ims[0])
    cb.set_label("Rate (Hz)")
    fig.tight_layout()
    

def stack(srPs, srSs, srP2s=np.empty(0)):
    A = [[] for _ in range(len(srPs[0]))]
    B = [[] for _ in range(len(srSs[0]))]
    if len(srP2s) != 0:
        C = [np.empty((0,3)) for _ in range(len(srP2s[0]))]
    else:
        C = []
    for i in range(len(srPs)):
        for j in range(len(srPs[0])):
            if len(C) != 0:
                C[j] = np.vstack((C[i], srP2s[i][j]))
            for k in range(len(srPs[i][j])):
                A[j].append(srPs[i][j][k])
                B[j].append(srSs[i][j][k])
    return A, B, C