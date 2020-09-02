# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 21:04:41 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import csv
from bisect import bisect_left
import matplotlib.pyplot as plt
from organized import velocity_filtering
from RateMap import takeClosest


#for R781 D4
def read_pos(path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R781 D4/", to_cm = False, file="pos.p.ascii"):
    if to_cm:
        a = np.array([1, 4.72, 4.72])
        ptsCm = 1
    else:
        a = np.array([1, 1, 1])
        ptsCm = 4.72
    with open(path + file, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        pos = [data for data in data_iter]
    with open(path+"sessionEpochInfo.txt", "r") as file:
        epochs = file.readlines()
    pos = np.array(pos[25:], dtype = "float64")
    start = bisect_left(pos[:, 0], float(epochs[0][:13]))
    end = bisect_left(pos[:, 0], float(epochs[0][14:27]))
    pos = pos[start:end, 0:3]/a[None,:]
    pos = pos[np.all(pos > np.array([0, 0, 0]), axis=1)]
    return velocity_filtering(pos, 3*ptsCm)

#for R781 D4
def getPosFromTs(yourts,position, path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R781 D4/"):
    '''position is the full list of ts, yourts are the times you'd like a pos for '''
    adjTs = [takeClosest(position[:,0],i) for i in yourts]
    target = np.empty((0,3))
    for i in adjTs:
        target = np.vstack((target,position[np.where(position[:,0]==i)][:,0:][0]))
    with open(path+"sessionEpochInfo.txt", "r") as file:
        epochs = file.readlines()
    start = float(epochs[0][:13])
    end = float(epochs[0][14:27])
    target = target[np.where(target[:,0] >= start) and target[:,0] < end]
    return target

def graph(title="", perim=[]):
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    colors = ["b","g","r","k","c","m","y"]
    for i in range(len(perim)):
        ax.plot(perim[i][:,0]/4.72, perim[i][:,1]/4.72, color=colors[i], 
                label=f"Subfield {i}")
    ax.legend()
