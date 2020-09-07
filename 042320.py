# -*- coding: utf-8 -*-
"""
Angle function by user unutbu from:
https://stackoverflow.com/questions/14631776/calculate-turning-points-pivot-points-in-trajectory-path

"""
from velocity_filtering import velocity_filtering
from RDP2_class import RDP2
import numpy as np
import matplotlib.pyplot as plt
import csv

def angle(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return np.arccos((dir1*dir2).sum(axis=1)/(
        np.sqrt((dir1**2).sum(axis=1)*(dir2**2).sum(axis=1))))

with open("pos.p.ascii", 'r') as csv_file:
    data_iter = csv.reader(csv_file)
    data = [data for data in data_iter]
data_np1 = np.array(data[50574:-17], dtype = "float64")
a = np.array([1, 4.72, 4.72])
data_np2 = data_np1[:,0:3]/a[None,:]
data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
data_np4 = velocity_filtering(data_np3)


#comparing epsilons
b = 28000
epsilons = [0.5, 0.75, 1]
min_angle = np.pi*0.25
fig, ax = plt.subplots(1, 4, figsize=(24, 6))
ax[0].plot(data_np4[b:b+300, 1], data_np4[b:b+300, 2])
ax[0].set_title("Original")
ax[0].scatter(data_np4[b,1], data_np4[b,2], marker = "+", color = "r", label = "first")
ax[0].scatter(data_np4[b+299,1], data_np4[b+299,2], marker = "x", color = "r", label = "last")
for i in range(len(epsilons)):
    RDP2.ResultList = data_np4[b, 1:3]
    RList = RDP2.RDPfunc(RDP2, data_np4[b:b+300, 1:3], epsilons[i])
    directions = np.diff(RList, axis=0)
    theta = angle(directions)
    idx = np.where(theta>min_angle)[0]+1
    percent = round(100 - len(RList)/3, 2)
    ax[i+1].plot(RList[:,0], RList[:,1])
    ax[i+1].set_title(f"epsilon = {round(epsilons[i], 2)}" "\n" f"{percent}% pts thrown away")
    ax[i+1].scatter(RList[0,0], RList[0,1], marker = "+", color = "r", label = "first")
    ax[i+1].scatter(RList[-1,0], RList[-1,1], marker = "x", color = "r", label = "last")
    ax[i+1].scatter(RList[idx,0], RList[idx,1], s = 16, color = "r", label='turning points')
plt.legend()
plt.show()

"""
b = 20000
epsilons = [0.5, 0.75, 1]
min_angle = np.pi*0.25
fig, ax = plt.subplots(1, 4, figsize=(24, 6))
ax[0].plot(data_np4[b:b+1800, 1], data_np4[b:b+1800, 2])
ax[0].set_title("Original")
ax[0].scatter(data_np4[b,1], data_np4[b,2], marker = "+", color = "r", label = "first")
ax[0].scatter(data_np4[b+1799,1], data_np4[b+1799,2], marker = "x", color = "r", label = "last")
for i in range(len(epsilons)):
    RDP2.ResultList = data_np4[b, 1:3]
    RList = RDP2.RDPfunc(RDP2, data_np4[b:b+300, 1:3], epsilons[i])
    
    for j in range(1, 6):
        RDP2.ResultList = data_np4[b+j*300, 1:3]
        RList = np.vstack([RList, RDP2.RDPfunc(RDP2, data_np4[b+j*300:b+j*300+300, 1:3], epsilons[i])])
        
    directions = np.diff(RList, axis=0)
    theta = angle(directions)
    idx = np.where(theta>min_angle)[0]+1
    percent = 100 - len(RList)/18
    ax[i+1].plot(RList[:,0], RList[:,1])
    ax[i+1].set_title("epsilon = %2f \n%f percent pts thrown away" % (epsilons[i], percent))
    ax[i+1].scatter(RList[0,0], RList[0,1], marker = "+", color = "r", label = "first")
    ax[i+1].scatter(RList[-1,0], RList[-1,1], marker = "x", color = "r", label = "last")
    ax[i+1].scatter(RList[idx], RList[idx], 'ro', markersize = 10, label='turning points')

plt.legend()
plt.show()
"""
