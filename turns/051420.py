# -*- coding: utf-8 -*-
"""

"""

from velocity_filtering import velocity_filtering
from RDP3_class import RDP3
from matplotlib.collections import LineCollection
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import csv
"""
def angle(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return np.arccos((dir1*dir2).sum(axis=1)/(
        np.sqrt((dir1**2).sum(axis=1)*(dir2**2).sum(axis=1))))
"""
def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)

with open("pos.p.ascii", 'r') as csv_file:
    data_iter = csv.reader(csv_file)
    data = [data for data in data_iter]
data_np1 = np.array(data[50574:-17], dtype = "float64")
a = np.array([1, 4.72, 4.72])
data_np2 = data_np1[:,0:3]/a[None,:]
data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
data_np4 = velocity_filtering(data_np3)


#sum of "turns" within a certain range
b = 48600
epsilons = [0.5, 0.75, 1]
i = 2
min_angle = np.pi*0.25
turn_range = 5 #angles within x cm of each "turn" shouldn't add up to be a circle

RDP3.ResultList = data_np4[b]
RList = RDP3.RDPfunc(RDP3, data_np4[b:b+300], epsilons[i]) 
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

fig, ax = plt.subplots(1, 1, figsize = (7, 5))
ax.set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax.set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax.scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax.scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")

directions = np.diff(RList[:, 1:3], axis=0)
theta2 = angle2(directions) #[0, 2pi]
theta = deepcopy(theta2)
for i in range(len(theta2)):
    if theta[i] > np.pi:
        theta[i] = 2*np.pi - theta2[i] #[0, pi]
idx = np.where(theta>min_angle)[0]+1 #in RList frame

theta_sum = np.empty(len(idx))
for i in range(len(idx)):
    dist = np.sqrt((RList[idx[i], 1]-RList[:, 1])**2 + (RList[idx[i], 2]-RList[:, 2])**2)
    idx2 = np.where(dist < turn_range)[0]-1 #which pts are within range of "turns", in theta frame
    
    idx2 = idx2[idx2 != -1]
    idx2 = idx2[idx2 != (len(RList)-2)]
    idx2_consec1 = idx2[1:] - idx2[:-1]
    idx2_consec2 = np.where(idx2_consec1 != 1)[0]
    if len(idx2_consec2) > 0:
        idx2 = idx2[:idx2_consec2[0]+1]
        
    theta_sum[i] = np.sum(theta2[idx2]) % (2*np.pi)
    if theta_sum[i] > np.pi:
        theta_sum[i] = 2*np.pi - theta_sum[i] #[0, pi]
idx3 = np.where(theta_sum>min_angle)
idx3 = idx[idx3]
ax.scatter(RList[idx3,1], RList[idx3,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax.add_collection(lc)
axcb = fig.colorbar(lc)
axcb.set_label('Instantaneous velocity')
ax.set_title("velocity threshold = 2 \n epsilon = 1")

plt.legend() #bbox_to_anchor=(0.05, 0.45), loc="center left"
plt.show()


#whole session theta > 1.5, velocity < 35
epsilon = 1

diffts, diffx, diffy = np.diff(data_np4[:,0]), np.diff(data_np4[:,1]), np.diff(data_np4[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v = np.sqrt((vx**2)+(vy**2))
v = v[1:]

directions = np.diff(data_np4[:, 1:3], axis=0)
theta2 = angle2(directions) #[0, 2pi]
theta = deepcopy(theta2)
for i in range(len(theta2)):
    if theta[i] > np.pi:
        theta[i] = 2*np.pi - theta2[i] #[0, pi]

index = np.intersect1d(np.where(theta>1.5)[0], np.where(v<35)[0])
fig, ax = plt.subplots()
ax.scatter(data_np4[index+1, 1], data_np4[index+1, 2])
plt.show()
"""