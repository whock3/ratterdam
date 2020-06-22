# -*- coding: utf-8 -*-
"""

"""

from velocity_filtering import velocity_filtering
from RDP4_class import RDP4
from matplotlib.collections import LineCollection
from turn3 import turn3
import numpy as np
import matplotlib.pyplot as plt
import csv

def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)
"""
with open("pos.p.ascii", 'r') as csv_file:
    data_iter = csv.reader(csv_file)
    data = [data for data in data_iter]
data_np1 = np.array(data[50574:-17], dtype = "float64")
a = np.array([1, 4.72, 4.72])
data_np2 = data_np1[:,0:3]/a[None,:]
data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
data_np4 = velocity_filtering(data_np3)
data_np5 = velocity_filtering(data_np3, 3)


#sum of "turns" within a certain range
b = 43621 #48600 when velocity_threshold = 2
epsilons = [0.5, 0.75, 1]
i = 2
min_angle1 = np.pi/12 #min angle for the first angle
min_angle2 = np.pi*0.25 #min angle for angles within a range to add up to
turn_range = 5 #angles within x cm of each "turn" shouldn't add up to be a circle

RDP3.ResultList = data_np5[b]
RList = RDP3.RDPfunc(RDP3, data_np5[b:43903], epsilons[i]) 
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

fig, ax = plt.subplots(1, 1, figsize = (7, 5))
ax.set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax.set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax.scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax.scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")

idx3, theta_sum2, idx2_1 = turn(RList[:, 1:3], min_angle1, min_angle2, turn_range)

ax.scatter(RList[idx3,1], RList[idx3,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax.add_collection(lc)
axcb = fig.colorbar(lc)
axcb.set_label('Instantaneous velocity')
ax.set_title("velocity threshold = 3 \n epsilon = 1")

plt.legend() #bbox_to_anchor=(0.05, 0.55), loc="center left"
plt.show()
"""

#writing egocentric turn direction into csv
b1 = 43621
b2 = 43903
min_angle1 = np.pi/12 #min angle for the first angle
min_angle2 = np.pi*0.25 #min angle for angles within a range to add up to
turn_range = 5 #angles within x cm of each "turn" shouldn't add up to be a circle

def classifyTurns(b1, b2, position, min_angle1=np.pi/12, min_angle2=np.pi/4, turn_range_start=5*4.72, vCutoff=15*4.72, vIncrease=1/7, epsilon=4.72):
    RList = RDP4(position[b1:b2], epsilon).ResultList
    idx3, theta_sum2, idx2_1 = turn3(RList, min_angle1, min_angle2, turn_range_start, vCutoff, vIncrease)
    
    ego_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if theta_sum2[i] < (3/4*np.pi):
            ego_turns[i] = 0
        elif theta_sum2[i] > (5/4*np.pi):
            ego_turns[i] = 2
        else:
            ego_turns[i] = 1
    ego_turns = np.reshape(ego_turns, (len(idx3), 1))
            
    directions = np.diff(RList[:, 1:3], axis=0)
    allo = np.arctan2(directions[idx2_1, 1], directions[idx2_1, 0]) #in RList frame, length of idx3
    allo_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if (np.pi/4) <= allo[i] < (3/4*np.pi):
            allo_turns[i] = 0
        elif (3/4*np.pi) <= allo[i] or (-3/4*np.pi) > allo[i]:
            allo_turns[i] = 3
        elif (-3/4*np.pi) <= allo[i] < (-1/4*np.pi):
            allo_turns[i] = 2
        else:
            allo_turns[i] = 1
    allo_turns = np.reshape(allo_turns, (len(idx3), 1))
    
    times = np.reshape(RList[idx3, 0], (len(idx3), 1))
    return np.hstack((times, ego_turns, allo_turns))

"""
TimeAndTurns = ClassifyTurns(b1, b2, data_np5)
with open("egocentric turns", "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Timestamp", "Egocentric", "Allocentric"])
    csvwriter.writerows(TimeAndTurns)
"""