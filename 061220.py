# -*- coding: utf-8 -*-
"""

"""
import csv
import numpy as np
from bisect import bisect
from velocity_filtering import velocity_filtering
from RateMap import read_clust, getPosFromTs
from UnitClass import Unit
from organized import classifyTurns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

def read_pos(pos="pos.p.ascii"):
    data_folder = Path("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Raw data/")
    file_to_open = data_folder / pos
    with open(file_to_open, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    data_np1 = np.array(data[50574:-17], dtype = "float64")
    #a = np.array([1, 4.72, 4.72])
    data_np2 = data_np1[:,0:3] #/a[None,:]
    data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
    #data_np4 = velocity_filtering(data_np3)
    return velocity_filtering(data_np3, 3*4.72)

s = getPosFromTs(read_clust("cl-maze1.1"), data_np5)
unit = Unit(s,data_np5)
#repUnit=RateMapClass.RateMap(unit) # unit is an class instance that has attributes unit.spikes, unit.position
"""

def turnsInOutFields(visits, position, subfield):
    ranges = [] #list of lists with 3s before and after visits to subfields 
                #(not 3s of activity)
    for i in range(len(visits[subfield][0])): #the ith visit
        #start and end are indexes of start and end of +- 3s of visit
        start = bisect(position[:,0], (visits[subfield][0][i][0] - 3e6))
        end = bisect(position[:,0], (visits[subfield][0][i][-1] + 3e6))
        ranges.append([start, end])
    
    turns = np.empty((0, 4)) #columns: ts, ego, allo, before/during/after visit
    for i in range(len(ranges)): #the ith visit
        turns1 = classifyTurns(ranges[i][0], ranges[i][1], position)
        bda = np.empty((0,1))
        for j in range(len(turns1)): #the jth turn
            if turns1[j][0] < visits[subfield][0][i][0]:
                bda = np.vstack((bda,0))
            elif turns1[j][0] > visits[subfield][0][i][-1]:
                bda = np.vstack((bda,2))
            else:
                bda = np.vstack((bda,1))
        turns1 = np.hstack((turns1, bda))
        turns = np.vstack((turns, turns1))
    return turns


#turns = turnsInOutFields(unit.visits, data_np5, 0)
"""
#3D histogram
#from https://matplotlib.org/3.1.0/gallery/mplot3d/hist3d.html#sphx-glr-gallery-mplot3d-hist3d-py
fig = plt.figure()
ax1 = fig.add_subplot(121, projection="3d")
ax2 = fig.add_subplot(122, projection="3d")
#ego turns
hist, xe, ye = np.histogram2d(turns[:,1], turns[:,3], bins=3, range=[[-0.5,2.5],[-0.5,2.5]])
xpos, ypos = np.meshgrid(xe[:-1] + 0.25, ye[:-1] + 0.25, indexing="ij")
xpos = xpos.ravel()
ypos = ypos.ravel()
zpos = 0
dx = dy = 0.5 * np.ones_like(zpos)
dz = hist.ravel()
ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
#allo turns
hist, xe, ye = np.histogram2d(turns[:,2], turns[:,3], bins=[4,3], range=[[-0.5,3.5],[-0.5,2.5]])
xpos, ypos = np.meshgrid(xe[:-1] + 0.25, ye[:-1] + 0.25, indexing="ij")
xpos = xpos.ravel()
ypos = ypos.ravel()
zpos = 0
dx = dy = 0.5 * np.ones_like(zpos)
dz = hist.ravel()
ax2.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
plt.show




#2D histogram with heatmap
#from  https://matplotlib.org/3.1.0/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
fig, ax = plt.subplots(1, 2)
bda_labels = ["before", "during", "after"]

#ego turns
hist, xe, ye = np.histogram2d(turns[:,1], turns[:,3], bins=3, range=[[-0.5,2.5],[-0.5,2.5]])
egoturns = ["left", "back", "right"]
im = ax[0].imshow(hist)
#We want to show all ticks...
ax[0].set_xticks(np.arange(len(bda_labels)))
ax[0].set_yticks(np.arange(len(egoturns)))
# ... and label them with the respective list entries
ax[0].set_xticklabels(bda_labels)
ax[0].set_yticklabels(egoturns)
# Rotate the tick labels and set their alignment.
plt.setp(ax[0].get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
# Loop over data dimensions and create text annotations.
for i in range(len(egoturns)):
    for j in range(len(bda_labels)):
        text = ax[0].text(j, i, hist[i, j],
                       ha="center", va="center", color="w")
ax[0].set_title("Egocentric turns")
ax[0].set_xlabel("Relative to visits")
ax[0].set_ylabel("Turn direction")

#allo turns
hist, xe, ye = np.histogram2d(turns[:,2], turns[:,3], bins=[4,3], range=[[-0.5,3.5],[-0.5,2.5]])
alloturns = ["north", "east", "south", "west"]
im = ax[1].imshow(hist)
#We want to show all ticks...
ax[1].set_xticks(np.arange(len(bda_labels)))
ax[1].set_yticks(np.arange(len(alloturns)))
# ... and label them with the respective list entries
ax[1].set_xticklabels(bda_labels)
ax[1].set_yticklabels(alloturns)
# Rotate the tick labels and set their alignment.
plt.setp(ax[1].get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
# Loop over data dimensions and create text annotations.
for i in range(len(alloturns)):
    for j in range(len(bda_labels)):
        text = ax[1].text(j, i, hist[i, j],
                       ha="center", va="center", color="w")
ax[1].set_title("Allocentric turns")
ax[1].set_xlabel("Relative to visits")
ax[1].set_ylabel("Turn direction")

fig.tight_layout()
plt.show()
"""


"""
#same graphs as before with vel segments, turn2 vs turn3
#data_np5 in camera units
b1 = 43621 #48600 when velocity_threshold = 2
b2 = 43903
min_angle1 = np.pi/12 #min angle for the first angle
min_angle2 = np.pi*0.25 #min angle for angles within a range to add up to
turn_range = 5 #angles within x cm of each "turn" shouldn't add up to be a circle

RDP3.ResultList = data_np5[b1]
RList = RDP3.RDPfunc(RDP3, data_np5[b1:b2], 4.72) 
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))/4.72

fig, ax = plt.subplots(1, 2, figsize = (12, 5), gridspec_kw={'width_ratios': [5, 6]})
ax[0].set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax[0].set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax[0].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax[0].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")

idx3, theta_sum2, idx2_1 = turn2(RList[:,1:], min_angle1, min_angle2, 5*4.72)

ax[0].scatter(RList[idx3,1], RList[idx3,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[0].add_collection(lc)
ax[0].set_title("turn2 (fixed turn range) \n velocity threshold = 3 cm/s \n epsilon = 1 cm")

ax[0].axis("equal")

#turn3
ax[1].set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax[1].set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax[1].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax[1].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")

idx3, theta_sum2, idx2_1 = turn3(RList, min_angle1, min_angle2, 5*4.72, 15*4.72, 1/7)

ax[1].scatter(RList[idx3,1], RList[idx3,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[1].add_collection(lc)
axcb = fig.colorbar(lc)
axcb.set_label('Instantaneous velocity (cm/s)')
ax[1].set_title("turn3 (variable turn range) \n velocity threshold = 3 cm/s \n epsilon = 1 cm")

plt.legend() #bbox_to_anchor=(0.05, 0.55), loc="center left"
ax[1].axis("equal")
plt.show()
"""