# -*- coding: utf-8 -*-
"""

"""
from velocity_filtering import velocity_filtering
from RDP3_class import RDP3
from matplotlib.collections import LineCollection
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

"""
#from https://matplotlib.org/3.2.1/gallery/lines_bars_and_markers/scatter_hist.html
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)
    plt.xlabel("instantaneous velocity")
    plt.ylabel("angle")

    # now determine nice limits by hand:
    binwidthx = 4
    binwidthy = np.pi/40
    xmax = np.max(np.abs(x))
    limx = (int(xmax/binwidthx) + 1) * binwidthx
    binsx = np.arange(0, limx + binwidthx, binwidthx)
    binsy = np.arange(0, np.pi, binwidthy)
    ax_histx.hist(x, bins=binsx)
    ax_histy.hist(y, bins=binsy, orientation='horizontal')

#2D histogram of angles and instantaneous velocities, different velocity thresholds
velocity = [2, 3, 4]
i = 0

data_np4 = velocity_filtering(data_np3, velocity[i])
diffts, diffx, diffy = np.diff(data_np4[:,0]), np.diff(data_np4[:,1]), np.diff(data_np4[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))
directions = np.diff(data_np4[:,1:3], axis=0)
theta = angle(directions)  

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)
v = v[1:len(v)]
idx = np.where(v < 500)[0]

scatter_hist(v[idx], theta[idx], ax, ax_histx, ax_histy)
plt.title(f"velocity threshold = {velocity[i]}")
plt.show()
"""

b = 48600
epsilons = [0.5, 0.75, 1]
i = 2
min_angle = np.pi*0.25

#with DRP and velocity filter
RDP3.ResultList = data_np4[b]
RList = RDP3.RDPfunc(RDP3, data_np4[b:b+300], epsilons[i]) 
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

fig, ax = plt.subplots(1, 3, figsize = (19, 5), gridspec_kw={'width_ratios': [6, 6, 7]})
ax[0].set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax[0].set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax[0].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax[0].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")
directions = np.diff(RList[:, 1:3], axis=0)
theta = angle(directions)
idx = np.where(theta>min_angle)[0]+1
ax[0].scatter(RList[idx,1], RList[idx,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[0].add_collection(lc)
axcb = fig.colorbar(lc)
axcb.set_label('Instantaneous velocity')
ax[0].set_title("velocity threshold = 2 \n epsilon = 1")

#sanity check with no RDP but has velocity filter
diffts, diffx, diffy = np.diff(data_np4[b:b+300,0]), np.diff(data_np4[b:b+300,1]), np.diff(data_np4[b:b+300,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

ax[1].set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax[1].set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax[1].scatter(data_np4[b,1], data_np4[b,2], marker = "+", color = "r", label = "first")
ax[1].scatter(data_np4[b+299,1], data_np4[b+299,2], marker = "x", color = "r", label = "last")

RList2 = np.reshape(data_np4[b:b+300, 1:3], (300, 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[1].add_collection(lc)
ax[1].set_title("velocity threshold = 2 \n No RDP")

#sanity check with no velocity filter
b1 = np.where(data_np3 == data_np4[b][0])[0]
b2 = np.where(data_np3 == data_np4[b+300][0])[0]
diffts, diffx, diffy = np.diff(data_np3[b1[0]:b2[0],0]), np.diff(data_np3[b1[0]:b2[0],1]), np.diff(data_np3[b1[0]:b2[0],2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

ax[2].set_xlim(np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
ax[2].set_ylim(np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
ax[2].scatter(data_np3[b1,1], data_np3[b1,2], marker = "+", color = "r", label = "first")
ax[2].scatter(data_np3[b2,1], data_np3[b2,2], marker = "x", color = "r", label = "last")

RList2 = np.reshape(data_np3[b1[0]:b2[0], 1:3], (b2[0]-b1[0], 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[2].add_collection(lc)
ax[2].set_title("No velocity filter")

plt.legend()
plt.show()
