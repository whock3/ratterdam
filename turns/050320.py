# -*- coding: utf-8 -*-
"""

"""
#after running top part of 042320
"""
#comparing velocity filters
b = 55600
epsilon = 1
min_angle = np.pi*0.25
fig, ax = plt.subplots(1, 3, figsize=(24, 6))

RDP2.ResultList = data_np4[b, 1:3]
RList = RDP2.RDPfunc(RDP2, data_np4[b:b+300, 1:3], epsilon)
directions = np.diff(RList, axis=0)
theta = angle(directions)
idx = np.where(theta>min_angle)[0]+1
percent = round(100 - len(RList)/3, 2)
ax[0].plot(RList[:,0], RList[:,1])
ax[0].set_title(f"epsilon = {round(epsilon, 2)}" "\n velocity threshold = 2 \n 300 points")
ax[0].scatter(RList[0,0], RList[0,1], marker = "+", color = "r", label = "first")
ax[0].scatter(RList[-1,0], RList[-1,1], marker = "x", color = "r", label = "last")
ax[0].scatter(RList[idx,0], RList[idx,1], s = 16, color = "r", label='turning points')

RDP2.ResultList = data_np5[50073, 1:3]
RList = RDP2.RDPfunc(RDP2, data_np5[50073:50365, 1:3], epsilon)
directions = np.diff(RList, axis=0)
theta = angle(directions)
idx = np.where(theta>min_angle)[0]+1
percent = round(100 - len(RList)/3, 2)
ax[1].plot(RList[:,0], RList[:,1])
ax[1].set_title(f"epsilon = {round(epsilon, 2)}" "\n velocity threshold = 3 \n 292 points")
ax[1].scatter(RList[0,0], RList[0,1], marker = "+", color = "r", label = "first")
ax[1].scatter(RList[-1,0], RList[-1,1], marker = "x", color = "r", label = "last")
ax[1].scatter(RList[idx,0], RList[idx,1], s = 16, color = "r", label='turning points')

RDP2.ResultList = data_np6[46447, 1:3]
RList = RDP2.RDPfunc(RDP2, data_np6[46447:46716, 1:3], epsilon)
directions = np.diff(RList, axis=0)
theta = angle(directions)
idx = np.where(theta>min_angle)[0]+1
percent = round(100 - len(RList)/3, 2)
ax[2].plot(RList[:,0], RList[:,1])
ax[2].set_title(f"epsilon = {round(epsilon, 2)}" "\n velocity threshold = 4 \n 269 points")
ax[2].scatter(RList[0,0], RList[0,1], marker = "+", color = "r", label = "first")
ax[2].scatter(RList[-1,0], RList[-1,1], marker = "x", color = "r", label = "last")
ax[2].scatter(RList[idx,0], RList[idx,1], s = 16, color = "r", label='turning points')

plt.legend()
plt.show()
"""




from velocity_filtering import velocity_filtering
from RDP3_class import RDP3
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
    binwidthx = 1
    binwidthy = np.pi/40
    xmax = np.max(np.abs(x))
    limx = (int(xmax/binwidthx) + 1) * binwidthx
    binsx = np.arange(0, limx + binwidthx, binwidthx)
    binsy = np.arange(0, np.pi, binwidthy)
    ax_histx.hist(x, bins=binsx)
    ax_histy.hist(y, bins=binsy, orientation='horizontal')



#2D histogram of angles and instantaneous velocities, different epsilons
b = 48600
epsilons = [0.5, 0.75, 1]
i = 2

RDP3.ResultList = data_np4[b]
RList = RDP3.RDPfunc(RDP3, data_np4[b:b+300], epsilons[i]) 
diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))
directions = np.diff(RList[:,1:3], axis=0)
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

scatter_hist(v[1:len(v)], theta, ax, ax_histx, ax_histy)
plt.title(f"epsilon = {epsilons[i]}")
plt.show()


