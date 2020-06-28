# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 18:12:25 2020

@author: Ruo-Yah Lai
"""
import numpy as np
from turn4 import turn4
from turn3 import turn3
import matplotlib.pyplot as plt

#from 061220: velSegments, turn3 vs turn4
#RList in camera units after RDP
RList = []
a = np.array([1,4.72,4.72])
RList = RList/a[None,:]

diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
v =  np.sqrt((vx**2)+(vy**2))

fig, ax = plt.subplots(1, 2, figsize = (12, 5), gridspec_kw={'width_ratios': [5, 6]})
ax[0].set_xlim(5,120)
ax[0].set_ylim(5, 120)
ax[0].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax[0].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")
ax[0].set_xlabel("x coordinates (cm)")
ax[0].set_ylabel("y coordinates (cm)")

idx3,_,_ = turn3(RList, np.pi/12, np.pi/4, 5, 15)

ax[0].scatter(RList[idx3,1], RList[idx3,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:,1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[0].add_collection(lc)
ax[0].set_title("turn3 (turn range depends on instantaneous velocity)\n" + f"{len(idx3)} turns")

#ax[0].axis("equal")

#turn4
ax[1].set_xlim(5,120)
ax[1].set_ylim(5,120)
ax[1].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
ax[1].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")
ax[1].set_xlabel("x coordinates (cm)")
ax[1].set_ylabel("y coordinates (cm)")

idx4,_,_,trs, percents = turn4(RList)

ax[1].scatter(RList[idx4,1], RList[idx4,2], s = 16, color = "r", label='turning points')

RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
lc = LineCollection(segments)
lc.set_array(v)
ax[1].add_collection(lc)
axcb = fig.colorbar(lc)
axcb.set_label('Instantaneous velocity (cm/s)')
ax[1].set_title("turn4 \n" + f"{len(idx4)} turns")

plt.legend() #bbox_to_anchor=(0.05, 0.55), loc="center left"
#ax[1].axis("equal")
plt.show()



#graphing trajectories after shifting and rotating them
fig,axs = plt.subplots(2,4)
titles = ["All turns", "Left", "Back", "Right", "North", "East", "South", "West"]
fig.suptitle("R859 D3 T6 1.1 subfield 0", y=1.04)
srP = 
    
for i in range(2):
    for j in range(4):
        axs[i][j].set_title(titles[i*4+j] + "\n" + f"n = {len(srP[i*4+j])}")
        for k in range(len(srP[i*4+j])):
            axs[i][j].plot(srP[i*4+j][k][:,1], srP[i*4+j][k][:,2])
            axs[i][j].axis("equal")
fig.tight_layout()

