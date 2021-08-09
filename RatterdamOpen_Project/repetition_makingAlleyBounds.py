# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 18:45:39 2021

@author: whockei1

Making alley bounds for R886 and R765

RYL did so with her own code (likely not saved just done by hand) for rest
"""


# Note for R886 when I initially went to look at the pos data it seemed flipped
# on the xaxis. Meaning we usually place rat in maze in upper left corner 
# wrt to the standard 'viewer position' which is below alley 14, looking up
# it looked as if i placed the rat in the upper right corner which i likely didnt
# the camera orientation had an x flip applied and when i removed it things
# looked normal. WH 8-6-21
#%%

import matplotlib.pyplot as plt
import numpy as np

#%% R886


v = np.asarray([87, 132, 218, 262, 352, 402, 490, 533])
h = np.asarray([83, 137, 210, 270, 350, 402])


plt.figure()
plt.scatter(unit.position[:,1], unit.position[:,2])
plt.title("R886 No cam flip")

for i in v:
    plt.vlines(i,0,480,color='k')
    
for i in h:
    plt.hlines(i,0,640, color='k')




#%% R765

v = np.array([97, 150, 240, 287, 383, 425, 515, 557])
h = np.array([75, 123, 211, 270, 350, 410])


plt.figure()
plt.scatter(pos[:,1], pos[:,2])
plt.title("R765")

for i in v:
    plt.vlines(i,0,480,color='k')
    
for i in h:
    plt.hlines(i,0,640, color='k')