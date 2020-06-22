# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from bisect import bisect_left
import scipy.ndimage
from pathlib import Path



def read_clust(clustfile):
    '''open a cl-mazeX.X file and read spike times into a list'''
    data_folder = Path("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/")
    file_to_open = data_folder / clustfile
    with open(file_to_open) as clust:
        raw_clust_data = clust.read().splitlines()
    spikes = []
    for line in raw_clust_data[13:]:
        spikes.append(float(line.split(',')[-1]))
    return spikes
 
def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def getPosFromTs(yourts,position):
    '''position is the full list of ts, yourts are the times you'd like a pos for '''
    adjTs = [takeClosest(position[:,0],i) for i in yourts]
    target = np.empty((0,3))
    for i in adjTs:
        target = np.vstack((target,position[np.where(position[:,0]==i)][:,0:][0]))
    return target

def weird_smooth(U,sigma):
    V=U.copy()
    V[U!=U]=0
    VV=scipy.ndimage.gaussian_filter(V,sigma=sigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=scipy.ndimage.gaussian_filter(W,sigma=sigma)

    Z=VV/WW
    return Z

def makeRM(spikes, position, bins=[30,50], smoothing_2d_sigma=2):
    """ this is 2d whole track. uses g"""
    rbins, cbins = wholeAlleyBins
    rows, cols = np.linspace(0, 480, rbins), np.linspace(0,640, cbins)
    hs,xe,ye = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rows, cols])
    ho = np.histogram2d(position[:,2],position[:,1],bins=[rows, cols])[0]
    n = (hs*np.reciprocal(ho))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(ho==0)] = np.nan
    n = weird_smooth(n,smoothing_2d_sigma)
    n[np.where(ho==0)] = np.nan
    return n

ptsCm = 4.72 # This is for Macaulay. The kriger number is 4.85
wholeAlleyBins = [round(int(480/ptsCm)),int(round(640/ptsCm))]

"""
#Run from here
spikes = read_clust("cl-maze1.1")
target = getPosFromTs(spikes, data_np5) #data_np5 without adjusting from camaera frame to cm
n = makeRM(target, data_np5)
fig, ax = plt.subplots(1, 1)
ax.set_title("R859  Ratterdam Open Day 3 \n Tetrode 6  cl-maze1.1")
plt.imshow(n)
"""