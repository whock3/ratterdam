# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 15:38:58 2022

@author: whockei1

Identify reward deliveries without using the TTL signal (that normally signals them)
This is in case where we think TTL wasn't working right, like R781

Look for 1) low speeds 2) near lickports 3) with a licking artifact 

"""
import numpy as np
import utility_fx as util
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import pandas as pd
import ratterdam_DataFiltering as Filt 
import utility_fx as util 
import newAlleyBounds as nab
from scipy.signal import butter, lfilter 


proxThresh = 75 # proximity threshold in camera frame coord units within which
                # the animal might be receiving a reward 
 

datafile = 'E:\\Ratterdam\\R859\\R859_RatterdamOpen_D2\\'
    
# data format for Neuralynx CSC data 
ncsType = np.dtype([
    ('ts', '<u8'),
    ('dwChannelNumber', '<u4'),
    ('fs', '<u4'),
   ('NumValidSamples', '<u4'),
    ('data', '<i2', (512, ))
])
                
def createCSCData(tt):
    """
    Given an electrode as input and the path to a recording data directory
    in the namespace:
    Load csc data as an (n,2) array of (ts, sample)
    Reminder that neuralynx stores data in chunks with one ts so we flatten those
    and interpolate the ts to get the final array returned 
    """

    fname = datafile + f"CSC{tt}.ncs"
    
    curFile = open(fname)
    curFile.seek(16*1024)
    curFile = np.fromfile(curFile, dtype=ncsType)
    
    dt = (1/32000)*1e6 # this is the timestep between each sample in packet. Only first sample is explicitly timestamped\n",
                         # by neuralynx at time of collection so increment this dt value repeatedly to get ts of all other samples in packet\n",
    cscdata = curFile['data'].flatten()
    cscts = []
    for packet in curFile:
        init_ts = packet['ts']
        cscts.append(init_ts)
        interp_ts = init_ts + dt # initialize the, basically, counter. this init is for the second point if that makes sense\n",
        for datapoint in packet['data'][1:]: # start with second bc first one is (the only one) timestamped
            cscts.append(interp_ts)
            interp_ts += dt # increment for next ts
            
    csc = np.column_stack((cscts, cscdata))
    
    return csc


def filter_data(data, low,high, sf=32000, order=2):
    """
    Bandpass filter [low, high] using a butterworth filter
    -If array input is 1d it should be the samples
    -If array input is 2d it should be (n,2) it should be (ts, samples)
     and return will be (n,2) like input
    """
    # Determine Nyquist frequency\n",
    nyq = sf/2
    # Set bands
    low = low/nyq
    high = high/nyq
    # Calculate coefficients
    b, a = butter(order, [low, high], btype='band')
    # Filter signal
    if data.ndim == 1:
        filter_input = data
    elif data.ndim == 2:
        filter_input = data[:,1]
    filtered_data = lfilter(b, a, filter_input)
    
    # if user passed (n,2) array with ts as first column, stack it back on before returning 
    if data.ndim == 2:
        filtered_data = np.column_stack((data[:,0], filtered_data))
    
    return filtered_data


def extractCSCSegment(csc, tsBeg, tsEnd):
    """
    Input - * csc is (n,2) array with whole session csc data\n",
            for one tt (sampled 1/32000 s). ts, sample value\n",
           * tsBeg, tsEnd - neuralynx ts for segment\n",
           start and end\n",
    Return: 
    """
    mask = (csc[:,0]>tsBeg)&(csc[:,0]<tsEnd)
    segment = csc[mask,:]
    return segment



rat, day = 'R859', 'D2'
ratborders = nab.loadAlleyBounds(rat, day)
aib = ratborders.alleyInterBounds
turns, unit = RepCore.loadTurns(rat, day)
smoothspeed = Filt.computeSpeed(unit.position)
rewards = RepCore.readinRewards(rat, day)
centers = []
for i in range(17):
    b = aib[str(i)]
    x,y = (b[0][0]+b[0][1])/2, (b[1][0]+b[1][1])/2
    centers.append([x,y])
    

cdists = []
for c in centers:
    cdists.append(np.linalg.norm(c-unit.position[:,1:],axis=1))


#highpass_csc = filter_data(csc, 20,100)
highpass_csc = csc
fig, ax = plt.subplots()
ax.vlines(rewards,0,50,color='g')
ax.plot(unit.position[:,0],smoothspeed)
for c in cdists:
    prox  = np.where(c < proxThresh)[0]
    proxts = unit.position[prox,0]
    proxspeed = smoothspeed[prox]
    ax.scatter(proxts,proxspeed,s=20)
    
ax2 = ax.twinx()
ax2.plot(highpass_csc[:,0], highpass_csc[:,1],color='r')