# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 16:10:36 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
from alleyTransitions import alleyTransitions
from string import ascii_uppercase


def directionFilter(pos, spikes):
    """
    Returns position and spikes filtered by which direction the rat was facing

    """
    directions = np.diff(pos[:, 1:3], axis=0)
    allo = np.arctan2(directions[:, 1], directions[:, 0])
    posN = pos[:-1][np.logical_and((np.pi/4)<=allo, allo<(3/4*np.pi))]
    posE = pos[:-1][np.logical_and((-np.pi/4)<=allo, allo<(np.pi/4))]
    posS = pos[:-1][np.logical_and((-3/4*np.pi)<=allo, allo<(-1/4*np.pi))]
    posW = pos[:-1][np.logical_or((3/4*np.pi)<=allo, allo<(-3/4*np.pi))]
    
    [spikesN, spikesE, spikesS, spikesW] = [[] for _ in range(4)]
    for spike in spikes:
        if spike[0] in posN[:,0]:
            spikesN.append(spike)
        elif spike[0] in posE[:,0]:
            spikesE.append(spike)
        elif spike[0] in posS[:,0]:
            spikesS.append(spike)
        elif spike[0] in posW[:,0]:
            spikesW.append(spike)
    [spikesN, spikesE, spikesS, spikesW] = [np.array(i) for i in [spikesN, spikesE, spikesS, spikesW]]
    return [posN, posE, posS, posW], [spikesN, spikesE, spikesS, spikesW]


def occupancy(pos, title=""):
    """
    Makes a bar graph of the amount of time spent in each alley/intersection
    """
    posNew, turns = alleyTransitions(pos)
    
    letters = [i for i in ascii_uppercase[:12]]
    times = np.empty(0)
    for i in np.arange(17).astype(str):
        inLoc = np.where(posNew[:,4] == i)
        entries = np.where(posNew[inLoc,3] == str(1))[1]
        exits = np.where(posNew[inLoc,3] == str(2))[1]
        time = sum(posNew[inLoc][exits,0].astype(float) - posNew[inLoc][entries,0].astype(float))
        times = np.hstack((times,time))
    for i in letters:
        inLoc = np.where(posNew[:,4] == i)
        entries = np.where(posNew[inLoc,3] == str(3))[1]
        exits = np.where(posNew[inLoc,3] == str(4))[1]
        time = sum(posNew[inLoc][exits,0].astype(float) - posNew[inLoc][entries,0].astype(float))
        times = np.hstack((times,time))
    
    totalTime = sum(times)
    times /= totalTime
    fig, ax = plt.subplots(figsize=(8,4))
    ax.bar(list(np.arange(17).astype(str))+letters, times)
    ax.set_title(title)
    ax.set_ylabel("Fraction time in the location")
    ax.set_xlabel("Alleys and intersections")