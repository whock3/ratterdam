# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 19:34:26 2020

@author: Ruo-Yah Lai
"""
import numpy as np


def skaggs(pos, spikes):
    hs = np.histogram2d(spikes[:,2],spikes[:,1],bins=[np.arange(480), np.arange(640)])[0]
    hp = np.histogram2d(pos[:,2],pos[:,1],bins=[np.arange(480), np.arange(640)])[0]
    n = (hs*np.reciprocal(hp))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(hp==0)] = np.nan
    
    hp /= np.nansum(hp)
    meanRate = len(spikes)/len(pos)*30
    print(meanRate, np.nansum(n * np.log2(n/meanRate) * hp)/meanRate)
    return np.nansum(n * np.log2(n/meanRate) * hp)