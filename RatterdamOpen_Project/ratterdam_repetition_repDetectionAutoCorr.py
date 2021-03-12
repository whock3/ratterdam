# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 18:11:43 2021

@author: whockei1

Place Field Repetition Project
Place Field Repetition Detection

This script will identify single units (done per cell or in bulk) that display
repeating fields according to a 2D autocorrelation grid peak detection.
"""

import numpy as np
import scipy.signal

from scipy.ndimage.filters import maximum_filter
#from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import skimage
import scipy.fftpack

import utility_fx as util
import os
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore

# %% 
rat = "859"
day = "D2"
df = 'E:\\Ratterdam\\R859\\R859_RatterdamOpen_D2\\'
clustname = "TT14\\cl-maze1.3"
unit = RepCore.loadRepeatingUnit(df, clustname, smoothing=1)


#%%
RepCore.plotRoutine_RepPF_TempDyn(unit, nf=99, time='time', save=False)

# %%
rm = unit.repUnit.rateMap2D

rm = np.nan_to_num(rm)
crm = scipy.signal.correlate2d(rm,rm)

#%% Find Peaks and Plot onto auto-correlation plot 
peaks = skimage.feature.peak_local_max(crm,min_distance=5)
plt.figure()
plt.imshow(crm,origin='lower',vmax=np.nanpercentile(crm,99))
plt.scatter(peaks[:,1], peaks[:,0],color='r',s=20) # cols then rows

#%%
# get manhattan dist between each pair of dists
manhattans = []
for a in peaks:
    for b in peaks:
        manhattans.append(scipy.spatial.distance.cityblock(a,b))
        
sm_manhattans = scipy.ndimage.gaussian_filter(manhattans,sigma=0)
hist, edges=np.histogram(sm_manhattans,bins=80)
plt.figure();plt.plot(hist)

# fft of smoothed histogram
x = np.linspace(0, 1/2, int(len(hist)/2))
y = scipy.fftpack.fft(hist)
plt.figure();plt.plot(x[:], 2/len(hist)*np.abs(y[:int(hist.shape[0]//2)]))
plt.title("FFT of (Smoothed) Histogram of Manhattan Distances")


 






