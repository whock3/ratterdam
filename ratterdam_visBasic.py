# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 15:18:39 2018

@author: whockei1
"""

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket
from scipy.stats import sem
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
from matplotlib.colors import LinearSegmentedColormap
import sys

if socket.gethostname() == 'Tolman':
    codeDirBase = 'C:\\Users\\whockei1\\Google Drive'
elif socket.gethostname() == 'DESKTOP-BECTOJ9':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'
elif socket.gethostname() == 'DESKTOP-0PH2B0B':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'
    
sys.path.insert(0, codeDirBase + '\\KnierimLab\\Ratterdam\\Code')
sys.path.insert(0, codeDirBase + '\\Python_Code\\KLab\\mts_analysis')
import utility_fx as util
import ratterdam_ParseBehavior as pBehav
from ratterdam_Defaults import *


class BasicRateMaps():
    """
    Class with plotting fx for:
    - 2d ratemaps
    - 1d ratemaps (computed differently than 2d)
    - spike maps (spikes over trajectory)
    """
    def __init__(self):
        self.subplotColNum = 4
        self.currArrayMax = 0 # current maximum of arrays being plotted so that
                               # acoss-plot entites like colorbars can have a common scale
        self.rmDims = {'long':30, 'short':15} #todo: dont hardcode this. lookup on init from ratterdam_Defaults.
        self.gs_lookup = {1:[0,1], 2:[2,1], 3:[1,0], 4:[1,2], 5:[0,3],
                          6:[1,4], 7:[0,5], 8:[1,6], 9:[2,5], 10:[3,6],
                          11:[4,5], 12:[3,4], 13:[2,3], 14:[4,3], 15:[3,2],
                          16:[4,1], 17:[3,0]}
        self.colormap = self.makeCustomColormap()
        
    
    def makeCustomColormap(self, nb=100,name='mymap',c=[]):

        if c ==[]:
            c = [(0,0,0.5),(1,1,0),(1,0,0)]
        mycm = LinearSegmentedColormap.from_list(name,c,N=nb)
        return mycm
    
    def generateSingleSpikeMap(self, spikes, pos, alley, **kwargs):
        """
        Spike/Pos data n,3 array ts,x,y
        """
        self.ax.plot(pos[:,1], pos[:,2],'k',alpha=0.5)
        self.ax.scatter(spikes[:,1], spikes[:,2], c='r')
        self.drawAlleyBounds(self.ax, alley)
        self.annotate(**kwargs)
        
    def drawAlleyBounds(self, ax, alley, **kwargs):
        """
        Assume 1-idx entry
        """
        alley = alley - 1
        a = alleyBounds[alley]
        x1,y1,x2,y2 = a[0][0], a[1][0], a[0][1], a[1][1]
        for x,y in zip([[x1, x1], [x1, x2], [x2, x2], [x1, x2]], [[y1, y2], [y2, y2], [y1, y2], [y1, y1]]):
            ax.plot(x, y, 'k')
            
    def compute_AvgTrace(self, unit, alley, txt):
        """
        Creates a rate map for each visit and returns
        avg trace and sem trace
        
        Because defining linear RM by avg of visits
        is only done for visualization purposes, the code
        to compute it (this fx) is in the plotting class 
        rather than the unit data class
        """
        lins = np.empty((0, self.rmDims['long']))
        for visit in unit.alleys[alley]:
            if visit['metadata']['stimulus'] == txt:
                spk, occ = visit['spikes'], visit['occs']
                if occ.shape[0] > 0: #yes all visits should, by def, have occupancy but you never know
                    lin = unit.computeSingleRM(spk, occ, alley, dim=1)
                    lins = np.vstack((lins, lin))
                    
        avg = np.nanmean(lins, axis=0) 
        err = sem(lins, axis=0, nan_policy='omit')
        return avg, err        
            
    def annotate(self, **kwargs):
        if 'title' in kwargs:
            self.ax.set_title(kwargs['title'])
            
    def plot_SpikemapAllVisits(self, unit, alley):
        """
        Generates subplots each of which a spike map of a visit to alley
        
        Reminder on unit structure. unit.alleys[alley] is a list of dicts
        Ea. has spikes, occ, rm, metadata
        """
        alleydata = unit.alleys[alley]
        self.fig, self.axes = plt.subplots(int(np.ceil(len(alleydata)/self.subplotColNum)), self.subplotColNum)
        for i in range(len(alleydata)):
            self.ax = self.fig.axes[i]
            spikes = alleydata[i]['spikes']
            occs   = alleydata[i]['occs']
            txt    = alleydata[i]['metadata']['stimulus']
            print(spikes.shape)
            print(occs.shape)
            print("-----------")
            self.generateSingleSpikeMap(spikes, occs, alley, **{'title':txt})
            

    def getMaxArrays(self, unit, stim, cutoff = 0.98):
        """
        Compute maximum value of colormap for 2d ratemap
        This max defined as the firing rate which, at
        or below this value, accounts for cutoff proportion
        of the total firing rate occurances (discretized in 100 bins)
        Essentially the histogram bin that accounts for cutoff proportion
        of the cumulative sum of the bin counts.
        """
        frs = []
        
        #This conditional is here because I want
        # to scale alleys in 'overall' plots to one colorbar
        # scale and the alleys in plots for A, B, C to a separate
        # common-to-them colorbar. So grab different sets of arrays to get max
        if stim  == 'textures':
            arrayTypes = ['A', 'B', 'C']
        elif stim == 'overall':
            arrayTypes = ['overall']
            
        for arrayType in arrayTypes:
            for i in range(1,18):
                rm = unit.alleyRMS[i][arrayType]
                if type(rm) == np.ndarray:
                    rm = np.ndarray.flatten(rm)
                    if not all(np.isnan(rm)):
                        frs.extend(rm)
        frs = np.asarray(frs)
        frs = frs[np.isfinite(frs)]
        h,b = np.histogram(frs, bins=100)
        frcum = np.cumsum(h)
        propExp = np.asarray([i/h.sum() for i in frcum])
        try:
            thresh = np.where(propExp < cutoff)[0][-1]
        except:
            thresh = np.where(b == np.median(b))
        return b[thresh]*2.5
    
    
    def getMaxLins(self, unit):
        """
        Get maximum y limit for plotting linear ratemaps
        Logic is substantially different from the 2d version
        so it's its own fx.
        Scans all rms across all alleys, finds max and sets the 
        max to be like .2* higher than that. 
        """
        mymax = 0
        for alley in range(1,18):
            for txt in ['A', 'B', 'C']:
                avg, err = self.compute_AvgTrace(unit, alley, txt)
                if np.nanmax(avg+err) > mymax:
                    mymax = np.nanmax(avg+err)
        mymax = mymax*1.2
        return mymax
    
    def calcGoodMax(self, array):
        frs = np.ndarray.flatten(array)
        frs = frs[np.isfinite(frs)]
        h,b = np.histogram(frs, bins=100)
        frcum = np.cumsum(h)
        try:
            thresh = np.where(propExp < 0.97)[0][-1]
        except:
            thresh = np.where(b == np.median(b))
        return b[thresh]*1.5
    
    
    def add_colorbar(self, fig, im):
        #left bottom width heigh
        cax = fig.add_axes([.2,.95,.6,.05])
        fig.colorbar(im, cax=cax, orientation='horizontal', cmap=self.colormap)
        
    def plot_2dAlleys(self, ax, unit, alley, data,**kwargs):
        '''Plot 2d (imshow) of alley
         Data == "overall" gives all visits
         Data == A,B,C gives you those visits, collapsed.
         If you say data == "texture" and dont give a kwarg['texture'] value ERROR
        '''
        im = False
        rm = unit.alleyRMS[alley][data]
        if type(rm) == np.ndarray:
            im = ax.imshow(rm, origin='lower', interpolation = 'None', vmax= kwargs['max'], cmap = self.colormap)
        
        # this ugly workaround is to get an imshow result
        # from at least one alley that i can use to
        # generate colorbar(). Wish i knew how to generate
        # a colobar without a plot 
        if kwargs and 'return' in kwargs:
            if kwargs['return'] == True:
                return im
            
    def plot_linAlley(self, ax, unit, alley, txt="All", **kwargs):
        """Will either just plot 1 txt if given
        or if txt='All' then all 3 traces. Traces
        are NOT linear RMS, but are avg/sem of individual trial
        RMs. For quantitative analyses eg perm tests
        lin RMS as defined by unit.linRMS are used. These
        are essentially collapsing acros trials to get trace. Here,
        and similar viz routines, use avging of trial rms. Methods yield different
        results. See Unit doc for details"""
        if txt.lower() == 'all':
            for txt,color in zip(['A', 'B', 'C'], ['r', 'b', 'g']):
                avg, err = self.compute_AvgTrace(unit, alley, txt)
                ax.plot(avg, f"{color}")
                ax.fill_between(range(len(avg)), avg+err, avg-err, color=f"{color}", alpha=0.5)
            ax.set_title(f"{alley}")
            if kwargs and 'ylim' in kwargs:
                ax.set_ylim([0, kwargs['ylim']])
        else:
            return "Further functionality not yet implemented. Please stick to defaults"
        
    def plot_2dWholeTrack(self, position, unit):
        """
        Will plot the whole track firing rate map.
        Does not take an axis arg - makes it's own sep figure.
        Does not break data into alleys. All spikes / all pos
        ratterdam_Defaults as bin sizes.
        Camera resolution of 640 (cols) x 480 (rows) is hardcoded.
        """
        fig, ax = plt.subplots(figsize=(6,6))
        col, row = np.linspace(0, 640, wholeAlleyBins[0]), np.linspace(0,480, wholeAlleyBins[1])
        hs = np.histogram2d(unit.spikes[:,2],unit.spikes[:,1],bins=[row, col])[0]
        ho = np.histogram2d(position[:,2],position[:,1],bins=[row, col])[0]
        n = (hs*np.reciprocal(ho))*33
        n = util.weird_smooth(n,1)
        n[np.where(ho==0)] = np.nan
        mymax = util.calcSmartMax(n)
        im = ax.imshow(n, aspect='auto', interpolation='None', origin='lower', vmax=mymax, cmap=self.colormap)
        self.add_colorbar(fig, im)
        ax.set_axis_off()
        #plt.tight_layout()
        ax.set_title(f"Overall ratemap for {unit.name}")
