# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 18:27:43 2020

@author: whockei1
"""


import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
from scipy.stats import sem
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
from matplotlib.backends.backend_pdf import PdfPages
import more_itertools, itertools
from sklearn.metrics import auc
from scipy.interpolate import splrep, splev
from scipy.spatial import ConvexHull
import scipy
import ratterdam_DataFiltering as Filt


def loadRepeatingUnit(df, clustName, smoothing=2, vthresh=3):
    """take a path to a data dir
    load spikes and position into two np arrays
    spikes is (n,1) and pos is typical (3,n) cols of ts,x,y
    use cameraOrientationInfo.txt to flip axes if needed
    use sessionEpochInfo.txt, specific for open Ratterdam exp
    to get session ts and clip spikes/pos"""
    
    with open(df+"sessionEpochInfo.txt","r") as f:
        lines = f.readlines()
    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
    pos = util.read_pos(df)
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = Parse.adjustPosCamera(df, pos, ts)
    position = np.column_stack((ts, posx, posy))
    position = position[(position[:,0]>=start) & (position[:,0]<=end)]
    position = position[np.logical_or(position[:,1]>0, position[:,2]>0)]
    position = Filt.velocity_filtering(position, vthresh)
    clust = np.asarray(util.read_clust(df+clustName))
    clust = Filt.unitVelocityFilter(ts, position, clust)
    clust = clust[(clust >= start) & (clust <= end)]
    spikexy = util.getPosFromTs(clust,position)
    spikes = np.column_stack((clust,spikexy))
    
    unit = Unit(spikes, position, clustName, smoothing)
    unit = Filt.filterFields(unit)
    
    return unit

class Unit():
    """
    Wrapper class because rep field ID algorithm looks
    for instance.spikes and instance.position
    """
    
    def __init__(self, s, p, clustname, smoothing):
        self.name = clustname
        self.spikes = s
        self.position = p
        self.fields = []
        self.visits = [] # nested list. each list is a subfield and values are themselves lists of points in visit
        self.perimeters = [] # has had the convex alg run on it
        self.colors = cnames
        self.smoothing = smoothing
        self.repUnit = RateMapClass.RateMap(self) # a different unit class from the pf alg someone else wrote
        # self.repUnit.PF is a list of pf objects. each object has pf.perimeter as an attribute which is the well-fitting but out of order [x,y] border lists
        self.alphaHullFactor = 1
        self.alpha = 0
        self.findFields()
        
    def PolyArea(self,x,y):
        """
        Found at https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
        """
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
    def findFields(self):
        self.fields = []
        self.visits = []
        for i,pf in enumerate(self.repUnit.PF[:]):
            self.visits.append([])
            #create boundary using alphahull alg which allows for concave hulls but does not work all that well as is
#             alpha = self.alphaHullFactor*alphashape.optimizealpha(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)))
#             hull = alphashape.alphashape(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)),alpha)
#             hxy = hull.exterior.coords.xy
#             contour = path.Path(list(zip(hxy[0],hxy[1])))

            #create boundary using convex hull
            points = np.asarray(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2))) # use this to go from 2d hist coords to camera coords
            hull = ConvexHull(points)
            vertices = np.append(hull.vertices, hull.vertices[0]) # add the first point to close the contour
            contour = path.Path(points[vertices])

            PinC = self.position[contour.contains_points(self.position[:,1:])]
            posVisits = getVisits(PinC[:,0])
            self.visits[-1].append(posVisits)
            field_FR = []
            field_TS = [] # take middle ts 
            for visit in posVisits:
                spk = self.spikes[np.logical_and(self.spikes[:,0] > visit[0], self.spikes[:,0] < visit[-1])]
                vdur = (visit[-1]-visit[0])/1e6 # do this instead of just multiplying number samples by fr because any lost frames from, e.g. occlusions will introduce a discrepancy
                field_FR.append(spk.shape[0]/vdur)
                field_TS.append(visit[0])
                
            field_FR = util.weird_smooth(np.asarray(field_FR), self.smoothing)
            
            totalRate = sum(field_FR)
            area = self.PolyArea(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)
            print(totalRate/area)
            if True:
                self.fields.append(np.column_stack((field_TS, field_FR)))
                self.perimeters.append(points[vertices])
                
def getVisits(data, maxgap=2*1e6):
    """
    """
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

cnames = ['black', 'blue', 'green', 'red', 'brown', 'purple', 'cornflowerblue', 'orchid', 'darkcyan', 'midnightblue', 'saddlebrown', 'darkviolet', 'seagreen', 'indianred', 'goldenrod', 'orange', 'olive']
binWidth = wmDef.binWidth
from matplotlib import path
cmap = util.makeCustomColormap()

def plotRoutine_RepPF_TempDyn(unit,smoothing=0, nf=99, time='time', save=False, savepath=[]):
    """
    Plotting routine using unit and ratemapclass classes in local namespace
    Create a ratemap of fields with detected fields outlined
    Below, plot visit FR over time for each field, color coded to match
    
    """
    clust = unit.name
    fig, ax = plt.subplots(2,1, figsize=(10,14))
    
    fig.axes[0].imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                       cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
    fig.axes[0].set_title(f"{clust}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=20)
    fig.axes[0].axis('equal')
    fig.axes[0].set_ylim([0,480])
    fig.axes[0].set_xlim([0, 640])
    fig.axes[0].set_xticks([])
    fig.axes[0].set_yticks([])
    fig.axes[0].spines['top'].set_visible(False)
    fig.axes[0].spines['right'].set_visible(False)
    fig.axes[0].spines['bottom'].set_visible(False)
    fig.axes[0].spines['left'].set_visible(False)

    #option to not visualize garbage fields. 99 is an 'infinity' value as no cell will have 99 fields.
    if nf == 99:
        end = len(unit.fields)
        contig=range(end)
    elif type(nf) == list:
        end = len(unit.fields)
        contig = nf
    else:
        pass
        
    for i, field in enumerate(unit.fields[:end]):
        f = util.weird_smooth(field[:,1],smoothing)
        
        if time == 'time':
            xval = field[:,0]
        elif time  == 'visit' or time == 'visits':
            xval = range(field.shape[0])
        
        fig.axes[1].plot(xval, f, color=unit.colors[i], marker='.',alpha=0.8)
        fig.axes[0].plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i])
        fig.axes[1].tick_params(axis='y', labelsize=14)
        fig.axes[1].tick_params(axis='x', labelsize=14)
        fig.axes[1].set_xlabel(f"Time in session ({(lambda x: 'min' if x == 'time' else 'visits')(time)})", fontsize=24)
        fig.axes[1].set_ylabel("Firing Rate (Hz, smoothed)", fontsize=24)
        fig.axes[1].spines['right'].set_visible(False)
        fig.axes[1].spines['top'].set_visible(False)
        fig.axes[1].set_title(f"gaussian smoothing sigma = {smoothing+unit.smoothing}", fontsize=12)
    
    if save:
        clustname = clust.replace("\\","_")
        plt.savefig(fname=savepath+clustname+".png", format='png')
        plt.close()
        
def interpolateField(unit,wnSize=5, wnStep=2,s=5,k=3,plot=True, ret=False):
    """
    Input: A unit object with attribute fields (list of [ts,fr] arrays])
           wnSize - size of sliding window in time units (minutes)
           wnStep - shift of sliding window in time units (minutes)
           s - smoothing of spline
           k - degree of spline
           
    Diagnostic function to view spline interpolation performance. Interp is also done in the analysis itself. 
    
    Uses scipy.interpolate.splrep and splev to create an interpolated representation
    of each field. This is because field visits are asynchronous (a rat cannot be in two places at once)
    and so to compare activity between fields at a given time point we need to interpolate (unless time window is large)
    
    The spline rep is contructed in each sliding window with the params above
    Returns (optional) xs, ys which are lists of the interpolated points (w overlap bc sliding window)
    """
    fmax = int(np.ceil(max([max(field[:,0]) for field in unit.fields])))
    wins = []
    for i in range(fmax):
        a,b = 0+(i*wnStep), wnSize+(i*wnStep)
        if b < np.ceil(fmax):
            wins.append((a,b))
    # For each field, get the spline params so two fields can be interpolated to same # pts within a window, allowing for a pearson R calc
    fieldFx = [splrep(d[:,0], d[:,1], k=k, task=0, s=s) for d in unit.fields]
    ## sample spline fx in wins as would be done in analysis and view
    nf = len(unit.fields)
    xs, ys = [], []
    for j in range(nf):
        xc, yc = [], []
        for w in wins:
            start, end = w
            x = np.linspace(start, end, 100)
            interp= splev(x, fieldFx[j])
            xc.append(x)
            yc.append(interp)
        xs.append(xc)
        ys.append(yc)
    for i in range(nf):
        xs[i] = [item for sublist in xs[i] for item in sublist]
        ys[i] = [item for sublist in ys[i] for item in sublist]
    if plot:
        plt.figure()
        for i,c in zip(range(nf),['b','r','g','k']):
            plt.plot(xs[i], ys[i],'.',color=c,markersize=4)
            plt.plot(unit.fields[i][:,0], unit.fields[i][:,1], color=c)
        plt.title(unit.name)
    if ret:
        return xs, ys
    
def plotMats(clust, mats, wins, mattype='diff',shuff='False',vthresh=[]):
    ncol=10
    fig, ax = plt.subplots(int(np.ceil(len(mats)/ncol)),ncol,figsize=(8,8))
    if vthresh == []:
        _max, _min = max([arr.max() for arr in mats]), min([arr.min() for arr in mats])
    else:
        _max, _min = vthresh[1], vthresh[0]
    
    
    for i in range(len(mats)):
        im = fig.axes[i].imshow(mats[i], aspect='auto',interpolation='None', cmap=cmap, vmin=_min, vmax=_max)
        fig.axes[i].set_title(f"Mins {wins[i][0]}-{wins[i][1]}")
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    for i in range(len(fig.axes)):
        fig.axes[i].set_xticks([])
        fig.axes[i].set_yticks([])
    plt.suptitle(f"{clust} Time Varying Unsigned Mean {mattype} Matrices, shuff = {shuff}", fontsize=20)
    
def plotCorrOfMats(mats):
    corrOfcorrs = np.empty((len(mats),len(mats)))

    for i in range(len(mats)):
        for j in range(len(mats)):
            corrOfcorrs[i,j] = scipy.stats.pearsonr(mats[i].flatten(),mats[j].flatten())[0]
        
    plt.figure()
    plt.imshow(corrOfcorrs, origin='lower', interpolation='None')
    plt.title("Autocorrelation Matrix of Difference Matrices",fontsize=22)
    plt.xlabel("Difference matrices",fontsize=16)
    plt.ylabel("Difference matrices", fontsize=16)
    
def makeSemaphores(fieldArray):
    # Semaphore Plot Analysis
    s=1
    k=3 # should be 3 usually
    fieldFx = [splrep(d[:,0], d[:,1], k=k, task=0, s=s) for d in fieldArray]
    fmax = int(np.ceil(max([max(field[:,0]) for field in fieldArray])))
    wnSize=5
    wnStep = 2
    wins = []
    for i in range(0,fmax):
        a,b = 0+(i*wnStep), wnSize+(i*wnStep)
        if b < np.ceil(fmax):
            wins.append((a,b))

    nf = len(fieldArray)
    diffmats = []
    for w in wins:
        start, end = w
        diffmat = np.zeros((nf,nf))
        for i in range(nf):
            for j in range(nf):
                x = np.linspace(start, end, 100)
                ainterp, binterp = splev(x,fieldFx[i]), splev(x, fieldFx[j])
                diff = np.mean(ainterp)-np.mean(binterp)
                diffmat[i,j] = diff
        diffmats.append(diffmat)
    diffmats = np.asarray(diffmats)
    return diffmats