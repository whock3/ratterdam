# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 18:27:43 2020

@author: whockei1
"""


import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import pandas as pd
from scipy.stats import sem
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import placeFieldBorders
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
from matplotlib.backends.backend_pdf import PdfPages
import more_itertools, itertools
from sklearn.metrics import auc
from scipy.interpolate import splrep, splev, PchipInterpolator as pchip
from scipy.spatial import ConvexHull
import scipy
import ratterdam_DataFiltering as Filt
import alleyTransitions as alleyTrans
import newAlleyBounds as nab


def loadRepeatingUnit(df, clustName, smoothing=2, vthresh=Def.velocity_filter_thresh):

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
    
    # On R886 D1 the track or camera got moved wrt where it was on D2 
    # which was the day used to make the standard track bound coords. -20 x fixes it
    day = df.split('\\')[-2] # last element is ''
    if day.split("_")[0] == 'R886' and day.split("_")[2] == 'D1':
        print("Readjusting position for R886 D1")
        position[:,1] = position[:,1] - 20
    
    clust = np.asarray(util.read_clust(df+clustName))
    clust = Filt.unitVelocityFilter(ts, position, clust)
    clust = clust[(clust >= start) & (clust <= end)]
    spikexy = util.getPosFromTs(clust,position)
    spikes = np.column_stack((clust,spikexy))
    

    unit = Unit(spikes,position, clustName, smoothing)
    unit = Filt.filterFields(unit)
    if smoothing:
        unit.smoothFieldsFx()
    
    return unit

def unitVelocityFilter(ts, position, clust):
        '''
        Remove spikes emitted when the animal's
        velocity was below threshold. 
        
        This class does not explicitly store
        the thresh, or even whether a velocity
        filter was applied. That occurs in the BehavioralData
        class and is reflected implicitly in the position array
        having been filtered. the ts dict is unfiltered and is used
        as a reference.
        '''
        
        allSpikeTs = np.asarray([util.takeClosest(ts, i) for i in clust])
        filtTs = clust[np.isin(allSpikeTs, position[:,0])]
        return filtTs

class Unit():
    """
    Wrapper class because rep field ID algorithm looks
    for instance.spikes and instance.position
    """
    
    def __init__(self, s, p, clustname, smoothing):
        self.name = clustname
        self.spikes = s
        self.position = p
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
    
    
    def filterVisits(self, visits, perim,thresh=0.2):
        """
        Take a list of visits. Each element is a list of ts corresponding to the
        pos samples from that visit.
        
        From there, apply a threshold to remove visits that are too short.
        7-26-2021 thresh will be based on max dist traveled thru field and size
        of field
        
        Thresh is 0-1 pct max dist a visit must travel
        """
        #calc bounding box and max distances
        maxx,minx, maxy, miny = max(perim[:,0]), min(perim[:,0]), max(perim[:,1]), min(perim[:,1])
        dista = np.sqrt(((minx-minx)**2)+((maxy-miny)**2))
        distb = np.sqrt(((minx-maxx)**2)+((miny-miny)**2))
    
        newvisits = []
        for v in visits:
            if len(v) >= 2:
                traj = self.position[(self.position[:,0]>v[0])&(self.position[:,0]<=v[-1])]
                #below,compute distance between 1st point and all others to get max dist inside field
                maxdist = max(np.linalg.norm(traj[0,1:]-traj[:,1:],axis=1))
                if maxdist >= thresh*dista or maxdist >= thresh*distb:
                    newvisits.append(v)
       
        return newvisits
    
    def shuffleFieldsFx(self):
        for i,field in enumerate(self.fields):
            self.fields[i] = np.column_stack((self.fields[i][:,0], np.random.permutation(self.fields[i][:,1])))
            
    def smoothFieldsFx(self):
        self.smoothedFields = []
        for i, field in enumerate(self.fields):
            self.smoothedFields.append(np.column_stack((self.fields[i][:,0], util.weird_smooth(self.fields[i][:,1],self.smoothing))))
        
    def findFields(self):
        self.fields = []
        self.visits = []
        self.perimeters = []
        for i,pf in enumerate(self.repUnit.PF[:]):
            border = placeFieldBorders.reorderBorder(pf.perimeter, i)
            self.perimeters.append(border)
            #border = np.append(border, border[0]) # add the first point to close the contour
            contour = path.Path(border)

            PinC = self.position[contour.contains_points(self.position[:,1:])]
            posVisits = getVisits(PinC[:,0])
            posVisits = self.filterVisits(posVisits, border, thresh=0.2)

            
            self.visits.append(posVisits)
            field_FR = []
            field_TS = [] # take middle ts 
            for visit in posVisits:
                spk = self.spikes[np.logical_and(self.spikes[:,0] > visit[0], self.spikes[:,0] < visit[-1])]
                vdur = (visit[-1]-visit[0])/1e6 # do this instead of just multiplying number samples by fr because any lost frames from, e.g. occlusions will introduce a discrepancy
                field_FR.append(spk.shape[0]/vdur)
                field_TS.append(visit[0])
                
            # stack and strip nans
            finalfield = np.column_stack((field_TS, field_FR))
            finalfield = finalfield[~np.isnan(finalfield[:,1])]
            
            self.fields.append(finalfield)
            
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

def plotRoutine_RepPF_TempDyn(unit, nf=99, time='time', save=False, savepath=[]):
    """
    Plotting routine using unit and ratemapclass classes in local namespace
    Create a ratemap of fields with detected fields outlined
    Below, plot visit FR over time for each field, color coded to match
    
    """
    clust = unit.name
    fig, ax = plt.subplots(1,2, figsize=(17,5))
    
    fig.axes[0].imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                       cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
    fig.axes[0].set_title(f"{clust}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=20)
    #fig.axes[0].axis('equal')
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
        fieldlist = range(99)
    elif nf < 99:
        fieldlist = range(nf)
    elif type(nf) == list:
        fieldlist = nf
    else:
        pass
    
    if unit.smoothedFields:
        fields = unit.smoothedFields
    else:
        fields = unit.fields
        
    for i, field in enumerate(fields):
        if i in fieldlist:
            
            if time == 'time':
                xval = field[:,0]
                xval = [((i-field[0,0])/1e6)/60 for i in xval]
            elif time  == 'visit' or time == 'visits':
                xval = range(field.shape[0])
            
            fig.axes[1].plot(xval, field[:,1], color=unit.colors[i], marker='.',alpha=0.8)
            fig.axes[0].plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i])
            fig.axes[1].text(xval[0]-0.1,field[0,1]-0.1,i)
            fig.axes[1].tick_params(axis='y', labelsize=14)
            fig.axes[1].tick_params(axis='x', labelsize=14)
            fig.axes[1].set_xlabel(f"Time in session ({(lambda x: 'min' if x == 'time' else 'visits')(time)})", fontsize=18)
            fig.axes[1].set_ylabel(f"Firing Rate Hz (sigma = {unit.smoothing})", fontsize=18)
            fig.axes[1].spines['right'].set_visible(False)
            fig.axes[1].spines['top'].set_visible(False)
            fig.axes[1].set_title("Place Field Dynamics", fontsize=20)
    
    if save:
        clustname = clust.replace("\\","_")
        plt.savefig(fname=savepath+clustname+".png", format='png')
        plt.close()
        
def interpolateField(fields,name,wnSize=5*1e6*60, wnStep=2*1e6*60,s=5,k=3,plot=True, ret=False):
    """
    
    DEPRECATED: This uses a spline (splev/splrep) which is prone to under and
                overshooting real curve. Instead using a PCHIP spline interpolation
                which is less sensitive to this in the makeSemaphores() 
                Implemented using pchip in scipy around 9/23/20 (but this deprecation notice added 3/19/21)
    
    
    ----------------------------
    
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
    fmax = int(np.ceil(max([max(field[:,0]) for field in fields])))
    fmin = int(np.ceil(min([min(field[:,0]) for field in fields])))

    wins = []
    begin = fmin
    stop = False
    
    while not stop:
        a,b = begin, begin + wnSize
        if b < np.ceil(fmax):
            wins.append((a,b))
            begin += wnStep
        else:
            stop = True
            
            
    # For each field, get the spline params so two fields can be interpolated to same # pts within a window, allowing for a pearson R calc
    fieldFx = [splrep(d[:,0], d[:,1], k=k, task=0, s=s) for d in fields]
    ## sample spline fx in wins as would be done in analysis and view
    nf = len(fields)
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
        for i,c in zip(range(nf),cnames):
            plt.plot(xs[i], ys[i],'.',color=c,markersize=4)
            plt.plot(fields[i][:,0], fields[i][:,1], color=c)
        plt.title(name)
    if ret:
        return xs, ys
    
def plotMats(clustname, mats, wins, mattype='diff',vthresh=[],cm=None):
    ncol=10
    if cm is None:
        cm = cmap
    fig, ax = plt.subplots(int(np.ceil(len(mats)/ncol)),ncol,figsize=(8,8))
    if vthresh == []:
        _max, _min = max([arr.max() for arr in mats]), min([arr.min() for arr in mats])
    else:
        _max, _min = vthresh[1], vthresh[0]
    
    
    for i in range(len(mats)):
        im = fig.axes[i].imshow(mats[i], aspect='auto',interpolation='None', cmap=cm, vmin=_min, vmax=_max)
        fig.axes[i].set_title(f"Mins {wins[i][0]}-{wins[i][1]}")
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    for i in range(len(fig.axes)):
        fig.axes[i].set_xticks([])
        fig.axes[i].set_yticks([])
    plt.suptitle(f"{clustname} Time Varying Unsigned Mean {mattype} Matrices", fontsize=16)
    
def corrOfMats(mats, plot=True, ret=False):
    corrOfcorrs = np.empty((len(mats),len(mats)))

    for i in range(len(mats)):
        for j in range(len(mats)):
            corrOfcorrs[i,j] = scipy.stats.pearsonr(mats[i].flatten(),mats[j].flatten())[0]
    if plot:
        plt.figure()
        plt.imshow(corrOfcorrs, origin='lower', interpolation='None')
        plt.title("Autocorrelation Matrix of Difference Matrices",fontsize=22)
        plt.xlabel("Difference matrices",fontsize=16)
        plt.ylabel("Difference matrices", fontsize=16)
    if ret:
        return corrOfcorrs
    
def makeSemaphores(fieldArray,  wnSize=5*1e6*60, wnStep=2*1e6*60):
    # Semaphore Plot Analysis
    fmax = int(np.ceil(max([max(field[:,0]) for field in fieldArray])))
    fmin = int(np.ceil(min([min(field[:,0]) for field in fieldArray])))
    pchip_fields = [pchip(field[:,0], field[:,1]) for field in fieldArray] 

    wins = []
    begin = fmin
    stop = False
    
    while not stop:
        a,b = begin, begin + wnSize
        if b < np.ceil(fmax):
            wins.append((a,b))
            begin += wnStep
        else:
            stop = True
    

    nf = len(fieldArray)
    diffmats = []
    for w in wins:
        start, end = w
        diffmat = np.zeros((nf,nf))
        for i in range(nf):
            for j in range(nf):
                x = np.linspace(start, end, 100)
                ainterp, binterp = pchip_fields[i](x), pchip_fields[j](x)
                if len(fieldArray) > 1:
                    diff = np.abs(np.mean(ainterp)-np.mean(binterp))
                else:
                    # i.e. if there's only one field then the data point is just that fields rate in that window
                    diff = np.mean(ainterp)
                diffmat[i,j] = diff
        diffmats.append(diffmat)
    diffmats = np.asarray(diffmats)
    return diffmats, wins


def loadRecordingSessionData(rat,day):
    """
    Load one day of recording from a rat
    Return  - the population as a dict of Unit() objects (see Unitclass in this file RepCoreFx)
            - turn dataframe from alleyTransitions.py 
    Units loaded with qual thresh >= 3 and peak FR above 1Hz
    """
    ratborders = {'R781':nab.R781, 'R808':nab.R808, 'R859':nab.R859, 'R765':nab.R765,'R886':nab.R886}[rat]
    datapath = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    clustList, clustQuals = util.getClustList(datapath)
    population = {}
    qualThresh = 3
    
    for i,clust in enumerate(clustList):
        if clustQuals[i] >= qualThresh:
       
            print(clust)
            unit = loadRepeatingUnit(datapath, clust, smoothing=1)                                   
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
            
    pos, turns = alleyTrans.alleyTransitions(unit.position, ratborders, graph=False)
    turns = pd.DataFrame(turns)
    turns.columns = ['Allo-','Ego','Allo+','Ts exit','Ts entry', 'Alley-', 'Inter','Alley+']
    
    turns = pd.DataFrame(data=turns)
    turns.dropna(inplace=True) 
    
    return population, turns
        