# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 21:25:28 2021

@author: whockei1

Script to decode temporal information from the field(s) of neurons recorded
in Ratterdam Open. 

This code was in ratterdam_PFRepetition_TemporalDecoding1.ipynb

From summer 2020 thru Feb 2021 the decoding used a kNN classifier to
decode time from the inter-field rate difference (IFD) matrices. See thesis
ch4 decoding section on details. 

Decoding parameters have not, as of file creation, been totally settled on
but usually decoding session time in thirds, with 60:40 test:train (yes test
has more data).  

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
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
from matplotlib import path
from matplotlib.backends.backend_pdf import PdfPages
import more_itertools, itertools
from sklearn.metrics import auc
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from scipy.interpolate import splrep, splev
from scipy.spatial import ConvexHull
import scipy

class GuttedUnit():   
    def __init__(self):
        self.fields = []

def kNN_semaphore_decoding(fieldArray, reps=50, n=2, epochs=3, stat='median', trainProp=0.67, wnSize=5*1e6*60, wnStep=2*1e6*60):
    """
    kNearest Neighbor Decoding
    input:  data - Unit class object 
            reps -  # decoding repeats
            n - num neighbors for kNN
            epochs - number of temporal epochs to divide data into to decode
            stat - measure of decoding perf for each set of decoding runs (# = reps).
                    choices are 'median' or 'mean'
            trainProp - proportion of data to use as train (1-trainProp is test)
            wnSize - (for semaphore plot) sliding window size to create a single flag plot. unit in mins.
            snStep - (for semaphore plot) sliding window slide amount. unit in mins.
            fieldType - "smoothed" or "raw". Option chooses whether to use smoothed or unsmoothed field attrs of Unit obj
                    
    Takes Unit class object. Creates series of semaphore plots using makeSemaphores()
    Then decodes temporal information from these. Semaphores divided into epochs,
    ntrain proportion used to find nearest neighbors of test in euclidian space and see
    if epoch labels match. Do multiple times and report a statistic of perf
    
    returns: scalar of decoding perf, as defined by 'stat'
    """
        
    mats,_ = RepCore.makeSemaphores(fieldArray, wnSize=wnSize, wnStep=wnStep)
    flatmats = np.asarray([i.flatten() for i in np.asarray(mats)])
    idx = np.ceil(np.linspace(0,len(flatmats),epochs+1)).astype(np.int)
    labels = []
    for i in range(len(idx)-1):
        l= idx[i+1]-idx[i]
        labels.extend([i]*l)
    labels = np.asarray(labels)
    
    perfs = []
    for rep in range(reps):
        train, test = train_test_split(list(range(len(flatmats))), train_size=trainProp)
        Xtrain, Xtest, Ytrain, Ytest = flatmats[train], flatmats[test], labels[train], labels[test]
        neigh = KNeighborsClassifier(n_neighbors=n)
        neigh.fit(Xtrain, Ytrain)
        perfs.append(sum(Ytest==neigh.predict(Xtest))/len(Xtest))  
        
    if stat == 'median':
        sumStat = np.median(perfs)
    elif stat == 'mean':
        sumStat = np.mean(perfs)
    else:
        print("Invalid statistic supplied. Choose 'median' or 'mean' ")
        sumStat = None
    
    return sumStat

def rotateFields(fields):
    """
    Input   - unit: Unit class object      
    selects a temporal shift for each field. Shift is rotational such that 
    whatever 'falls off' end of vector 'comes back' to the front.
    Returns - GuttedField object with shifted fields in gunit.fields attribute
    """
    gunit = GuttedUnit()
    sfields = []
    for field in fields:
        nvisits = field.shape[0]
        shift = np.random.randint(0,nvisits)
        sf = np.column_stack((field[:,0], np.roll(field[:,1], shift)))
        sfields.append(sf)
    gunit.fields = sfields
    return gunit

def shuffleFields(fields):
    """
    Input   - unit: Unit class object   
    Shuffle visits within a field.
    Returns - GuttedField object with shifted fields in gunit.fields attribute
    """
    gunit = GuttedUnit()
    newfields = []
    for field in fields:
        newfields.append(np.column_stack((field[:,0],np.random.permutation(field[:,1]))))
    gunit.fields = newfields
    return gunit


if __name__ == '__main__':

    rat = "R808"
    day = "D6"
    savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
    df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
    clustList = util.getClustList(df)
    population = {}
    for clust in clustList:
        try:
            print(clust)
            unit = RepCore.loadRepeatingUnit(rat, day, clust, smoothing=1)
            rm = util.makeRM(unit.spikes, unit.position)
            if np.nanpercentile(rm, 95) > 1.:
                population[clust] = unit
                print(f"{clust} included")
            else:
                print(f"{clust} is not included")
        except:
            pass
        
        # Setup parameters and setup data
    parmdict = {
        'epochs':3, # number of temporal epochs to decode
        'nshuffles':1000, # shuffles creating new data each time
        'nreps':10, # repeated decoding given a set of data
        'trainProp':0.4, # test is the complement
        'nNeighbors':2,
        'wnSize':1*1e6*60,
        'wnStep':1*1e6*60,
        'smoothing':1,
        'stat':'median',
        'fieldRateThresh':0.3,
        'fieldpctThresh':20
        
    }
    data = {}
    ts = util.genTimestamp()
    for clust in population.keys():
        print(clust)
        try:
            unit = population[clust]
            shuffunit = RepCore.loadRepeatingUnit(rat, day, clust, smoothing=parmdict['smoothing'])
    
    
            realperf = kNN_semaphore_decoding(unit.smoothedFields, reps=parmdict['nreps'], n=parmdict['nNeighbors'], epochs=parmdict['epochs'], stat=parmdict['stat'], 
                                              trainProp=parmdict['trainProp'], wnSize=parmdict['wnSize'], wnStep=parmdict['wnStep'])
    
            rsperfs = [] # each data point is the median or mean (depending on user choice) of multiple decodings of a single shuffle
            vsperfs = []
            for n in range(parmdict['nshuffles']):
                shuffunit.shuffleFieldsFx()
                shuffunit.smoothFieldsFx()
    
                #rgunit = rotateFields(unit.smoothedFields)
                #sperf = kNN_semaphore_decoding(rgunit.fields, reps=parmdict['nreps'], n=parmdict['nNeighbors'], epochs=parmdict['epochs'], stat=parmdict['stat'], 
                   #                           trainProp=parmdict['trainProp'], wnSize=parmdict['wnSize'], wnStep=parmdict['wnStep'])
                #rsperfs.append(sperf)
    
                sperf = kNN_semaphore_decoding(shuffunit.smoothedFields, reps=parmdict['nreps'], n=parmdict['nNeighbors'], epochs=parmdict['epochs'], stat=parmdict['stat'], 
                                              trainProp=parmdict['trainProp'], wnSize=parmdict['wnSize'], wnStep=parmdict['wnStep'])
                vsperfs.append(sperf)
            if realperf > np.percentile(vsperfs, 95):
                data[unit.name] = {'real':realperf, 'null95':np.percentile(vsperfs, 95), 'pass':True}
            else:
                data[unit.name] = {'real':realperf, 'null95':np.percentile(vsperfs, 95), 'pass':False}
        except:
            print(f"{clust} did not run successfully")