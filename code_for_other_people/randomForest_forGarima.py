# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 17:48:11 2020

@author: whockei1
"""

# Imports 
import sklearn as skl
import matplotlib.pyplot as plt
from sklearn import svm, preprocessing, metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score
from sklearn.metrics import auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
from scipy.ndimage import center_of_mass
import numpy as np, random, json, pickle, datetime, copy, socket, os, sys, scipy
from scipy.stats import sem
import matplotlib.colors as colors
from importlib import reload

def createVisitSummaryFeatures(unit, alley, visit, features):
    """
    For a given pass for a given unit summarize the 1d ratemap into a simpler,
    explicit vector of attributes. Which attributes to use are given by
    the 'features' list. Visit is the visitnum not the data itself
    """
    feats = np.empty((0))
    rm = unit.alleys[alley][visit]['ratemap1d']
    # i dont know of a better way of doing this other than to just check param name in list and add it if present
    if 'rm' in features:
        feats = np.append(feats, rm)
    if 'time' in features:
        feats = np.append(feats, visit)
    if 'max95' in features:
        maximum = np.nanpercentile(rm, 95)
        feats = np.append(feats, maximum)
    if 'locmax95' in features:
        locmax = np.searchsorted(np.sort(rm), np.percentile(rm, 95))
        feats = np.append(feats, locmax)
    if 'mean' in features:
        mean = np.nanmean(rm)
        feats = np.append(feats, mean)
    if 'median' in features:
        feats = np.append(feats, np.nanmedian(rm))
    if 'auc' in features:
        auc = simps(rm)
        feats = np.append(feats, auc)
    if 'avgdrds' in features:
        avgdrds = np.mean(np.abs(np.diff(rm))) # avg dr/ds change in rate / change in pos. 
        feats = np.append(feats, avgdrds)
    if 'maxdrds' in features:
        maxdrds = np.percentile(np.abs(np.diff(rm)), 95)
        feats = np.append(feats, maxdrds)
    if 'com' in features:
        try:
            com = center_of_mass(np.nan_to_num(rm))[0]
            feats = np.append(feats, com)
        except:
            com = int(round(Def.singleAlleyBins[1]-1)/2)
            feats = np.append(feats, com)
    if 'comval' in features:
        try:
            comval = rm[int(np.round(com))]
            feats  = np.append(feats, comval)
        except:
            comval = np.nanpercentile(rm, 95)
            feats  = np.append(feats, comval)
    if 'boundMaxs' in features:
        # think there may be something going on at entrace/exit to alley so get the max val 
        # within each trisector of alley. NB real intersection ends with alleybounds_manuallyshifted2
        # gives a 3:6:3 ratio of approach L:alley:approach/exit R but I want to squeeze in the bounds to give
        # more space to the flanks (4:4:4 ratio) to capture whats happening at boundary itself as well.
        max1,max2,max3 = np.nanmax(rm[:4]), np.nanmax(rm[4:8]), np.nanmax(rm[9:]) # assumes 12bin rm. make generalized later
        feats = np.append(feats,(max1,max2,max3))
    if 'isi' in features:
        begin, end, bsize = 0, 0.075e6, 5000
        bins = np.arange(begin,end,bsize)
        spikes = unit.alleys[alley][visit]['spikes']
        hist = np.histogram(np.diff(spikes[:,0]),bins=bins)[0]
        feats = np.append(feats, hist)
        
    if 'gamma_params':
        rm = rm[~np.isnan(rm)]
        rm = rm/rm.sum()
        a,loc,b = scipy.stats.gamma.fit(rm)
        feats = np.append(feats, [a,loc,b])
        
    return feats

def parmString(parmdict, features):
    string = ''
    for k,v in parmdict.items():
        string += f"{k}:{v}\n"
    for f in features:
        string +=f"{f}\n"
    return string

def setupAlleyData(unit, alley, repFx, features):
    """
    Create a matrix (n,b) where n is the number of trials at that alley 
    (usually with rewards removed, but that's done in the unit.loadRawData fx)
    and b are the number of spatial bins in the 1d ratemap of each trial at that alley
    """
    X = [] # dont know how size of feature vec ahead of time so array it later
    Y = np.empty((0))
    
    for visitNum,visit in enumerate(unit.alleys[alley]):
        reprm = repFx(unit, alley, visitNum, features)
        X.append(reprm)
        Y = np.append(Y, unit.alleys[alley][visitNum]['metadata']['stimulus'])
    
    X = np.asarray(X)
    X[np.where(~np.isfinite(X))] = 0
    #X = preprocessing.StandardScaler().fit_transform(X)
    
    return X, Y


def runRandomForest(X, Y, parmdict):
    oobs = []
#     fimps = []
#     paths = []
    for i in range(parmdict['nRuns']):
        clf = RandomForestClassifier(n_estimators=parmdict['nTrees'], 
                                     oob_score=True,
                                     max_features = parmdict['Max Feats'],
                                     max_depth = parmdict['Max Depth']
                                    )       
        clf.fit(X,Y)
        oobs.append(clf.oob_score_)
#         fimps.append(clf.feature_importances_)
#         paths.append(clf.decision_path(X))
        
    return oobs

# Load data into population dict. Each cell will be decoded separately. Within each cell each alley will be decoded separately.
rat = 'R859'
expCode = "BRD5"
datafile = f"E:\\Ratterdam\\{rat}\\{rat}{expCode}\\"

alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
population = {}
for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            unit = Core.UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
            unit.loadData_raw()
            rm = util.makeRM(unit.spikes, unit.position)            
            if np.nanpercentile(rm,Def.wholetrack_imshow_pct_cutoff) >= 1.:
                print(clustname)
                population[unit.name] = unit
                

parmdict = {
    'nRuns':3, # reps of decoding for a given dataset. tech replicates. 
    'nTrees':10000, 
    'Max Depth':None, 
    'Max Feats':'auto',
    'Cell inclusion in population': '95pctile >=1Hz overall',
    'Visit inclusion in data matrix': '12x visits Mean alley activity >=1 Hz, 3 contig bins >=20% max avg field',
    'Bootstrap': 'None',
    'Shuffle': 200
    }

features = [
    'rm'
           ]


clust = 'TT6cl-maze1.8'
unit = population[clust]

string = parmString(parmdict, features)
stamp = util.genTimestamp()

fig, axes = plt.subplots(3,3,figsize=(10,10))
plt.suptitle(f"{clust} RF Decoding by Alley")
plt.text(0.005, 0.2, string, fontsize=6, transform=plt.gcf().transFigure)

for i,alley in enumerate(Def.beltwayAlleys):
    valid = Filt.checkMinimumPassesActivity(unit, alley)
    if valid:
        ax = fig.axes[i]
        print(alley)
        X,Y = setupAlleyData(unit, alley, createVisitSummaryFeatures, features)
        realoobs = runRandomForest(X,Y, parmdict)
        realmean = np.mean(realoobs)
        nulloobs = np.zeros((0,1))
        for i in range(parmdict['Shuffle']):
            Ys = np.random.permutation(Y)
            ssoobs = runRandomForest(X,Ys, parmdict)
            nulloobs = np.vstack((nulloobs, np.mean(ssoobs)))
        
        ax.hist(nulloobs)
        ax.vlines(np.percentile(nulloobs, 95),0,150,'k')
        ax.vlines(realmean,0,100,'r')
        ax.set_title(f"Alley {alley}, mean {realmean} vs {np.percentile(nulloobs, 95)} 95% null pct")
    else:
        fig.axes[i].set_title(f"Insufficient Activity in Alley {alley} for Decoding")
