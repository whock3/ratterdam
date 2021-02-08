# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 21:23:40 2019

@author: whockei1
"""

# ML libraries
import sklearn as skl
import matplotlib.pyplot as plt
from sklearn import svm, preprocessing, metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score
from sklearn.decomposition import PCA


#Mostly General scientific computing libraries, some general purpose libraries
import numpy as np 
from datetime import datetime
from scipy.stats import sem
from scipy.ndimage import center_of_mass
from scipy.integrate import simps
from importlib import reload
import argparse, multiprocessing

# My libraries
import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_CoreDataStructures as Core
import ratterdam_PermutationTests as Perm
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt

from itertools import cycle

# Variables
beltwayAlleys = [16, 17, 3, 1, 5, 7, 8, 10, 11] # beltway alley IDs in terms of their full track, 17-alley ID
nbins = Def.singleAlleyBins[0]-1


def checkRM(ratemap):
    """Utility function to take a 1-d linear ratemap and see if it is valid.
    Dec 2019: If it's not all nans it's valid. Ie. don't ignore any passes
    by virture of lack of firing (as long as there is sampling)
    """
    if type(ratemap) == np.ndarray and np.nanmax(ratemap) >= 0:
        return True
    else:
        return False
    
def calcStatPair(unit, alley, pair):
    x,y = pair #unpack the pair of stimuli eg AB
    rma, rmb = unit.linRMS[alley][x], unit.linRMS[alley][y]
    stat = np.abs(rma-rmb)
    return stat

def generateLabel(target, alley, stimulus):
    if target == 'Alley':
        label = str(alley)
    elif target == 'Stimulus':
        if len(stimulus) == 1:
            label = stimulus
        elif len(stimulus) == 2:
            label = stimulus[0] # if label is a pair of stims for a test stat label
    elif target == 'AlleyXStimulus':
        if len(stimulus) == 1:
            label = f"{alley}{stimulus}"
        elif len(stimulus) == 2:
            label = f"{alley}{stimulus[0]}"
    
    return label

def calcRMdeviation(unit, alley, visit, window=5, deviation='diff'):
    """
    Compute the deviation of a visit's ratemap
    for a cell from the average within a window
    of all visits' rms. Window var is one-sided
    so whole window size is 2*window. 
    
    Deviation is 
    - 'diff' for binwise difference
    """
                    
    rmvisit = unit.alleys[alley][visit]['ratemap1d']
    rms_in_window = np.empty((0, nbins))
    
    if visit < window:
        for i in range(0, (window-visit)+1):
            rm = unit.alleys[alley][i]['ratemap1d']
            rms_in_window = np.vstack((rms_in_window, rm))
            
    elif (visit+window) >= len(unit.alleys[alley]):
        for i in range(visit-window, visit+(len(unit.alleys[alley])-visit)):
            rm = unit.alleys[alley][i]['ratemap1d']
            rms_in_window = np.vstack((rms_in_window, rm))     
            
    else:
        for i in range(visit-window, visit+window+1):
            rm = unit.alleys[alley][i]['ratemap1d']
            rms_in_window = np.vstack((rms_in_window, rm))
            
    mask = np.ma.masked_invalid(rms_in_window)
    avg = np.median(mask.data,axis=0) # ignores inf and nan
    if deviation == 'diff':
        stat = np.abs(rmvisit-avg)
    if deviation == 'grad': 
        stat = np.gradient(rmvisit)
    elif deviation.lower() == 'none':
        stat = rmvisit # this passthrough is so you dont have to change anything but a toggle in setup fx

    else:
        return "Incorrect deviation argument"
    
    return stat
        
def calcRMSummary(unit, alley, visit):
    """
    For a given pass for a given unit
    summarize the 1d ratemap into a simpler,
    explicit vector of attributes
    - period within session (half, third?,)
    - max
    - min
    - mean
    - loc max
    - loc min
    """
    rm = unit.alleys[alley][visit]['ratemap1d']
    epoch = np.digitize(visit,[0,30,60])
    maximum, minimum, mean = np.nanpercentile(rm, 98), np.nanmin(rm), np.nanmean(rm)
    locmax, locmin = np.nanargmax(rm), np.nanargmin(rm)
    auc = simps(rm)
    avgdrds = np.mean(np.abs(np.diff(rm))) # avg dr/ds change in rate / change in pos. 
    maxdrds = np.max(np.abs(np.diff(rm)))
    com = center_of_mass(rm)
    return np.array((epoch, maximum, minimum, mean, locmax, locmin, auc, avgdrds, maxdrds, com[0]))

def setupRFdata(target, population, shuffle=False):

    """
    This block collects data into data matrix X and label matrix Y.
    Each entry is an averaged linear ratemap for a cell under a given condition (alley/txt).
    Entry is the test statistic comparing AvB (labeled A), B vs C (labeled B) and CvA (labeled C)

    """

    X = np.empty((0, (10)*len(population)))
    Y = []


    for alley in beltwayAlleys:

        for visit in range(len(population[list(population.keys())[0]].alleys[alley])):

            stim = population[list(population.keys())[0]].alleys[alley][visit]['metadata']['stimulus']
            label = generateLabel(target, alley, stim)
            dataRow = np.empty((0))
            
            for unitname, Unit in population.items():
                #rm = calcStatPair(Unit, alley, pair)
                #rmA = calcRMdeviation(Unit,alley,visit, deviation='diff')
                rm = calcRMSummary(Unit, alley, visit)
                #rm = np.concatenate((rmA,rmB))
                if checkRM(rm)== True:
                    #rm = np.hstack((unitID, rm))
                    dataRow = np.hstack((dataRow, rm))
                else:
                    dataRow = np.hstack((dataRow, np.zeros((nbins)) ))


            Y.append(label)
            X = np.vstack((X, dataRow))


    X[np.where(~np.isfinite(X))] = 0
    X = preprocessing.StandardScaler().fit_transform(X)


    if shuffle is True:
        if target == 'AlleyXStimulus':
            txt = [i[-1] for i in Y]
            np.random.shuffle(txt)
            a = [i[:-1] for i in Y]
            Y = [f"{x}{y}" for x,y in zip(a,txt)]
        elif target == 'Stimulus':
            np.random.shuffle(Y)
        elif target == 'Alley':
            np.random.shuffle(Y)

    return X, Y

def runRandomForest(X, Y, runs=500, trees=1500, avgType='macro'):
    oobs, precisions, recalls, f1s, accuracies = [], [], [], [],[]

    for i in range(runs):
        print(i)
        clf = RandomForestClassifier(n_estimators=trees, oob_score=True)
        Xtrain, Xtest, ytrain, ytest = train_test_split(X,Y,shuffle=True,random_state=0)
        clf.fit(Xtrain,ytrain)
        oobs.append(clf.oob_score_)
        yfit = clf.predict(Xtest)
        p = precision_score(ytest, yfit, average=avgType)
        r = recall_score(ytest, yfit, average=avgType)
        f1 = f1_score(ytest, yfit, average=avgType)
        acc = accuracy_score(ytest,yfit)
        precisions.append(p)
        recalls.append(r)
        f1s.append(f1)
        accuracies.append(acc)
    return oobs, precisions, recalls, f1s, accuracies


if __name__ == '__main__':
    
    print("Running Random Forest Decoding Routine")
    parser = argparse.ArgumentParser()
    parser.add_argument("Xalley")
    parser.add_argument("Yalley")
    
    parser.add_argument("Xstim")
    parser.add_argument("Ystim")
    
    parser.add_argument("Xalleyxstim")
    parser.add_argument("Yalleyxstim")
    
    args = parser.parse_args()
    
    pool = multiprocessing.Pool(processes=3)
    data = pool.starmap(runRandomForest, zip([args.Xalley, args.Xstim, args.Xalleyxstim],[args.Yalley, args.Ystim, args.Yalleyxstim]))
    
    