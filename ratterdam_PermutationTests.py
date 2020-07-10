# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 19:44:34 2018

@author: whockei1
"""

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket, math
from scipy.stats import sem
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys, os, csv

if socket.gethostname() == 'Tolman':
    codeDirBase = 'C:\\Users\\whockei1\\Google Drive'
elif socket.gethostname() == 'DESKTOP-BECTOJ9':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'
    
sys.path.insert(0, codeDirBase + '\\KnierimLab\\Ratterdam\\Code')
import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_Defaults as Def
import ratterdam_CoreDataStructures as Core
import ratterdam_DataFiltering as Filt

def poolTrials(unit, alley, labels, txt):
    """
    Pool all trials that will form a group.
    Group defined as linear RM (computed differently from viz. lin rm)
    from all visits to a given alley when it harbored a given texture.
    
    This does not subsample to approx. balance group sizes. That is done after.
    
    Labels is a list of texture labels, either real or shuffled prior to this fx
    """
    rms = []
    idx = []
    visits = unit.alleys[alley]
    for i,visit in enumerate(visits):
        if labels[i] == txt:
            rm = visit['ratemap1d']
            if type(rm) == np.ndarray:
                rm = np.nan_to_num(rm)
                rms.append(rm)
                idx.append(i)
    rms = np.asarray(rms)
    return idx, rms

def computeTestStatistic_Diffs(groupX, groupY):
    """
    Takes two arrays. Each of which is a stack
    of single trial {RM or avg? decide}. 
    
    Avgs them to a summary trace and returns their bin-wise diff
    """    
    
    maskX= np.ma.masked_invalid(groupX)
    avgX = maskX.mean(axis=0) # ignores inf and nan
    maskY= np.ma.masked_invalid(groupY)
    avgY = maskY.mean(axis=0) # ignores inf and nan
    return avgX-avgY

def getLabels(unit, alley):
    """
    Get actual trial labels for a group
    Group defined as visits to a given txt at given alley
    """
    visits = unit.alleys[alley]
    labels = []
    for visit in visits:
        labels.append(visit['metadata']['stimulus'])
    return labels

def genSingleNullStat(unit, alley, txtX, txtY, labels):
    """
    Generate a single null test statistic (diff x-y here)
    Shuffle labels, recompute means and take diff. 1x
    """
    shuffLabels = np.random.permutation(labels)
    idxX, rmsX = poolTrials(unit, alley, shuffLabels, txtX)
    idxY, rmsY = poolTrials(unit, alley, shuffLabels, txtY)
    null = computeTestStatistic_Diffs(rmsX, rmsY)
    return null

def genRealStat(unit, alley, txtX, txtY):
    labels = getLabels(unit, alley)
    idxX, rmsX = poolTrials(unit, alley, labels, txtX)
    idxY, rmsY = poolTrials(unit, alley, labels, txtY)
    stat = computeTestStatistic_Diffs(rmsX, rmsY)
    return stat

def computeBandThresh(nulls, alpha, side):
    '''Given a list of null array traces, find ordinate at 
    at each point that  admits a proportion of nulls equal to cutoff'''
    
    if side == 'upper':
        isReversed = True
    elif side == 'lower':
        isReversed = False
        
    propNull = int(((alpha / 2) * len(nulls)) + 1)
    datarange = range(len(nulls[0]))
    significanceBand = []
    for point in datarange:
        nullOrdinates = nulls[:,point]
        sortedVals = list(sorted(nullOrdinates, reverse=isReversed))
        significanceBand.append(sortedVals[propNull - 1]) #explicitly +1 to cutoff and -1 here to keep clear where thresh is and how 0idx works

    significanceBand = np.asarray(significanceBand)
    return significanceBand

def computeGlobalCrossings(nulls, lowerBand, upperBand):
    """
    Given an array of null test statistics, compute 
    the number of crossings *anywhere* given the supplied
    significance bands. Return proportion (obs. p-value)
    """
    
    passBools = [any(np.logical_or(probe > upperBand, probe < lowerBand)) for probe in nulls] # eg [T,F,F,T..etc]
    return sum(passBools)/len(passBools)

def global_FWER_alpha(nulls, unit, alpha=0.05): # fwerModifier should be 3 (txts) x n alleys. 9 in beltway task. But below adjust by how many alleys actually have activity so 9 may become smaller
    """
    Calculates the global, FWER corrected p-value at each bin of the data trace
    Returns the actual global P and the bands of test statistic ordinates that
    are the thresholds. 
    
    """
    
    FWERalphaSelected = None
    globalLower, globalUpper = None, None
    
    validalleys = []
    for a in [16, 17, 3, 1, 5, 7, 8, 10, 11]:
        valid = Filt.checkMinimumPassesActivity(unit, a)
        validalleys.append(valid)
    multCompFactor = sum(validalleys)
    
    if multCompFactor == 0:
        pass

    else:
        fwerModifier = 3*multCompFactor
        FWERalpha = (alpha / fwerModifier)  # nb this is a proportion (decimal) not a list cutoff (integer)
        alphaIncrements = np.linspace(0.017, 1e-6, 50) # start at 0.017 because thats the largest the adj p value could be: 0.05/(3*1)
        fwerSatisfied = False
        for adjustedAlpha in alphaIncrements:
            if not fwerSatisfied:
                lowerBand, upperBand = computeBandThresh(nulls, adjustedAlpha, 'lower'), computeBandThresh(nulls, adjustedAlpha, 'upper')
                propCrossings = computeGlobalCrossings(nulls, lowerBand, upperBand)
                if propCrossings < FWERalpha: 
                    fwerSatisfied = True
                    FWERalphaSelected = adjustedAlpha
                    globalLower, globalUpper = lowerBand, upperBand
                


    return FWERalphaSelected, globalLower, globalUpper

def genNNulls(n, unit, alley, txtX, txtY):
    """
    Generates n null test statistics, hard coded
    now to be the binwise diff of avg(txtA) - avg(txtB)
    Returns np array nXl where l is length of 1d RM in bins
    """
    nulls = np.empty((0,Def.singleAlleyBins[0]-1)) # by convention long dim is first
    labels = getLabels(unit, alley)
    for i in range(n):
        null = genSingleNullStat(unit, alley, txtX, txtY, labels)

        nulls = np.vstack((nulls, null))
    return nulls

def unitPermutationTest_SinglePair(unit, alley, txtX, txtY, nnulls, plot=False, returnInfo=True):
    """
    Wrapper function for global_FWER_alpha() that plots results
    """
    
    nulls = genNNulls(nnulls,unit,alley,txtX,txtY)
    FWERalphaSelected, glowerBand, gupperBand = global_FWER_alpha(nulls, unit)
    if FWERalphaSelected == None:
        glowerBand, gupperBand, pwAlphaLower, pwAlphaUpper = None, None, None, None
        globalCrossings, pointwiseCrossings, bounds, stat = None, None, None, None
    else:
        stat = genRealStat(unit, alley, txtX, txtY)
    
        #Below, calculate the pw alpha bc significantly modulated regions are defined
        # as those that pass the global band somewhere but then their extent is defined
        # as the whole region where they pass the pointwise band. See Buzsaki paper. 
        pwAlphaUpper, pwAlphaLower = computeBandThresh(nulls, 0.05, 'upper'), computeBandThresh(nulls, 0.05, 'lower')
        globalCrossings = np.where(np.logical_or(stat > gupperBand, stat < glowerBand))[0]
        
        if globalCrossings.shape[0] > 0:
            pointwiseCrossings = np.where(np.logical_or(stat > pwAlphaUpper, stat < pwAlphaLower))[0]
        else:
            globalCrossings, pointwiseCrossings = None, None
                
        if plot:
            plt.plot(nulls.T, 'k', alpha=0.4)
            plt.plot(stat,'g')
            plt.xlabel("Linearized Position, Long Axis of Alley")
            plt.ylabel("Difference in Firing Rate")
            plt.title(f"Permutation Test Results for Texture {txtX} vs {txtY} on Alley {alley}")
            for band, style in zip([glowerBand, gupperBand, pwAlphaLower, pwAlphaUpper], ['r', 'r', 'r--', 'r--']):
                plt.plot(band, style)
                
    if returnInfo:
        bounds = glowerBand, gupperBand, pwAlphaLower, pwAlphaUpper
        return globalCrossings, pointwiseCrossings, bounds, stat
    
    
def permutationResultsLogger(d,fname):
    
    doesPass = False
    for alley in [1,3,5,7,8,10,11,16,17]:
        for pair in ["AB", "BC", "CA"]:
            for crossType in ["global", "pointwise"]:
                if d[alley][pair][crossType] != 'XXX':
                    doesPass = True
    
    if doesPass:
        savename = fname + "_PASS"
    else:
        savename = fname
    
    with open(savename+'.csv', "w") as f:
        w = csv.writer(f, delimiter = ' ')
        for alley in [1,3,5,7,8,10,11,16,17]:
            w.writerow([alley])
            for pair in ["AB", "BC", "CA"]:
                w.writerow([pair])
                for crossType in ["global", "pointwise"]:
                    w.writerow([crossType, d[alley][pair][crossType]])
        f.close()
        
def unitPermutationTest_AllPairsAllAlleys(unit, nnulls,fpath, logger=True, plot='sepFile'):
    """
    Wrapper function to complete permutation tests for a unit
    across all alleys and all pairwise stim (A,B,C) combinations
    
    Pointwise p-value is set to 0.05
    Global p-value is set to 0.00098 (0.05/(3*17))
    
    All perm tests can be saved to a file for later use, depending on option:
    Plots will be in a 17x3 grid where each row is an alley 1-17
    and each column is a test stat in order AB, BC, CA
        
    plot = False -> Do not plot
    plot = sepFile -> Plot all test results to it's own file in the fpath dir
    plot = addFile -> Do not save as this plot is an addon to another file's
           plotting routines (which will save the file itself)
    """
    
    if plot != False:
        fig, axes = plt.subplots(9, 3, figsize=(12,12), dpi=200) #bigger plot, bigger dpi
    
    pairs = ["AB", "BC", "CA"]
    fname = fpath + f"{stamp}_{unit.name}_{Def.singleAlleyBins[0]-1}bins_{Def.smoothing_1d_sigma}smooth_permutationResults"
    crossings = {i:{pair:{'global':"XXX", 'pointwise':"XXX"} for pair in pairs} for i in [1,3,5,7,8,10,11,16,17]}
    
    axCounter = 0 
    for alley in [1,3,5,7,8,10,11,16,17]:
        
        valid = Filt.checkMinimumPassesActivity(unit, alley)
        if valid:
        
            print(f"Running Permutation test in alley {alley}")
        
            for pair in pairs:
                
                txtX, txtY = pair[0], pair[1]
                globalCrossings, pointwiseCrossings, bounds, stat = unitPermutationTest_SinglePair(unit, alley, txtX, txtY, nnulls, 
                                                                                           plot=False, returnInfo=True)
                if globalCrossings is not None:
                    crossings[alley][pair]['global'] = globalCrossings
                    crossings[alley][pair]['pointwise'] = pointwiseCrossings
                    
                conditionName = unit.name + "_" + str(alley) + "_" + pair
                    
                if plot != False and bounds[0] is not None:
                    # the plot keyword will tell plotting fx whether to save sep or leave live for sep file to save
                    plotPermutationResults(unit, bounds, stat, conditionName, globalCrossings, pointwiseCrossings, fig.axes[axCounter])
                
                axCounter += 1 # increment to get the next subplot next iteration.
                
        else:
            print(f"Insufficient activity in alley {alley}")
                
    plt.suptitle(f"Permutation Test Results for {unit.name}")
                         
                    
    if logger == True:
        permutationResultsLogger(crossings, fname)
        
    if plot == 'sepFile':      
        # just in case this is buggy in future: when sep script is saving the fpath var is ''
        plt.savefig(fname + ".svg")
        plt.close()
        
    elif plot == 'addFile':
        pass # just to be explicit that if another script is saving this plot
             # to its own set of plots (e.g the ratemap routine) then leave open
    
    
def plotPermutationResults(unit, bounds, stat, conditionName, globalCrossings, pointwiseCrossings, ax):
    """
    If the observed test statistic passes the test, plot bounds.
    Plot test statistic and original linear ratemaps on top
    Does not save, that is done in the wrapper fx for all pairs/alleys (or 
    in sep script calling it)
    """
    colorLookup = {'A':'r', 'B':'b', 'C': 'g'} # keep color coordination
    
    # Get the real traces. Should refactor so I don't need to do this here and in test itself. 
    txtX, txtY = conditionName.split("_")[2]
    alley = int(conditionName.split("_")[1])
    labels = getLabels(unit, alley)
    _, rmsX = poolTrials(unit, alley, labels, txtX)
    _, rmsY = poolTrials(unit, alley, labels, txtY)
    traceX, traceY = np.mean(rmsX, axis=0), np.mean(rmsY, axis=0)
    
    g_upper, g_lower, pw_upper, pw_lower = bounds
    
    ax.fill_between(range(len(g_upper)), g_upper, g_lower, color='cornflowerblue')
    ax.fill_between(range(len(pw_upper)), pw_upper, pw_lower, color='darkblue')
    
    ax.plot(stat, 'k')
    ax.plot(traceX, colorLookup[txtX])
    ax.plot(traceY, colorLookup[txtY])
    
    # Were plotting all test results so if it failed, no crossings to highlight
    if globalCrossings is not None:
        ax.scatter(globalCrossings, stat[globalCrossings], c='cornflowerblue', marker='^')
        ax.scatter(pointwiseCrossings, stat[pointwiseCrossings], c='darkblue', marker='^')
    
    if globalCrossings is not None:
        ax.set_title(f"{conditionName.split('_')[1:]}", color='r')
    else:
        ax.set_title(f"{conditionName.split('_')[1:]}", color='k')

    
       
if __name__ == '__main__':
    
    rat = "R859"
    expCode = "BRD3"
    datafile = f"E:\\Ratterdam\\{rat}\\{rat}{expCode}\\"
    fpath = f"E:\\Ratterdam\{rat}\\permutation_tests\\{expCode}\\"
    stamp = util.genTimestamp()
    
    alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
    
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
        
    print(expCode)
    
    for subdir, dirs, fs in os.walk(datafile):
        for f in fs:
            if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f:
                clustname = subdir[subdir.index("TT"):] + "\\" + f
                unit = Core.UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
                unit.loadData_raw()
                rm = util.makeRM(unit.spikes, unit.position)            
                if np.nanpercentile(rm,Def.wholetrack_imshow_pct_cutoff) >= 0.75:
                    print(clustname)
                    unitPermutationTest_AllPairsAllAlleys(unit, 5000, fpath)






