# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:08:20 2020

@author: whockei1
"""

#### Imports 
import numpy as np
import utility_fx as util
import os
from matplotlib import pyplot as plt
import scipy
import bisect
import csv, sys
import itertools

import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt



def createNeuralSymbols(symbolType):
    '''
    Hardcoded based on neural fr alphabet you want. All possible bins or bitstrings
    depending on what kind of alphabet.
    '''
    if symbolType == 'max' or symbolType == 'mean':
        mat  = np.empty((0, Def.singleAlleyBins[0]-1))
        for i,visit in enumerate(unit.alleys[alley]):
            mat = np.vstack((mat, visit['ratemap1d']))
        flatmat = np.ndarray.flatten(mat)
    
        
        freedman_diaconis_binWidth = 2*(scipy.stats.iqr(flatmat, nan_policy='omit')/np.cbrt(flatmat.shape[0]))
        numbins = round(np.nanmax(flatmat)/freedman_diaconis_binWidth)
        try:
            symbolBins = np.linspace(0, np.nanmax(flatmat), int(numbins)) 
        except:
            symbolBins = np.linspace(0, np.nanmax(flatmat), 40)
        return symbolBins
    
    elif symbolType == '3binMeanComp':
        symbolBins = list(itertools.product([0,1],repeat=3))
        return symbolBins

def convertActivitytoSymbols(unit, alley, symbolType, symbolBins, extent):
    
    if symbolType == 'mean':
        vec  = np.empty((0,1))
        for i,visit in enumerate(unit.alleys[alley]):
            c = bisect.bisect_left(symbolBins, np.nanmean(visit['ratemap1d'][extent[0]:extent[1]]))
            vec = np.vstack((vec, c))
        return vec
    
    elif symbolType == 'max':
        vec  = np.empty((0,1))
        for i,visit in enumerate(unit.alleys[alley]):
            c = bisect.bisect_left(symbolBins, np.nanmean(visit['ratemap1d'][extent[0]:extent[1]]))
            vec = np.vstack((vec, c))
        return vec
    
    elif symbolType == '3binMeanComp':
        # Divide a 12bin ratemap into 3 sections. Compute overall mean in each.
        # Then for each visit find mean in each secion and say if its above (1) or below/equal (0) mean
        rmstack = collectAllRMs(alley)
        mat = np.empty((0,3))
        for visit in unit.alleys[alley]:
            data = visit['ratemap1d']
            means = np.nanmean(np.nanmean(rmstack[:,:4],axis=0)), np.nanmean(np.nanmean(rmstack[:,4:8],axis=0)), np.nanmean(np.nanmean(rmstack[:,8:],axis=0))
            bindata = np.nanmean(data[:4]), np.nanmean(data[4:8]), np.nanmean(data[8:])
            symbol = np.greater(bindata, means)*1 # i think its arbitrary whether we say greater or less than
            mat = np.vstack((mat, symbol))
        return mat
    
    elif symbolType == '3binStdInterval':
        # Divide 12bin ratemap into 3 sections. Find mean+n*Std in each section for 3 std vals
        # For each visit, find mean in each section and see what interval it is in
        pass
        

def collectAllRMs(alley):
    """
    Collects into a matrix all linear ratemaps from all visits to the alley
    w/o regard to stimulus. Returns a NxB matrix of N trials by B bins in lin RM
    """
    stack = np.empty((0, Def.singleAlleyBins[0]-1))
    for trial in unit.alleys[alley]:
        stack = np.vstack((stack, trial['ratemap1d']))
    return stack


        

def calculateMutualInformation(neural_symbols_set, converted_array, shuffle):

    trialStimSymbols = np.asarray([visit['metadata']['stimulus'] for visit in unit.alleys[alley]])
    if shuffle == True:
        trialStimSymbols = np.random.permutation(trialStimSymbols)
    joint_table = np.zeros((len(neural_symbols_set),0))
    for i,stim in enumerate(['A','B','C']):
        p_txty = np.zeros((0,1))
        for j,sym in enumerate(neural_symbols_set):
            c = np.intersect1d(np.where((converted_array == sym).all(axis=1))[0],np.where(trialStimSymbols==stim)[0]).shape[0]
            p_txty = np.vstack((p_txty,c))
        joint_table  = np.hstack((joint_table, p_txty))

    mi = 0
    for i in range(3):
        for j in range(len(neural_symbols_set)):
            px = joint_table.sum(axis=0)[i]/joint_table.sum()
            py = joint_table.sum(axis=1)[j]/joint_table.sum()
            pxy = joint_table[j,i]/joint_table.sum()
            if pxy != 0:
                mi += pxy * np.log2(pxy/(px*py))
    return mi


### Subplots of all alleys for a cell - whole track
rat = sys.argv[1]
expCode = sys.argv[2]
datafile = f'E:\\Ratterdam\\{rat}\\{rat}{expCode}\\'
figurePath = f'E:\\Ratterdam\\{rat}\\mutual_information\\{expCode}\\'
alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
population = {}
allRealIs = {i:[] for i in Def.beltwayAlleys}
stamp = util.genTimestamp()
n=1000 # num shuffles
symbolType = '3binMeanComp' # Options: "Mean", "Max", "3binMeanComp"
for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            unit = Core.UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
            unit.loadData_raw()
            rm = util.makeRM(unit.spikes, unit.position)            
            if np.nanpercentile(rm,Def.wholetrack_imshow_pct_cutoff) >= 1.:
                print(clustname)

                passFlag = False # bool to toggle whether any alley has passed criteria
                passingalleys = [] # will populate with (alley, percentile of nulls) for any passing alleys
                
                # Check how many alleys have some appreciable activity and adjust your Bonf. correction by that
                validalleys = []
                for a in Def.beltwayAlleys:
                    valid = Filt.checkMinimumPassesActivity(unit, a)
                    validalleys.append(valid)
                multCompFactor = sum(validalleys) # we will adjust for multiple comparisons by considering the # alleys where there is 
                                                # appreciable activity (avg 1Hz overall) and therefore a valid test to even attempt
                if multCompFactor >= 1:
                    thresh = round((0.05/multCompFactor),2)*n


                    fig, ax = plt.subplots(3,3,figsize=(8,8))
                    for i, alley in enumerate(Def.beltwayAlleys):
                        axis = fig.axes[i]
                        if  validalleys[i] == True:
                            print(f"Alley {alley}")
                            stimulus_symbols = ['A','B','C']
                            neural_symbols_set = createNeuralSymbols(symbolType)
                            converted_array = convertActivitytoSymbols(unit, alley, symbolType, neural_symbols_set)
                            neural_symbols_set = range(len(neural_symbols_set)) 
                            realI = calculateMutualInformation(neural_symbols_set, converted_array, False)
                            allRealIs[alley].append(realI)

                            nulls = []
                            for i in range(n):
                                nullI = calculateMutualInformation(neural_symbols_set, converted_array, True)
                                nulls.append(nullI)
                            nulls = sorted(nulls)

                                
                            # check if the alley passes and save if so.
                            realPct = round(scipy.stats.percentileofscore(nulls,realI, kind='strict'),2) # what percentile of nulls is real mi score?
                            if realPct > (1-(0.05/multCompFactor))*100:
                                passFlag = True
                                passingalleys.append([alley, realPct])
                            
                            
                            axis.hist(nulls,50)
                            axis.vlines(realI,0,150,'r')
                            axis.vlines(nulls[int(-1*thresh)],0,150,'grey')
                            axis.set_title(f"Alley {Def.beltwayAlleyLookup[alley]}, {realPct}")
                        else:
                             axis.set_title(f"Insufficient Activity Alley {Def.beltwayAlleyLookup[alley]}")
                             print(f"{Def.beltwayAlleyLookup[alley]} not included in MI calc")

                    plt.suptitle(f"Whole Alley MI Unit {expCode} {unit.name}")
                    plt.savefig(figurePath+f"{stamp}_{unit.name}_MIPlot_{passFlag}.png", format='png',dpi=300)
                    plt.close()
                        
                        
                        
                    with open(figurePath+f"{stamp}_{unit.name}_MIData_{passFlag}.csv", "w") as file:
                        writer = csv.writer(file)

                        if len(passingalleys) > 0:
                            writer.writerows(passingalleys)
                        else:
                            writer.writerow("No alleys passed")
                                
                                