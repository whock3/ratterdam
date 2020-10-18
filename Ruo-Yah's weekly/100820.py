# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:49:45 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as path
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import csv
import ast
from bisect import bisect_left
from copy import deepcopy

import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import utility_fx as util
from placeFieldBorders import reorderBorder
from alleyTransitions import possibleTurnsA, possibleTurnsI


def graphRM(position, spikes, title, percentile=95, perim=[], mycm="jet"):
    n = util.makeRM(spikes, position)
    vmax = np.nanpercentile(n, percentile)
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title + "\n" + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    im = ax.imshow(n, cmap=mycm, origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    colors = ["b","g","r","k","c","m","y","orange"]
    for i in range(len(perim)):
        ax.plot(perim[i][:,0]/4.72, perim[i][:,1]/4.72, color=colors[i], 
                label=f"Subfield {i}")
    ax.legend()
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")

def filterField(unit, index, rateThresh=0.2, pctThresh=10):
    """
    Inputs: unit - Unit class object
            index - index of field to filter
            rateThresh - pixel rate that will be used to see % pixels below it
            pctThresh -  % of pixels (0-100) with rate below rateThresh above which field should be discarded
    Returns: True/False, as to whether field passes. True = pass
    Function to filter detected place fields by what percent of the field is 'too quiet'. Count pct
    of pixels with rate below rateThresh. if this pct exceeds pctThresh, return False.
    """
    perim = unit.perimeters[index]
    contour = path.Path(perim)
    spkIn = unit.spikes[contour.contains_points(unit.spikes[:,1:])]
    occIn = unit.position[contour.contains_points(unit.position[:,1:])]
    rm=util.makeRM(spkIn,occIn)
    area = np.sum(~np.isnan(rm.flatten()))
    binsBelowThresh = np.where(rm.flatten()<=rateThresh)[0].shape[0]
    pct = (binsBelowThresh/area)*100
    if pct > pctThresh:
        return False
    else:
        return True
    
    
def filterFields(unit):
    """
    Test whether fields are good and give tight fitting borders to good fields
    """
    perims = []
    for i in range(len(unit.perimeters)):
        if filterField(unit, i):
            perims.append(reorderBorder(unit.repUnit.PF[i].perimeter, i))
    unit.perimeters = perims
    return unit


#def loadSpikes(df, clustName, position):
#    with open(df+"sessionEpochInfo.txt","r") as f:
#        lines = f.readlines()
#    start, end = int(lines[0].split(',')[0]), int(lines[0].split(',')[1])
#    clust = np.asarray(util.read_clust(df+clustName))
#    clust = clust[(clust >= start) & (clust <= end)]
#    spikexy = util.getPosFromTs(clust,position)
#    return np.column_stack((clust,spikexy))

            
#modified from turnsInFieldIO in alleyTransitions
def turnsInFieldIO(turns, unit, position, filename, df, subfield, suptitle="", percentile=98):
    """
    Graphs the turns associated with a subfield, separated by into/inside/out of
    turns: from alleyTransitions
    position: velocity filtered
    filename: the file with which locations a field is in
    """
    with open(df+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    fieldLoc = ast.literal_eval(data[subfield+1][1])
    
    turns2 = np.empty((0,6)) #allo before turn, ego, allo after turn, ts of exit, ts of entry, into/inside/out of/through
    rates = np.empty(0)
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            #if neither the alley exited nor the alley entered is in the field, it's a through turn
            if len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])) and len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],3))))
            #elif alley exited isn't in the field, it's an into turn
            elif len(fieldLoc + [turn[5]]) == len(set(fieldLoc + [turn[5]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],0))))
            #elif alley entered isn't in the field, it's an out of turn
            elif len(fieldLoc + [turn[7]]) == len(set(fieldLoc + [turn[7]])):
                turns2 = np.vstack((turns2, np.hstack((turn[:5],2))))
            #else (alley exited and alley entered are both in the field), it's an inside turn
            else:
                turns2 = np.vstack((turns2, np.hstack((turn[:5],1))))
            start = bisect_left(unit.spikes[:,0], float(turn[3]))
            end = bisect_left(unit.spikes[:,0], float(turn[4]))
            rate = (end-start) / (float(turn[4])-float(turn[3])) * 1e6
            rates = np.hstack((rates,rate))
            
            
    ebins = np.arange(1,5)
    abins = np.arange(1,5)
    meanRates = np.zeros((4,4,4)) #4 x number of ebins x number of abins
    Ns = np.empty((4,4,4))
    for ebin in ebins:
        for abin in abins:
            #into
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 0)
            meanRates[0,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[0,ebin-1,abin-1] = np.sum(inBin)
            
            #inside
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 1)
            meanRates[1,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[1,ebin-1,abin-1] = np.sum(inBin)
            
            #out of
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 2)
            meanRates[2,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[2,ebin-1,abin-1] = np.sum(inBin)
            
            #through
            inBin = (turns2[:,1].astype(float) == ebin) & (turns2[:,2].astype(float) == abin) \
                    & (turns2[:,5].astype(float) == 3)
            meanRates[3,ebin-1,abin-1] = np.nanmean(rates[inBin])
            Ns[3,ebin-1,abin-1] = np.sum(inBin)
    
    possTurns = np.full((4,4,4), False)
    for loc in fieldLoc:
        if loc.isdigit():
            possTurns = np.logical_or(possTurns, possibleTurnsA[:,int(loc),:,:])
        else:
            possTurns = np.logical_or(possTurns, possibleTurnsI[:,ord(loc)-65,:,:])
        
    fig = plt.figure(figsize=(12,8))
    gs = GridSpec(2, 3, figure=fig)
    vmax = np.nanpercentile(meanRates,95)
    titles = ["Into the field", "Inside the field", "Out of the field", "Through the field"]
    fig.suptitle(f"{suptitle} in {fieldLoc}"+"\nVelocity filtered with threshold = 3 cm/s", y=1.05)
    norm = Normalize(vmax=vmax)
    m = cm.ScalarMappable(norm=norm, cmap="jet")
    for i in range(4):
        ax = fig.add_subplot(gs[i//3, i%3])
        ax.set_title(titles[i])
        ax.set_xticks(np.arange(1.5, 5.5, 1))
        ax.set_xticklabels(["N", "E", "S", "W"])
        ax.set_yticks(np.arange(1.5, 5.5, 1))
        ax.set_yticklabels(["S", "R", "B", "L"])
        ax.set_xlim(1, 5)
        ax.set_ylim(1, 5)
        ax.set_xlabel("Allocentric turn direction")
        ax.set_ylabel("Egocentric turn direction")
    
        meanRates2 = deepcopy(meanRates[i])
        meanRates2[np.where(meanRates[i] != meanRates[i])] = 0
        colors = m.to_rgba(meanRates2)
        #alphas = Ns/np.percentile(Ns,90) *0.8 + 0.2
        #alphas[np.where(alphas > 1)] = 1
        #colors[:,:,3] = alphas
        colors[np.where(meanRates[i] != meanRates[i])] = 0
        for j in range(colors.shape[0]):
            for k in range(colors.shape[1]):
                patch = patches.Rectangle((abins[k],ebins[j]), 1, 1, color=colors[j,k,:])
                ax.annotate(int(Ns[i,j,k]), (k+1.5,j+1.5), color=(0.5,0.5,0.5))
                if not possTurns[i,j,k] and Ns[i,j,k]==0: #if the turn in "not possible" and the rat didn't make the turn
                    patch = patches.Rectangle((abins[k],ebins[j]), 1, 1, color=(0.5,0.5,0.5))
                ax.add_patch(patch)
    axcb = fig.colorbar(m, ax=ax)
    axcb.set_label("Rate (Hz)")
    
    n = util.makeRM(unit.spikes, position)
    vmax = np.nanpercentile(n, percentile)
    ax = fig.add_subplot(gs[1, 1:])
    ax.set_title(f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    colors = ["b","g","r","k","c","m","y","orange","b","g","r","k","c","m","y","orange"]
    for i in range(len(unit.perimeters)):
        ax.plot(unit.perimeters[i][:,0]/4.72, unit.perimeters[i][:,1]/4.72, color=colors[i], 
                label=f"Subfield {i}")
    ax.legend()
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Rate (Hz)")
    
    fig.tight_layout()


def bulkGraphs(turns, pos, timestamp, units, filenames, df, ratDayTetrode):
    """
    Makes multiple graphs using turnsInFieldIO
    timestamp: current timestamp
    units and filenames: lists, unit.perimeters filtered with filterField
    """
    for i in range(len(units)):
        for j in range(len(units[i].perimeters)):
            print(i,j)
            turnsInFieldIO(turns, units[i], pos, filenames[i], df, j, f"{ratDayTetrode} 1.{i+1} subfield {j}")
            plt.savefig("C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/Graphs/"
                        + timestamp + " - " + f"1.{i+1} subfield {j}" + ".png",
                        bbox_inches="tight")
            plt.close()
            
#bulkGraphs(turns,posFilt,"20201008-213157",[unit1,unit2,unit3,unit4,unit5,unit6,unit7,unit8,unit9],["100820 - unit1 locations","100820 - unit2 locations","100820 - unit3 locations","100820 - unit4 locations","100820 - unit5 locations","100820 - unit6 locations","100820 - unit7 locations","100820 - unit8 locations","100820 - unit9 locations"],"C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/","R859 D3 T6")