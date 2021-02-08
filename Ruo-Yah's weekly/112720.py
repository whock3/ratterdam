# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 16:57:49 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
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
from alleyTransitions import possibleTurnsA, possibleTurnsI, alleyTransitions
import ratterdam_RepetitionCoreFx as RepCore


def turnsInFieldIO2(turns, unit, filename, df, subfieldDim, suptitle=""):
    """
    Graphs the turns associated with each subfield, separated by into/inside/out of
    Collapsed across egocentric directions
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    subfieldDim: list, dimensions (rows,columns) for putting all subfields on the same plot
    """
    with open(df+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    
    fig, axs = plt.subplots(subfieldDim[0], subfieldDim[1], figsize=(subfieldDim[1]*4,subfieldDim[0]*3.5))
    fig.suptitle(f"{suptitle}\nvthresh = 3 cm/s", y=1.07)
    axs = axs.flatten()
    for i in range(len(unit.perimeters)):
        fieldLoc = ast.literal_eval(data[i+1][1])
        
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
                
                
        abins = np.arange(1,5)
        meanRates = np.zeros((4,4)) #4 x number of abins
        Ns = np.empty((4,4))
        for abin in abins:
            #into
            inBin = (turns2[:,2].astype(float) == abin) & (turns2[:,5].astype(float) == 0)
            meanRates[0,abin-1] = np.nanmean(rates[inBin])
            Ns[0,abin-1] = np.sum(inBin)
                
            #inside
            inBin = (turns2[:,2].astype(float) == abin) & (turns2[:,5].astype(float) == 1)
            meanRates[1,abin-1] = np.nanmean(rates[inBin])
            Ns[1,abin-1] = np.sum(inBin)
                
            #out of
            inBin = (turns2[:,2].astype(float) == abin) & (turns2[:,5].astype(float) == 2)
            meanRates[2,abin-1] = np.nanmean(rates[inBin])
            Ns[2,abin-1] = np.sum(inBin)
                
            #through
            inBin = (turns2[:,2].astype(float) == abin) & (turns2[:,5].astype(float) == 3)
            meanRates[3,abin-1] = np.nanmean(rates[inBin])
            Ns[3,abin-1] = np.sum(inBin)
        
        possTurns = np.full((4,4), False)
        for loc in fieldLoc:
            if loc.isdigit():
                possTurns = np.logical_or.reduce((possTurns, possibleTurnsA[:,int(loc),0,:], possibleTurnsA[:,int(loc),1,:], possibleTurnsA[:,int(loc),2,:], possibleTurnsA[:,int(loc),3,:]))
            else:
                possTurns = np.logical_or.reduce((possTurns, possibleTurnsI[:,ord(loc)-65,0,:], possibleTurnsI[:,ord(loc)-65,1,:], possibleTurnsI[:,ord(loc)-65,2,:], possibleTurnsI[:,ord(loc)-65,3,:]))
            
        
        vmax = np.nanpercentile(meanRates,95)
        
        norm = Normalize(vmax=vmax)
        m = cm.ScalarMappable(norm=norm, cmap="jet")
        axs[i].set_title(f"Subfield {i} in {fieldLoc}\nCutoff = 95th %tile, {round(vmax,1)} Hz")
        axs[i].set_xticks(np.arange(1.5, 5.5, 1))
        axs[i].set_xticklabels(["N", "E", "S", "W"])
        axs[i].set_yticks(np.arange(1.5, 5.5, 1))
        axs[i].set_yticklabels(["Into", "Inside", "Out of", "Through"])
        axs[i].set_xlim(1, 5)
        axs[i].set_ylim(1, 5)
        axs[i].set_xlabel("Allocentric turn direction")
        
        meanRates2 = deepcopy(meanRates)
        meanRates2[np.where(meanRates != meanRates)] = 0
        colors = m.to_rgba(meanRates2)
        #alphas = Ns/np.percentile(Ns,90) *0.8 + 0.2
        #alphas[np.where(alphas > 1)] = 1
        #colors[:,:,3] = alphas
        colors[np.where(meanRates != meanRates)] = 0
        for j in range(colors.shape[0]):
            for k in range(colors.shape[1]):
                patch = patches.Rectangle((abins[k],j+1), 1, 1, color=colors[j,k,:])
                axs[i].annotate(int(Ns[j,k]), (k+1.5,j+1.5), color=(0.5,0.5,0.5))
                if not possTurns[j,k] and Ns[j,k]==0: #if the turn in "not possible" and the rat didn't make the turn
                    patch = patches.Rectangle((abins[k],j+1), 1, 1, color=(0.5,0.5,0.5))
                axs[i].add_patch(patch)
        axcb = fig.colorbar(m, ax=axs[i])
        axcb.set_label("Rate (Hz)")
        
        #n = util.makeRM(unit.spikes, unit.position)
        #vmax = np.nanpercentile(n, percentile)
        #ax = fig.add_subplot(gs[0, 1:])
        #ax.set_title(f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz")
        #im = ax.imshow(n, cmap="jet", origin="lower", vmin=0, vmax=vmax)
        #ax.set_xlabel("x coordinates (cm)")
        #ax.set_ylabel("y coordinates (cm)")
        #ax.set_xlim(0,160)
        #for i in range(len(unit.perimeters)):
        #    ax.plot(unit.perimeters[i][:,0]/4.72, unit.perimeters[i][:,1]/4.72, color=unit.colors[i], 
        #            label=f"Subfield {i}")
        #ax.legend(loc="lower right")
        #cb = fig.colorbar(im, ax=ax)
        #cb.set_label("Rate (Hz)")
    fig.tight_layout()


def bulkGraphs(file, title, timestamp, rat, df, df2, df3, df4):
    """
    Makes multiple graphs using turnsInFieldIO2
    
    file: file name of the file with repetition tabulations
    title: part of the title of the generated files (include rat and day)
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    df: path to data
    df2: path to repetition tabulations
    df3: path to the files with which locations a field is in
    df4: where to save the graphs
    """
    cells = util.readRepeatingCells(file, df2)
    unit = RepCore.loadRepeatingUnit(df, cells[0])
    _, turns = alleyTransitions(unit.position, rat)
    tooManyFields = []
    for cell in cells:
        cellTitle = cell.replace("\\cl-maze", " ")
        unit = RepCore.loadRepeatingUnit(df, cell)
        if len(unit.perimeters) <= 2:
            subfieldDim = [2,1]
        elif 2 < len(unit.perimeters) <= 4:
            subfieldDim = [2,2]
        elif 4 < len(unit.perimeters) <= 6:
            subfieldDim = [2,3]
        elif 6 < len(unit.perimeters) <= 9:
            subfieldDim = [3,3]
        elif 9 < len(unit.perimeters) <= 12:
            subfieldDim = [3,4]
        elif 12 < len(unit.perimeters) <= 16:
            subfieldDim = [4,4]
        else:
            tooManyFields.append([cell, len(unit.perimeters)])
            continue
        turnsInFieldIO2(turns, unit, f"20201127 - {title} {cellTitle} locations", df3, subfieldDim, f"{title} {cellTitle}")
        plt.savefig(df4 + f"{timestamp} - {title} {cellTitle}.png",
                    bbox_inches="tight")
        plt.close()