# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 20:57:51 2021

@author: Ruo-Yah Lai
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal
import csv
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.path as path
import scipy
from math import ceil, sqrt
from confounds import directionFilterS
import ratterdam_Defaults as Def
import utility_fx as util
from ratterdam_RepetitionCoreFx import loadRepeatingUnit
from confounds import graphDirRatemaps


def graphVisits(unit, subfield, title=""):
    """
    Graphs all visits to a field
    """
    cm = np.array([1, Def.ptsCm, Def.ptsCm])
    pos = unit.position/cm[None,:]
    perim = unit.perimeters[subfield]/Def.ptsCm
    
    fig, ax = plt.subplots()
    ax.plot(perim[:,0], perim[:,1])
    for i in range(len(unit.visits[subfield])):
        start = np.where(pos[:,0] == unit.visits[subfield][i][0])[0]
        end = np.where(pos[:,0] == unit.visits[subfield][i][-1])[0]
        ax.plot(pos[int(start):int(end+1),1], pos[int(start):int(end+1),2])
    ax.set_title(title)
    ax.axis('equal')


def graphVisitsDir(unit, subfield, title=""): 
    """
    Graph all visits to a field color coded by direction
    """
    fig, ax = plt.subplots()
    ax.plot(unit.perimeters[subfield][:,0], unit.perimeters[subfield][:,1])
    for visit in range(len(unit.visits[subfield])):
        visitPos = np.logical_and(unit.position[:,0] >= unit.visits[subfield][visit][0],
                                  unit.position[:,0] <= unit.visits[subfield][visit][-1])
        visitSpk = np.logical_and(unit.spikes[:,0] >= unit.visits[subfield][visit][0],
                                  unit.spikes[:,0] <= unit.visits[subfield][visit][-1])
        posDir, spikesDir = directionFilterS(unit.position[visitPos], unit.spikes[visitSpk])
        colors = ["black","red","green","blue"]
        for i in range(4):
            ax.scatter(posDir[i][:,1], posDir[i][:,2], c=colors[i], s=1)
    ax.set_title(title)
    ax.axis('equal')
            

def rateDir(unit, title, test):
    """
    Plots a bar graph of the average rate vs NESW for each subfield. Calculates
    stats using Kruskal Wallis and Tukey
    
    test: "kruskal" or "anova"
    """
    #calculates rates for each direction using visits as samples
    avgRate = np.empty((len(unit.perimeters), 4)) #rows=subfields, columns=NESW
    visitRates = [] #levels of lists: subfields, directions, visits
    for i in range(len(unit.perimeters)):
        visitRates.append([[], [], [], []])
        contour = path.Path(unit.perimeters[i])
        PinC = unit.position[contour.contains_points(unit.position[:,1:])]
        SinC = unit.spikes[contour.contains_points(unit.spikes[:,1:])]
        for visit in range(len(unit.visits[i])):
            visitPos = np.logical_and(PinC[:,0] >= unit.visits[i][visit][0],
                                      PinC[:,0] <= unit.visits[i][visit][-1])
            visitSpk = np.logical_and(SinC[:,0] >= unit.visits[i][visit][0],
                                      SinC[:,0] <= unit.visits[i][visit][-1])
            
            posDir, spikesDir = directionFilterS(PinC[visitPos], SinC[visitSpk])
            for j in range(4):
                if len(posDir[j]) == 0:
                    visitRate = np.nan
                else:
                    visitRate = len(spikesDir[j]) / len(posDir[j]) * 30
                visitRates[i][j].append(visitRate)
        
        for j in range(4):
            avgRate[i,j] = np.nanmean(visitRates[i][j])
    
    #dimensions for subplots
    rows = round(sqrt(len(unit.perimeters)))
    subfieldDim = [rows, ceil(len(unit.perimeters)/rows)]
    
    #graphing and calculating stats
    fig, axs = plt.subplots(subfieldDim[0], subfieldDim[1], figsize=(subfieldDim[1]*3,subfieldDim[0]*3))
    if len(unit.perimeters) > 1:
        axs = axs.flatten()
    else: #if there is only 1 field so only one axis for graphing
        axs = [axs] #make axs subscriptable
    tukeys = []
    for i in range(len(visitRates)): #i = subfield number
        sig = ""
        if test == "kruskal":
            Hstats, p = kruskal(visitRates[i][0], visitRates[i][1], visitRates[i][2], visitRates[i][3], nan_policy="omit")
        elif test == "anova":
            noNan = []
            for j in range(4):
                noNan.append(np.array(visitRates[i][j])[~np.isnan(visitRates[i][j])])
            f, p = scipy.stats.f_oneway(noNan[0], noNan[1], noNan[2], noNan[3])
            
        if p < 0.05:
            sig = "*"
            groups = []
            allVisitRates = np.empty(0)
            directions = [["N"], ["E"], ["S"], ["W"]]
            for j in range(4):
                rates = np.array(visitRates[i][j])[~np.isnan(visitRates[i][j])]
                allVisitRates = np.hstack((allVisitRates, rates))
                groups = groups + directions[j]*len(rates)
            tukey = pairwise_tukeyhsd(allVisitRates, groups)
            tukeys.append(tukey)
        else:
            tukeys.append(False)
        axs[i].bar(np.arange(4), avgRate[i])
        axs[i].set_title(f"Subfield {i} {sig}")
        axs[i].set_ylabel("Average firing rate (Hz)")
        axs[i].set_xticks(np.arange(4))
        axs[i].set_xticklabels(["N", "E", "S", "W"])
    for i in range(len(unit.perimeters), subfieldDim[0]*subfieldDim[1]):
        axs[i].axis("off")
    fig.suptitle(title, y=1.08)
    fig.tight_layout()
    
    return fig, tukeys


def readCells(file, df):
    """
    Returns tetrode\\cell for all cells according to the tabulations file
    """
    with open(df+file,"r") as f:
        lines = f.readlines()
        tabulations = [line.split(",") for line in lines]
    tabulations = np.array(tabulations)
    cells = tabulations[1:, 0]
    return cells


def rateDirs(file, title, timestamp, dfData, test, df2="W:/Ratterdam/R859/ratterdam_tabulations/", dfGenerated="R:/072720-/Directionality/"):
    """
    file: file name of the file with repetition tabulations
    title: part of the title of file being generated (include rat and day)
    dfData: path to data
    test: "kruskal" or "anova"
    df2: path to repetition tabulations
    dfGenerated: path to where to save the generated files
    
    Make rate vs direction bar graphs, directional ratemaps, and csv's of 
    directional tuning of fields
    """
    vthresh = Def.velocity_filter_thresh
    cells = readCells(file, df2)
    for cell in cells:
        cellTitle = cell.replace("\\cl-maze", " ")
        unit = loadRepeatingUnit(dfData, cell)
        fig, tukeys = rateDir(unit, f"{title} {cellTitle}\nvthresh = {vthresh} cm/s  {test}", test)
        
        fig.savefig(dfGenerated + f"{timestamp} - {title} {cellTitle} rate vs direction.png",
                    bbox_inches="tight")
        plt.close()
        
        #graphDirRatemaps(unit, f"{title} {cellTitle}\nvthresh = {vthresh} cm/s")
        #plt.savefig(dfGenerated + f"{timestamp} - {title} {cellTitle} directional ratemap.png",
        #            bbox_inches="tight")
        #plt.close()
        
        with open(f"{dfGenerated}{timestamp} - {title} {cellTitle} directionality.csv", "w", newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Field", test, "EN","ES","EW","NS","NW","SW"])
            for i in range(len(unit.perimeters)):
                if tukeys[i]:
                    csvwriter.writerow([i, True, tukeys[i].reject[0], tukeys[i].reject[1], tukeys[i].reject[2]
                                        , tukeys[i].reject[3], tukeys[i].reject[4], tukeys[i].reject[5]])
                else:
                    csvwriter.writerow([i, False])
                