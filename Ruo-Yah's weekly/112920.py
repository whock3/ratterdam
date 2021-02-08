# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 15:57:29 2020

@author: Ruo-Yah Lai
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import ast
import sys
sys.path.insert(1, "C:/Users/Ruo-Yah Lai/Documents/GitHub/ratterdam")
import utility_fx as util
import ratterdam_Defaults as Def
from alleyTransitions import alleyTransitions


def makeRM2(spikes, position):
    """ makeRM adapted for avg RM for turns"""
    #camXMax, camXMin = [round(np.max(position[:,1])), round(np.min(position[:,1]))]
    #camYMax, camYMin = [round(np.max(position[:,2])), round(np.min(position[:,2]))]
    #rbins = int((camXMax-camXMin)/4.72)
    #cbins = int((camYMax-camYMin)/4.72)
    rows, cols = np.linspace(-20*Def.ptsCm, 20*Def.ptsCm, 21), np.linspace(-20*Def.ptsCm, 20*Def.ptsCm, 21)
    hs,xe,ye = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rows, cols])
    ho = np.histogram2d(position[:,2],position[:,1],bins=[rows, cols])[0]
    n = (hs*np.reciprocal(ho))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(ho==0)] = np.nan
    #cam = [floor(camXMin/4.72),ceil(camXMax/4.72),floor(camYMin/4.72),ceil(camYMax/4.72)]
    return n, ho


def graphRM2(position, pos2, spikes, suptitle, percentile=98, smoothing_2d_sigma=1):
    """ 
    graphRM adapted for avg RM for turns
    position, pos2, and spikes from shiftPos
    """
    rows, cols = np.linspace(-20*Def.ptsCm, 20*Def.ptsCm, 21), np.linspace(-20*Def.ptsCm, 20*Def.ptsCm, 21)
    n = []
    titles2 = []
    for i in range(8):
        if len(position[i]) == 0:
            n.append(np.full([20,20], np.nan))
        else:
            n1 = []
            ho = np.histogram2d(pos2[i][:,2], pos2[i][:,1], bins=[rows,cols])[0]
            for j in range(len(position[i])):
                n1.append(makeRM2(spikes[i][j], position[i][j])[0])
            n2 = np.nanmean(n1, axis=0)
            
            n2 = util.weird_smooth(n2,smoothing_2d_sigma)
            n2[np.where(ho==0)] = np.nan
            n.append(n2)
        titles2.append(f"n = {len(position[i])}")
    
    fig, axs = plt.subplots(2,4)
    vmax = np.nanpercentile(n[0], percentile)
    fig.suptitle(f"{suptitle}\nvthresh = 3 cm/s\nCutoff = {percentile}th percentile, {round(vmax,1)} Hz", y=1.12)
    titles = ["Straight", "Right", "Back", "Left", "North", "East", "South", "West"]

    ims = []
    axs = axs.flatten()
    for i in range(8):
        axs[i].set_title(titles[i] + "\n" + titles2[i])
        ims.append(axs[i].imshow(n[i], cmap="jet", origin="lower", vmin=0, vmax=vmax))#, extent=(-7,7,-7,7)))
        axs[i].axis("equal")
        axs[i].set_xticks([0,10,20])
        axs[i].set_xticklabels([-20,0,20])
        axs[i].set_yticks([0,10,20])
        axs[i].set_yticklabels([-20,0,20])
    cb = fig.colorbar(ims[0])
    cb.set_label("Rate (Hz)")
    fig.tight_layout()


def shiftPos(turns, unit, subfield, filename, df):
    """
    turns: from alleyTransitions
    filename: the file with which locations a field is in
    
    Returns: 
        shifted_pos: list of 8 lists each with arrays where each array is 1 turn
    pos shifted based on the 1st point of turns
        shifted_pos2: list of 8 arrays
        shifted_s: list of 8 lists of spikes
    """
    with open(df+filename, "r") as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    fieldLoc = ast.literal_eval(data[subfield+1][1])
    turns2 = np.empty((0,5)) #turns in the field
    for turn in turns:
        allLoc = fieldLoc + list(turn[5:8])
        if len(set(allLoc)) < len(allLoc):
            #exclude back around turns that are not in the field
            if len(set(turn[5:8])) == 2 and len(set(allLoc)) == len(allLoc)-1:
                continue
            turns2 = np.vstack((turns2, turn[:5].astype(float))) #allo before turn, ego, allo after turn, ts of exit, ts of entry
    
    shifted_pos = [[] for _ in range(8)]
    shifted_pos2 = [np.empty((0,3)) for _ in range(8)]
    shifted_s = [[] for _ in range(8)]
    for turn in turns2:
        start = np.where(unit.position[:,0] == turn[3])[0][0] #1st point of the turn (when alley is exited)
        end = np.where(unit.position[:,0] == turn[4])[0][0] #last point of the turn (when alley is entered)
        
        shiftedP = np.subtract(unit.position[start:end+1], np.array([0,unit.position[start,1],unit.position[start,2]]))
        srP = rotate(shiftedP, turn[0])
        
        if end < len(unit.position)-1:
            spike_idx = np.where(np.logical_and(unit.spikes[:,0]>unit.position[start,0], unit.spikes[:,0]<unit.position[end+1,0]))[0]
        else:
            spike_idx = np.where(np.logical_and(unit.spikes[:,0]>unit.position[start,0], unit.spikes[:,0]<unit.position[end,0]))[0]
        shiftedS = np.subtract(unit.spikes[spike_idx], np.array([0,unit.position[start,1],unit.position[start,2]]))
        if len(shiftedS) > 0:
            srS = rotate(shiftedS, turn[0])
        else:
            srS = np.empty((0,3))
        
        #shifted_pos[0].append(shiftedP)
        #shifted_pos2[0] = np.vstack((shifted_pos2[0],shiftedP))
        #shifted_s[0].append(shiftedS)
        #ego turns
        if turn[1] == 1: #S
            shifted_pos[0].append(srP)
            shifted_pos2[0] = np.vstack((shifted_pos2[0],srP))
            shifted_s[0].append(srS)
        elif turn[1] == 2: #R
            shifted_pos[1].append(srP)
            shifted_pos2[1] = np.vstack((shifted_pos2[1],srP))
            shifted_s[1].append(srS)
        elif turn[1] == 3: #B
            shifted_pos[2].append(srP)
            shifted_pos2[2] = np.vstack((shifted_pos2[2],srP))
            shifted_s[2].append(srS)
        elif turn[1] == 4: #L
            shifted_pos[3].append(srP)
            shifted_pos2[3] = np.vstack((shifted_pos2[3],srP))
            shifted_s[3].append(srS)
        
        #allo turns
        if turn[2] == 1: #N
            shifted_pos[4].append(shiftedP)
            shifted_pos2[4] = np.vstack((shifted_pos2[4],shiftedP))
            shifted_s[4].append(shiftedS)
        elif turn[2] == 2: #E
            shifted_pos[5].append(shiftedP)
            shifted_pos2[5] = np.vstack((shifted_pos2[5],shiftedP))
            shifted_s[5].append(shiftedS)
        elif turn[2] == 3: #S
            shifted_pos[6].append(shiftedP)
            shifted_pos2[6] = np.vstack((shifted_pos2[6],shiftedP))
            shifted_s[6].append(shiftedS)
        elif turn[2] == 4: #W
            shifted_pos[7].append(shiftedP)
            shifted_pos2[7] = np.vstack((shifted_pos2[7],shiftedP))
            shifted_s[7].append(shiftedS)
    return shifted_pos, shifted_pos2, shifted_s


def rotate(position, allo):
    """
    position: [ts,x,y]
    allo: allocentric direction before turning; 1=N, 2=E, 3=S, 4=W
    """
    theta = (allo-1)*np.pi/2
    x = position[:, 1]
    y = position[:, 2]
    rx = x*np.cos(theta) - y*np.sin(theta)
    ry = x*np.sin(theta) + y*np.cos(theta)
    return np.column_stack((position[:,0],rx,ry))


def trajectories(position, suptitle=""):
    """
    graph trajectories after shifting and rotating them
    """
    fig,axs = plt.subplots(2,4)
    titles = ["Straight", "Right", "Back", "Left", "North", "East", "South", "West"]
    fig.suptitle(suptitle, y=1.04)
        
    for i in range(2):
        for j in range(4):
            axs[i][j].set_title(titles[i*4+j] + "\n" + f"n = {len(position[i*4+j])}")
            for k in range(len(position[i*4+j])):
                axs[i][j].plot(position[i*4+j][k][:,1]/4.72, position[i*4+j][k][:,2]/4.72)
                axs[i][j].axis("equal")
    fig.tight_layout()
    

def bulkGraphs(unit, title, timestamp, rat, df, df2="R:\\072720-\\"):
    """
    Makes multiple graphs using graphRM2, which makes average rate maps for turns
    
    file: file name of the file with repetition tabulations
    title: part of the title of the generated files 
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    df: where to save the graphs
    df2: path to the files with which locations a field is in
    """
    _, turns = alleyTransitions(unit.position, rat)
    for i in range(len(unit.perimeters)):
        pos, pos2, spikes = shiftPos(turns, unit, i, f"20201127 - {title} locations", df2)
        graphRM2(pos, pos2, spikes, f"{title} subfield {i}")
        plt.savefig(df + f"{timestamp} - {title} subfield {i}.png",
                    bbox_inches="tight")
        plt.close()