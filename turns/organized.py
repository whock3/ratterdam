# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 20:06:17 2020

@author: Ruo-Yah Lai
"""
import csv
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from RDP4_class import RDP4
from matplotlib.collections import LineCollection
from turn3 import turn3
from bisect import bisect
from matplotlib.colors import LinearSegmentedColormap
from RateMap import makeRM, weird_smooth
from math import ceil, floor
from copy import deepcopy


def velocity_filtering(position, threshold = 3, winsz=50):
    """
    from Will's ratterdam_DataFiltering
    """
    gradts, gradx, grady = np.gradient(position[:,0]), np.gradient(position[:,1]), np.gradient(position[:,2])
    gradx = [np.mean(gradx[0+i:winsz+i]) for i in range(len(gradx))]
    grady = [np.mean(grady[0+i:winsz+i]) for i in range(len(grady))]
#    gradx = np.asarray([i/Def.ptsCm for i in gradx])
#    grady = np.asarray([i/Def.ptsCm for i in grady])
    
    vx = np.asarray([1e6*(a/b) for a,b in zip(gradx,gradts)])
    vy = np.asarray([1e6*(a/b) for a,b in zip(grady,gradts)])
    v =  np.sqrt((vx**2)+(vy**2))  
    
    sv = [np.mean(v[0+i:winsz+i]) for i in range(len(v))]
    sv = np.asarray(sv)
    
    vf_pos = position[sv > threshold]
    
    return vf_pos


#originally from 041620
def read_pos(path="C:/Users/Ruo-Yah Lai/Desktop/My folder/College/Junior/K lab research/R859 OD3/", pos="pos.p.ascii", to_cm = True):
    data_folder = Path(path)
    file_to_open = data_folder / pos
    if to_cm == True:
        a = np.array([1, 4.72, 4.72])
        ptsCm = 1
    else:
        a = np.array([1,1,1])
        ptsCm = 4.72
    with open(file_to_open, 'r') as csv_file:
        data_iter = csv.reader(csv_file)
        data = [data for data in data_iter]
    data_np1 = np.array(data[50574:-17], dtype = "float64")
    data_np2 = data_np1[:,0:3]/a[None,:]
    data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
    return velocity_filtering(data_np3, 3*ptsCm)


def angle(dir):
    """
    from https://stackoverflow.com/questions/14631776/calculate-turning-points-pivot-points-in-trajectory-path
    """
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return np.arccos((dir1*dir2).sum(axis=1)/(
        np.sqrt((dir1**2).sum(axis=1)*(dir2**2).sum(axis=1))))


#from 051420
def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)


def compareEpsilons(position, b=28000):
    """
    from 042320: compare epsilons + primitive turn algorithm
    """
    epsilons = [0.5, 0.75, 1]
    min_angle = np.pi*0.25
    fig, ax = plt.subplots(1, 4, figsize=(24, 6))
    ax[0].plot(position[b:b+300, 1], position[b:b+300, 2])
    ax[0].set_title("Original")
    ax[0].scatter(position[b,1], position[b,2], marker = "+", color = "r", label = "first")
    ax[0].scatter(position[b+299,1], position[b+299,2], marker = "x", color = "r", label = "last")
    for i in range(len(epsilons)):
        RList = RDP4(position[b:b+300], epsilons[i]).ResultList
        directions = np.diff(RList[:,1:3], axis=0)
        theta = angle(directions)
        idx = np.where(theta>min_angle)[0]+1
        percent = round(100 - len(RList)/3, 2)
        ax[i+1].plot(RList[:,1], RList[:,2])
        ax[i+1].set_title(f"epsilon = {round(epsilons[i], 2)}" "\n" f"{percent}% pts thrown away")
        ax[i+1].scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
        ax[i+1].scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")
        ax[i+1].scatter(RList[idx,1], RList[idx,2], s = 16, color = "r", label='turning points')
    plt.legend()
    plt.show()
    

#from 050320: historical, comparing velocity filters
    

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    """
    from https://matplotlib.org/3.2.1/gallery/lines_bars_and_markers/scatter_hist.html
    """
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)
    plt.xlabel("instantaneous velocity")
    plt.ylabel("angle")

    # now determine nice limits by hand:
    binwidthx = 1
    binwidthy = np.pi/40
    xmax = np.max(np.abs(x))
    limx = (int(xmax/binwidthx) + 1) * binwidthx
    binsx = np.arange(0, limx + binwidthx, binwidthx)
    binsy = np.arange(0, np.pi, binwidthy)
    ax_histx.hist(x, bins=binsx)
    ax_histy.hist(y, bins=binsy, orientation='horizontal')


def anglesVelocities(position, epsilon, b1=48600, b2=48900):
    """
    from 050320: 
        2D histogram (not heatmap) of angles and instantaneous velocities
        allows different epsilons
    """
    RList = RDP4(position[b1:b2], epsilon).ResultList
    diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
    vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
    vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
    v =  np.sqrt((vx**2)+(vy**2))
    directions = np.diff(RList[:,1:3], axis=0)
    theta = angle(directions)  
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    
    # start with a square Figure
    fig = plt.figure(figsize=(8, 8))
    
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)
    
    scatter_hist(v[1:len(v)], theta, ax, ax_histx, ax_histy)
    plt.title(f"epsilon = {epsilon}")
    plt.show()
    

#from 050820: historical, 2D histogram of angles and instantaneous velocities, different velocity thresholds
    

def velSegments(RList, idx, title=""): 
    """
    originally from 050820: modified for single plot and external turn algorithm
    RList = after RDP; idx = indexes for turns; title = graph title
    """
    diffts, diffx, diffy = np.diff(RList[:,0]), np.diff(RList[:,1]), np.diff(RList[:,2]) 
    vx = np.asarray([1e6*(a/b) for a,b in zip(diffx,diffts)])
    vy = np.asarray([1e6*(a/b) for a,b in zip(diffy,diffts)])
    v =  np.sqrt((vx**2)+(vy**2))
    
    fig, ax = plt.subplots(1, 1) #1, 3, figsize = (19, 5), gridspec_kw={'width_ratios': [6, 6, 7]})
    ax.set_xlim(0, 120) #np.min(RList[:, 1])-3, np.max(RList[:, 1])+3)
    ax.set_ylim(0, 100) #np.min(RList[:, 2])-3, np.max(RList[:, 2])+3)
    ax.scatter(RList[0,1], RList[0,2], marker = "+", color = "r", label = "first")
    ax.scatter(RList[-1,1], RList[-1,2], marker = "x", color = "r", label = "last")
     
    ax.scatter(RList[idx,1], RList[idx,2], s = 16, color = "r", label='turning points')
    
    RList2 = np.reshape(RList[:, 1:3], (len(RList), 1, 2))
    segments = np.concatenate([RList2[0:-1], RList2[1:]], axis=1)
    lc = LineCollection(segments)
    lc.set_array(v)
    ax.add_collection(lc)
    axcb = fig.colorbar(lc)
    axcb.set_label('Instantaneous velocity (cm/s)')
    ax.set_title(title)
    ax.axis("equal")
    ax.legend()


#from 051420: historical, first version of turn algorithm (with summing)
#from 051420: historical, whole session theta > 1.5, velocity < 35
#from 052320: historical, another version of turn algorithm


def classifyTurns(b1, b2, position, min_angle1=np.pi/12, min_angle2=np.pi/4, turn_range_start=5*4.72, vCutoff=15*4.72, vIncrease=1/7, epsilon=4.72):
    """
    from 052320
    min_angle1 = min angle for the first angle
    min_angle2 = min angle for angles within a range to add up to
    turn_range = angles within x cm of each "turn" shouldn't add up to be a circle
    """
    RList = RDP4(position[b1:b2], epsilon).ResultList
    idx3, theta_sum2, idx2_1 = turn3(RList, min_angle1, min_angle2, turn_range_start, vCutoff, vIncrease)
    
    ego_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if theta_sum2[i] < (3/4*np.pi):
            ego_turns[i] = 0
        elif theta_sum2[i] > (5/4*np.pi):
            ego_turns[i] = 2
        else:
            ego_turns[i] = 1
    ego_turns = np.reshape(ego_turns, (len(idx3), 1))
            
    directions = np.diff(RList[:, 1:3], axis=0)
    allo = np.arctan2(directions[idx2_1, 1], directions[idx2_1, 0]) #in RList frame, length of idx3
    allo_turns = np.empty(len(idx3))
    for i in range(len(idx3)):
        if (np.pi/4) <= allo[i] < (3/4*np.pi):
            allo_turns[i] = 0
        elif (3/4*np.pi) <= allo[i] or (-3/4*np.pi) > allo[i]:
            allo_turns[i] = 3
        elif (-3/4*np.pi) <= allo[i] < (-1/4*np.pi):
            allo_turns[i] = 2
        else:
            allo_turns[i] = 1
    allo_turns = np.reshape(allo_turns, (len(idx3), 1))
    
    times = np.reshape(RList[idx3, 0], (len(idx3), 1))
    return np.hstack((times, ego_turns, allo_turns))


def writeToCsv(TimeAndTurns, title): 
    """
    from 052320
    get TimeAndTurns from classifyTurns
    """
    with open(title, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Timestamp", "Egocentric", "Allocentric"])
        csvwriter.writerows(TimeAndTurns)


def turnsInOutFields(visits, position, subfield):
    """
    from 061220: returns turns within 3s of visits to subfields (from Unit.visits)
        columns: ts, ego, allo, before/during/after visit
    """
    ranges = [] #list of lists with 3s before and after visits to subfields 
                #(not 3s of activity)
    for i in range(len(visits[subfield][0])): #the ith visit
        #start and end are indexes of start and end of +- 3s of visit
        start = bisect(position[:,0], (visits[subfield][0][i][0] - 3e6))
        end = bisect(position[:,0], (visits[subfield][0][i][-1] + 3e6))
        ranges.append([start, end])
    
    turns = np.empty((0, 4)) #columns: ts, ego, allo, before/during/after visit
    for i in range(len(ranges)): #the ith visit
        turns1 = classifyTurns(ranges[i][0], ranges[i][1], position)
        bda = np.empty((0,1))
        for j in range(len(turns1)): #the jth turn
            if turns1[j][0] < visits[subfield][0][i][0]:
                bda = np.vstack((bda,0))
            elif turns1[j][0] > visits[subfield][0][i][-1]:
                bda = np.vstack((bda,2))
            else:
                bda = np.vstack((bda,1))
        turns1 = np.hstack((turns1, bda))
        turns = np.vstack((turns, turns1))
    return turns


#from 061220: historical, 3D histogram and 2D histogram with heatmap
    

def makeCustomColormap(nb=100,name='mymap',c=[]):
    """
    From Will's utility_fx
    """
    if c ==[]:
        c = [(0,0,0.5),(1,1,0),(1,0,0)] #dark blue, yellow, red
    mycm = LinearSegmentedColormap.from_list(name,c,N=nb)
    return mycm


def graphRM(position, spikes, title, percentile=95, perim=[], mycm="jet"):
    n = makeRM(spikes, position)
    vmax = np.nanpercentile(n, percentile)
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title + "\n" + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz")
    im = ax.imshow(n, cmap=mycm, origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x coordinates (cm)")
    ax.set_ylabel("y coordinates (cm)")
    colors = ["b","g","r","k","c","m","y"]
    for i in range(len(perim)):
        ax.plot(perim[i][:,0]/4.72, perim[i][:,1]/4.72, color=colors[i], 
                label=f"Subfield {i}")
    ax.legend()
    cb = fig.colorbar(im)
    cb.set_label("Rate (Hz)")


def adjustPosCamera(datafile, pos):
    """
    Camera can be oriented differently w.r.t track across days/rats
    if things get rotated. This reads in a .txt file that should be 
    in each day's data dir saying whether x or y should be flipped
    format is e.g. x:640\ny:None
    
    datafile: path to the data file
    pos: position array with ts,x,y
    """
    with open(datafile+"cameraOrientationInfo.txt","r") as cam:
        info = cam.readlines()
    info = [i.rstrip() for i in info]
    info = {info[0][0]:info[0][2:], info[1][0]:info[1][2:]}
    if info['x'] != 'None':
        posx = [int(info['x'])-pos[i][1] for i in range(len(pos))]
    else:
        posx = [pos[i][1] for i in range(len(pos))]

    if info['y'] != 'None':
        posy = [int(info['y'])-pos[i][2] for i in range(len(pos))]
    else:
        posy = [pos[i][2] for i in range(len(pos))]
    
    ts = np.reshape(pos[:,0], (len(pos),1))
    posx = np.reshape(posx, (len(pos), 1))
    posy = np.reshape(posy, (len(pos), 1))
    posTsXY = np.hstack((ts,posx,posy))
    return posTsXY


#from 061820: historical, bar graphs of firing rates before and after turns
    

def shiftPos(pos, timeAndTurns, spikes, RMrange=7*4.72):
    """
    pos: [ts,x,y] before RDP
    Rpos: [ts,x,y] after RDP
    idx: index of the first point of a turn
    spikes: after getPosFromTs
    RMrange: +/- n unit distance of the 1st point of turns to make rate map
    
    Returns: pos (before RDP) and spikes within RMrange of the 1st point of 
    turns, shifted based on the 1st point of turns
    """
    shifted_pos = [np.empty((0,3)) for _ in range(8)]
    shifted_s = [np.empty((0,3)) for _ in range(8)]
    for i in range(len(timeAndTurns)): #adapted from turn3
        k = np.where(pos == timeAndTurns[i,0])[0] #1st point of the turn, in pos frame
        j = deepcopy(k)
        dist = 0
        #idx2: which pts are within range of 1st pt of turns, in pos frame
        while dist < RMrange: #checking points after the turn
            idx2_end = deepcopy(j)
            j += 1
            if j > len(pos)-1:
                break
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
        
        j = k - 1 #one point before the turn
        if j >= 0:
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
            idx2_start = j+1
        else:
            dist = RMrange+10
            idx2_start = 0
        while dist < RMrange: #checking points before the turn
            idx2_start = deepcopy(j)
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((pos[k, 1]-pos[j, 1])**2 + (pos[k, 2]-pos[j, 2])**2)
        
        shiftedP = np.subtract(pos[int(idx2_start):int(idx2_end)], np.array([0,pos[k,1],pos[k,2]]))
        #i_s = np.empty((len(shiftedP),1))
        #i_s.fill(i)
        #shiftedP = np.hstack((shiftedP, i_s))
        #shifted_pos = np.vstack((shifted_pos, shiftedP))
        
        spike_idx = np.where(np.logical_and(spikes>pos[idx2_start,0], spikes<pos[idx2_end,0]))[0]
        shiftedS = np.subtract(spikes[spike_idx], np.array([0,pos[k,1],pos[k,2]]))
        #i_s = np.empty((len(shiftedS),1))
        #i_s.fill(i)
        #shiftedS = np.hstack((shiftedS, i_s))
        #shifted_s = np.vstack((shifted_s, shiftedS))
        
        shifted_pos[0] = np.vstack((shifted_pos[0], shiftedP))
        shifted_s[0] = np.vstack((shifted_s[0], shiftedS))
        #ego turns
        if timeAndTurns[i,1] == 0:
            shifted_pos[1] = np.vstack((shifted_pos[1], shiftedP))
            shifted_s[1] = np.vstack((shifted_s[1], shiftedS))
        elif timeAndTurns[i,1] == 1:
            shifted_pos[2] = np.vstack((shifted_pos[2], shiftedP))
            shifted_s[2] = np.vstack((shifted_s[2], shiftedS))
        elif timeAndTurns[i,1] == 2:
            shifted_pos[3] = np.vstack((shifted_pos[3], shiftedP))
            shifted_s[3] = np.vstack((shifted_s[3], shiftedS))
        
        #allo turns
        if timeAndTurns[i,2] == 0:
            shifted_pos[4] = np.vstack((shifted_pos[4], shiftedP))
            shifted_s[4] = np.vstack((shifted_s[4], shiftedS))
        elif timeAndTurns[i,2] == 1:
            shifted_pos[5] = np.vstack((shifted_pos[5], shiftedP))
            shifted_s[5] = np.vstack((shifted_s[5], shiftedS))
        elif timeAndTurns[i,2] == 2:
            shifted_pos[6] = np.vstack((shifted_pos[6], shiftedP))
            shifted_s[6] = np.vstack((shifted_s[6], shiftedS))
        elif timeAndTurns[i,2] == 3:
            shifted_pos[7] = np.vstack((shifted_pos[7], shiftedP))
            shifted_s[7] = np.vstack((shifted_s[7], shiftedS))
    return shifted_pos, shifted_s


def makeRM2(spikes, position, smoothing_2d_sigma=2):
    """ makeRM adapted for avg RM for turns"""
    #camXMax, camXMin = [round(np.max(position[:,1])), round(np.min(position[:,1]))]
    #camYMax, camYMin = [round(np.max(position[:,2])), round(np.min(position[:,2]))]
    #rbins = int((camXMax-camXMin)/4.72)
    #cbins = int((camYMax-camYMin)/4.72)
    rows, cols = np.linspace(-7*4.72, 7*4.72, 15), np.linspace(-7*4.72, 7*4.72, 15)
    hs,xe,ye = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rows, cols])
    ho = np.histogram2d(position[:,2],position[:,1],bins=[rows, cols])[0]
    n = (hs*np.reciprocal(ho))*30 #adjust for camera 30 frames/sec to get Hz
    n[np.where(ho==0)] = np.nan
    n = weird_smooth(n,smoothing_2d_sigma)
    n[np.where(ho==0)] = np.nan
    #cam = [floor(camXMin/4.72),ceil(camXMax/4.72),floor(camYMin/4.72),ceil(camYMax/4.72)]
    return n

def graphRM2(position, spikes, suptitle, percentile=95, mycm="jet"):
    """ 
    graphRM adapted for avg RM for turns
    position and spikes from shiftPos
    """
    n = []
    for i in range(8):
        n.append(makeRM2(spikes[i], position[i]))
    
    fig, axs = plt.subplots(2,4)
    vmax = np.nanpercentile(n[0], percentile)
    fig.suptitle(suptitle + "\n" + f"Cutoff = {percentile}th percentile, {round(vmax,1)} Hz", y=1.08)
    titles = ["All turns", "Left", "Back", "Right", "North", "East", "South", "West"]
    ims = []
    for i in range(2):
        for j in range(4):
            axs[i][j].set_title(titles[i*4+j])
            ims.append(axs[i][j].imshow(n[i*4+j], cmap=mycm, origin="lower", vmin=0, vmax=vmax))
            axs[i][j].axis("equal")
    cb = fig.colorbar(ims[0])
    cb.set_label("Rate (Hz)")
    fig.tight_layout()