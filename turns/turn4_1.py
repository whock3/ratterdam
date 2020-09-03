# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 2020

@author: Ruo-Yah Lai

1. idx = Points that are > min_angle1
2. idx2 = Points within turn_range of a point in idx
      a. While taking out points in idx2 from idx
      b. Turn range depends on instantaneous velocity and how
      many pts are removed with RDP (min turn range = 2cm)
3. theta_sum = Add up angles of groups of points in idx2
4. idx3 = where theta_sum > min_angle2
5. NEW - If there is back and forth movement both before and after the 
   possible turn, then it's not a turn
6. Returns first point in a group
"""

import numpy as np
from copy import deepcopy
from bisect import bisect

def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)

def turn4(position, min_angle1=np.pi/12, min_angle2=np.pi/4):
    """
    position: [ts, x, y] after RDP
    
    """
    sAhead = 20 #how many seconds ahead to look when calculating %pts kept after RDP
    raw_radius = 1/20 #in Hz, from 30 frames/s and 3 frames for a turn
    
    diffts, diffx, diffy = np.diff(position[:,0]), np.diff(position[:,1]), np.diff(position[:,2]) 
    vx = 1e6*diffx/diffts
    vy = 1e6*diffy/diffts
    v =  np.sqrt((vx**2)+(vy**2))

    directions = np.diff(position[:,1:], axis=0)
    theta2 = angle2(directions) #[0, 2pi]
    theta = deepcopy(theta2)
    for i in range(len(theta2)):
        if theta[i] > np.pi:
            theta[i] = 2*np.pi - theta2[i] #[0, pi]
    idx = np.where(theta>min_angle1)[0]+1 #in Rpos frame
    idx_opp = np.where(theta<=min_angle1)[0]+1 #in Rpos frame
    idx, idx_opp = list(idx), list(idx_opp)

    theta_sum2 = np.empty(0) #[0, 2pi]
    idx2_1 = np.empty(0)
    idx_1 = np.empty(0)
    trs = np.empty(0)
    percents = np.empty(0)
    while len(idx) > 0:
        i = idx[0] #in Rpos frame
        #percentage of points kept after RDP
        a = bisect(position[:,0], (position[i,0] + sAhead*1e6))
        if a == len(position):
            a = a-1
        percent = (a-i) / ((position[a,0]-position[i,0])/1e6*30)
        turn_range = raw_radius * v[i] / percent
        if turn_range < 2*4.72:
            turn_range = 2*4.72
            
        j = i-1 #the "turn" in question, in theta frame
        dist = 0
        idx2 = np.array([]) #which pts are within range of "turns", in theta frame
        while dist < turn_range: #checking points after the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j += 1
            if j > len(position)-3:
                break
            dist = dist + np.sqrt((position[i, 1]-position[j+1, 1])**2 + 
                                  (position[i, 2]-position[j+1, 2])**2)
        
        j = i-2 #one point before the "turn", in theta frame
        if j >= 0:
            dist = np.sqrt((position[i, 1]-position[j+1, 1])**2 + 
                           (position[i, 2]-position[j+1, 2])**2)
        else:
            dist = turn_range+10
        while dist < turn_range: #checking points before the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j -= 1
            if j < 0:
                break
            dist = dist + np.sqrt((position[i, 1]-position[j+1, 1])**2 + 
                                  (position[i, 2]-position[j+1, 2])**2)
            
        idx2 = idx2.astype(int)
        for k in idx2:
            if k+1 in idx:
                idx.remove(k+1)
            elif k+1 in idx_opp:
                idx_opp.remove(k+1)
            else:
                idx2 = list(idx2)
                idx2.remove(k)
                idx2 = np.array(idx2)
        
        idx_1 = np.hstack([idx_1, np.min(idx2)+1]) #1st pt in each "turn"
        theta_sum2 = np.hstack([theta_sum2, np.sum(theta2[idx2]) % (2*np.pi)])
        idx2_1 = np.hstack([idx2_1, np.max(idx2)+1]) #last point in each group of "turns", in Rpos frame
        trs = np.hstack((trs, turn_range))
        percents = np.hstack((percents, percent))
        #print(idx_1[-1]-1, idx2, round(turn_range/4.72,1), round(v[i]/4.72,1), round(theta_sum2[-1]/np.pi*180,1))
        
                
    theta_sum = np.empty(len(theta_sum2)) #[0, pi]
    for i in range(len(theta_sum2)):
        if theta_sum2[i] > np.pi:
            theta_sum[i] = 2*np.pi - theta_sum2[i] #[0, pi]
        else:
            theta_sum[i] = theta_sum2[i]
    idx3 = np.where(theta_sum>min_angle2)[0]
    
    window = 20*4.72
    bTheta_sum, aTheta_sum = np.empty(0), np.empty(0)
    for i in idx3: #checking for back and forth movements
        j = idx2_1[i] #the 1st pt after the last pt of "turn" in question, in theta frame
        if j < len(position)-2:
            dist = np.sqrt((position[idx2_1[i], 1]-position[j+1, 1])**2 + 
                           (position[idx2_1[i], 2]-position[j+1, 2])**2)
        else:
            dist = window+10
        aidx2 = np.array([]) #which pts are within window of and after "turn", in theta frame
        while dist < window: #checking points after the "turn"
            aidx2 = np.concatenate((aidx2, np.array([j])))
            j += 1
            if j > len(position)-3:
                break
            dist = np.sqrt((position[idx2_1[i], 1]-position[j+1, 1])**2 + 
                           (position[idx2_1[i], 2]-position[j+1, 2])**2)
        
        
        j = idx_1[i]-2 #the 1st pt before the 1st pt of "turn", in theta frame
        if j >= 0:
            dist = np.sqrt((position[idx_1[i], 1]-position[j+1, 1])**2 + 
                           (position[idx_1[i], 2]-position[j+1, 2])**2)
        else:
            dist = window+10
        bidx2 = np.array([]) #which pts are within window of and before "turn", in theta frame
        while dist < window: #checking points before the "turn"
            bidx2 = np.concatenate((bidx2, np.array([j])))
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((position[idx_1[i], 1]-position[j+1, 1])**2 + 
                           (position[idx_1[i], 2]-position[j+1, 2])**2)
        
        bTheta_sum = np.hstack([bTheta_sum, np.sum(theta[bidx2])])
        aTheta_sum = np.hstack([aTheta_sum, np.sum(theta[aidx2])])
    
    #large angles both before and after, in idx frame
    widx = idx3[np.where(np.logical_and(bTheta_sum > np.pi, aTheta_sum > np.pi))[0]] 
        
    return idx_1[idx3].astype(int), theta_sum2[idx3], idx2_1[idx3].astype(int), trs[idx3], percents
