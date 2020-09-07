# -*- coding: utf-8 -*-
"""
1. idx = Points that are > min_angle1
2. idx2 = Points within turn_range of a point in idx
      a. While taking out points in idx2 from idx
      b. NEW - Turn range depends on instantaneous velocity
3. theta_sum = Add up angles of groups of points in idx2
4. idx3 = where theta_sum > min_angle2
      a. NEW - returns first point in a group
"""
import numpy as np
from copy import deepcopy

def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)

def turn3(position, min_angle1, min_angle2, turn_range_start=5, vCutoff=15, trIncrease=1/7):
    #position has ts, x, y after RDP\
    #turn_range_start in cm(or camera), vCutoff in cm/s(or camera/s), trIncrease in s
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
    idx = np.where(theta>min_angle1)[0]+1 #in RList/position frame
    idx_opp = np.where(theta<=min_angle1)[0]+1 #in RList/position frame
    idx, idx_opp = list(idx), list(idx_opp)

    theta_sum2 = np.empty(0) #[0, 2pi]
    idx2_1 = np.empty(0)
    idx_1 = np.empty(0)
    idx_2 = np.empty(0)
    while True:
        if len(idx) < 1:
            break
        i = idx[0]
        
        if v[i] > vCutoff:
            turn_range = turn_range_start + (v[i]-vCutoff)*trIncrease
        else:
            turn_range = deepcopy(turn_range_start)
            
        j = i-1 #the "turn" in question, in theta frame
        dist = 0
        idx2 = np.array([]) #which pts are within range of "turns", in theta frame
        while dist < turn_range: #checking points after the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j += 1
            if j > len(position)-3:
                break
            dist = np.sqrt((position[i, 1]-position[j+1, 1])**2 + (position[i, 2]-position[j+1, 2])**2)
        
        j = i-2 #one point before the "turn", in theta frame
        if j >= 0:
            dist = np.sqrt((position[i, 1]-position[j+1, 1])**2 + (position[i, 2]-position[j+1, 2])**2)
        else:
            dist = turn_range+10
        while dist < turn_range: #checking points before the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((position[i, 1]-position[j+1, 1])**2 + (position[i, 2]-position[j+1, 2])**2)
        
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
        
        idx_1 = np.hstack([idx_1, np.min(idx2)+1])
        theta_sum2 = np.hstack([theta_sum2, np.sum(theta2[idx2]) % (2*np.pi)])
        idx2_1 = np.hstack([idx2_1, np.max(idx2)+1]) #last point in each group of "turns", in RList frame
        idx_2 = np.hstack([idx_2, np.min(idx2)]) #1 pt before 1st pt of turn
        #print(idx_1[-1]-1, idx2, round(turn_range/4.72,1), round(v[i]/4.72,1), round(theta_sum2[-1]/np.pi*180,1))
        
                
    theta_sum = np.empty(len(theta_sum2)) #[0, pi]
    for i in range(len(theta_sum2)):
        if theta_sum2[i] > np.pi:
            theta_sum[i] = 2*np.pi - theta_sum2[i] #[0, pi]
        else:
            theta_sum[i] = theta_sum2[i]
                        
    idx3 = np.where(theta_sum>min_angle2)[0]
        
    return idx_1[idx3].astype(int), theta_sum2[idx3], idx2_1[idx3].astype(int), idx_2[idx3].astype(int)