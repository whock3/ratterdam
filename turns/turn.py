# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from copy import deepcopy

def angle2(dir):
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return (np.arctan2(dir2[:, 1], dir2[:, 0]) - np.arctan2(dir1[:, 1], dir1[:, 0])) % (2*np.pi)

def turn(position, min_angle1, min_angle2, turn_range):
    directions = np.diff(position, axis=0)
    theta2 = angle2(directions) #[0, 2pi]
    theta = deepcopy(theta2)
    for i in range(len(theta2)):
        if theta[i] > np.pi:
            theta[i] = 2*np.pi - theta2[i] #[0, pi]
    idx = np.where(theta>min_angle1)[0]+1 #in RList/position frame

    theta_sum2 = np.empty(len(idx)) #[0, 2pi]
    theta_sum = np.empty(len(idx)) #[0, pi]
    idx2_1 = np.empty(len(idx))
    for i in range(len(idx)):
        
        j = idx[i]-1 #the "turn" in question, in theta frame
        dist = 0
        idx2 = np.array([]) #which pts are within range of "turns", in theta frame
        while dist < turn_range: #checking points after the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j += 1
            if j > len(position)-3:
                break
            dist = np.sqrt((position[idx[i], 0]-position[j+1, 0])**2 + (position[idx[i], 1]-position[j+1, 1])**2)
        
        j = idx[i]-2 #one point before the "turn", in theta frame
        if j >= 0:
            dist = np.sqrt((position[idx[i], 0]-position[j+1, 0])**2 + (position[idx[i], 1]-position[j+1, 1])**2)
        else:
            dist = turn_range+10
        while dist < turn_range: #checking points before the "turn"
            idx2 = np.concatenate((idx2, np.array([j])))
            j -= 1
            if j < 0:
                break
            dist = np.sqrt((position[idx[i], 0]-position[j+1, 0])**2 + (position[idx[i], 1]-position[j+1, 1])**2)
        
        idx2 = idx2.astype(int)
        theta_sum2[i] = np.sum(theta2[idx2]) % (2*np.pi)
        if theta_sum2[i] > np.pi:
            theta_sum[i] = 2*np.pi - theta_sum2[i] #[0, pi]
        else:
            theta_sum[i] = theta_sum2[i]
        
        idx2_1[i] = np.max(idx2)+1 #last point in each group of "turns", in RList frame
    idx3 = np.where(theta_sum>min_angle2)
    
    return idx[idx3], theta_sum2[idx3], idx2_1[idx3].astype(int)

#for j in len(idx2):
    #idx2_1 = np.where(idx2 == idx[i]-1)[0]
    #idx2_1 = np.concatenate((idx2_1, np.where(abs(idx2 - (idx[i]-1)) == 1)[0]))