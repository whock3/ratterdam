# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:23:49 2020

@author: Ruo-Yah Lai
"""
import numpy as np
from timeit import timeit

def PerpDist(p1, p2, p3):
    return abs(np.cross(p2 - p1, p3 - p1)/np.linalg.norm(p1 - p2))

def RDPnew(PointList, epsilon, ResultList):
        dmax = 0
        
        d = PerpDist(PointList[0, 1:3], PointList[-1, 1:3], PointList[:, 1:3])
        for i in range(1, len(PointList)-1):
            if d[i] > dmax:
                dmax = d[i]
                index = i

        if dmax > epsilon:
            ResultList = RDPnew(PointList[:index+1], epsilon, ResultList)
            ResultList = RDPnew(PointList[index:], epsilon, ResultList)
        else:
            ResultList = np.vstack([ResultList, PointList[-1]])

        return ResultList