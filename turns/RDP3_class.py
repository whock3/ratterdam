# -*- coding: utf-8 -*-
"""
Created on Sun May  3 20:23:01 2020

@author: Ariel Lai
"""
import numpy as np

def PerpDist(p1, p2, p3):
    return abs(np.cross(p2 - p1, p3 - p1)/np.linalg.norm(p1 - p2))

#keeps timestamp
class RDP3:
    
    ResultList = [0,0,0]
    
    def RDPfunc(self, PointList, epsilon, Results = []):
        dmax = 0
        
        d = PerpDist(PointList[0, 1:3], PointList[-1, 1:3], PointList[:, 1:3])
        for i in range(1, len(PointList)-1):
            if d[i] > dmax:
                dmax = d[i]
                index = i

        if dmax > epsilon:
            self.RDPfunc(self, PointList[:index+1], epsilon, self.ResultList)
            self.RDPfunc(self, PointList[index:], epsilon, self.ResultList)
        else:
            self.ResultList = np.vstack([Results, PointList[-1]])

        return self.ResultList

