# -*- coding: utf-8 -*-
"""

"""
import numpy as np

def PerpDist(p1, p2, p3):
    return abs(np.cross(p2 - p1, p3 - p1)/np.linalg.norm(p1 - p2))

class RDP2:
    
    ResultList = [0,0]
    
    def RDPfunc(self, PointList, epsilon, Results = []):
        dmax = 0
        
        d = PerpDist(PointList[0], PointList[-1], PointList)
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

