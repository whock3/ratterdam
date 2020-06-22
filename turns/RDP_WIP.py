# -*- coding: utf-8 -*-
"""
PointList = ([[x1, y1],
              [x2, y2],
              ...
              [xn, yn]])
"""

import numpy as np
from pandas import read_excel

#import data from excel
data = read_excel(r"C:\Users\Ariel Lai\Desktop\My folder\College\Junior\Sample trajectory.xlsx", headers=None)
data_x = np.array(data)

#perpendicular distance between point p3 and line with points p1 and p2
def PerpDist(p1, p2, p3):
        return abs(np.cross(p2 - p1, p3 - p1)/np.linalg.norm(p1 - p2))

ResultList = data_x[0]

def RDP(PointList, epsilon, Results = []):
    dmax = 0
    global ResultList
      
    for i in range(1, len(PointList)-1):
        d = PerpDist(PointList[0], PointList[-1], PointList[i])
        if d > dmax:
            dmax = d
            index = i
    
    if dmax > epsilon:
        RDP(PointList[:index+1], epsilon, ResultList)
        RDP(PointList[index:], epsilon, ResultList)
    else:
        ResultList = np.vstack([Results, PointList[-1]])

    return ResultList



