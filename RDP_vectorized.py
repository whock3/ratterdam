# -*- coding: utf-8 -*-
"""
data = ([[x1, y1],
         [x2, y2],
         ...
         [xn, yn]])
"""

import numpy as np
from velocity_filtering import velocity_filtering
import csv
with open("pos.p.ascii", 'r') as csv_file:
    data_iter = csv.reader(csv_file)
    data = [data for data in data_iter]
data_np1 = np.array(data[50574:-17], dtype = "float64")
a = np.array([1, 4.72, 4.72])
data_np2 = data_np1[:,0:3]/a[None,:]
data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]
data_np4 = velocity_filtering(data_np3)

#perpendicular distance between point p3 and line with points p1 and p2
def PerpDist(p1, p2, p3):
    return abs(np.cross(p2 - p1, p3 - p1)/np.linalg.norm(p1 - p2))

ResultList = data_np4[500, 1:3]

def RDP(PointList, epsilon, Results = []):
    dmax = 0
    global ResultList
        
    d = PerpDist(PointList[0], PointList[-1], PointList)
    for i in range(1, len(PointList)-1):
        if d[i] > dmax:
            dmax = d[i]
            index = i

    if dmax > epsilon:
        RDP(PointList[:index+1], epsilon, ResultList)
        RDP(PointList[index:], epsilon, ResultList)
    else:
        ResultList = np.vstack([Results, PointList[-1]])

    return ResultList

