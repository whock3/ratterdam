# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 20:29:14 2018

@author: whock
"""

alley = alleyTracking[x]

def getVisits(data, maxgap=3*1e6):
    data.sort()
    groups = [[data[0]]]
    idx = [[0]]
    for i,x in enumerate(data[1:]):
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
            idx[-1].append(i)
        else:
            groups.append([x])
            idx.append([i])
            
    for v in range(len(idx)):
        idx[v] = np.asarray(idx[v])
            
    return idx,groups

pvs = getVisits(alley[:,2],1e6)
v1 = alley[pvs[0]]

bounds = Dir.extractCorners(alleyBounds[x-1])
ul,ll,ur,lr = bounds
def orien = checkOrientation(alley):
    if abs(ll[0]-lr[0]) > abs(ll[1] - ul[1]):
        mp = (ul[0]+ur[0])/2
    else:
        mp = mp = (ul[1]+ur[1])/2
    return mp



mins = [abs(mp-i) for i in v[:,1]]
if min(mins) < epsilon:
    visitPass = 1
    
# check if it reaches other end too, less chance occulusion will burn you
# how to deal with that in general (occlusions or lack of smapling?) can interpolate :( 