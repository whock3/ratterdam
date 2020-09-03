# -*- coding: utf-8 -*-
"""

"""
import RateMapClass_William_20190308 as RateMapClass
from williamDefaults import binWidth
from scipy.spatial import ConvexHull
from matplotlib import path
import numpy as np
from RateMap import weird_smooth


def getVisits(data, maxgap=2*1e6):
    """
    """
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

class Unit():
    """
    Wrapper class because rep field ID algorithm looks
    for instance.spikes and instance.position
    """
    
    def __init__(self, s, p):
        self.spikes = s
        self.position = p
        self.fields = []
        self.visits = [] # nested list. each list is a subfield and values are themselves lists of points in visit
        self.perimeters = []
        #self.colors = cnames
        self.smoothing = 2
        self.repUnit = RateMapClass.RateMap(self)
        self.alphaHullFactor = 1
        self.alpha = 0
        self.findFields()
        
    def PolyArea(self,x,y):
        """
        Found at https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
        """
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
    def findFields(self):
        self.fields = []
        self.visits = []
        for i,pf in enumerate(self.repUnit.PF[:]):
            
            #create boundary using alphahull alg which allows for concave hulls but does not work all that well as is
#             alpha = self.alphaHullFactor*alphashape.optimizealpha(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)))
#             hull = alphashape.alphashape(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)),alpha)
#             hxy = hull.exterior.coords.xy
#             contour = path.Path(list(zip(hxy[0],hxy[1])))

 

            #create boundary using convex hull
            points = np.asarray(list(zip(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)))
            hull = ConvexHull(points)
            vertices = np.append(hull.vertices, hull.vertices[0]) # add the first point to close the contour
            contour = path.Path(points[vertices])

 

            PinC = self.position[contour.contains_points(self.position[:,1:])]
            posVisits = getVisits(PinC[:,0])
            self.visits.append(posVisits)
            field_FR = []
            field_TS = [] # take middle ts 
            for visit in posVisits:
                spk = self.spikes[np.logical_and(self.spikes[:,0] > visit[0], self.spikes[:,0] < visit[-1])]
                vdur = (visit[-1]-visit[0])/1e6 # do this instead of just multiplying number samples by fr because any lost frames from, e.g. occlusions will introduce a discrepancy
                field_FR.append(spk.shape[0]/vdur)
                field_TS.append(((visit[0]-self.position[0,0])/1e6)/60)
                
            field_FR = weird_smooth(np.asarray(field_FR), self.smoothing)
            
            totalRate = sum(field_FR)
            area = self.PolyArea(pf.perimeter[1]*binWidth+binWidth/2, pf.perimeter[0]*binWidth+binWidth/2)
            print(i, totalRate/area)
            if True:
                self.fields.append(np.column_stack((field_TS, field_FR)))
                self.perimeters.append(points[vertices])
