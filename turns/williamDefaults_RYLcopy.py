#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: William Snider
"""
import numpy as np
import pandas as pd

## Default values to be used for all ratterdam scripts, coded here for consistency

## Binning values and size of arena
# Values are now hardcoded, commented out lines calculated more precisely for individual plots
binWidth = 10
xmax = 650 #np.max([np.max(unit.spikes[:,1]),np.max(unit.position[:,1])])
ymax = 500 #np.max([np.max(unit.spikes[:,2]),np.max(unit.position[:,2])])
xedges = np.arange(0, xmax+binWidth, binWidth) #np.arange(0, int(np.ceil(xmax/binWidth)*(binWidth+1)), binWidth)
yedges = np.arange(0, ymax+binWidth, binWidth) #np.arange(0, int(np.ceil(ymax/binWidth)*(binWidth+1)), binWidth)

## Gaussian Filter for 2D Rate map
sigma = 1.5

velocity_filter_thresh = 1 #TODO: play this with value #units are cm/s

# %%format is [x1, x2], [y1, y2] #??? Adjusted intersection 2,3,5,7 from y2=315 to y2=310 to fix bin alignment
# Responsive coordinates: v1, v2, ... are the 1st, 2nd vertical lines defining the bounds, same for h0,h1,h2...
v1, v2, v3, v4, v5, v6, v7, v8 = [100, 150, 240, 290, 374, 426, 505, 554]
h1, h2, h3, h4, h5, h6 = [70, 130, 205, 266, 340, 400]


cityBlocks = np.array([[[v2,h2], [v3,h2], [v3,h3], [v2,h3], [v2,h2]],
                       [[v2,h4], [v3,h4], [v3,h5], [v2,h5], [v2,h4]],
                       [[v4,h2], [v5,h2], [v5,h3], [v4,h3], [v4,h2]],
                       [[v4,h4], [v5,h4], [v5,h5], [v4,h5], [v4,h4]],
                       [[v6,h2], [v7,h2], [v7,h3], [v6,h3], [v6,h2]],
                       [[v6,h4], [v7,h4], [v7,h5], [v6,h5], [v6,h4]]])


alleyLines = np.array([[[v1,h1], [v1,h6]],
                       [[v2,h1], [v2,h6]],
                       [[v3,h1], [v3,h6]],
                       [[v4,h1], [v4,h6]],
                       [[v5,h1], [v5,h6]],
                       [[v6,h1], [v6,h6]],
                       [[v7,h1], [v7,h6]],
                       [[v8,h1], [v8,h6]],
                       [[v1,h1], [v8,h1]],
                       [[v1,h2], [v8,h2]],
                       [[v1,h3], [v8,h3]],
                       [[v1,h4], [v8,h4]],
                       [[v1,h5], [v8,h5]],
                       [[v1,h6], [v8,h6]],])


alleyInterBounds = {'0':[[v2, v3],[h5, h6]],
              '1':[[v2, v3],[h3, h4]],
              '2':[[v1, v2],[h4, h5]],
              '3':[[v3, v4],[h4, h5]],
              '4':[[v4, v5],[h5, h6]],
              '5':[[v5, v6],[h4, h5]],
              '6':[[v6, v7],[h5, h6]],
              '7':[[v7, v8],[h4, h5]],
              '8':[[v6, v7],[h3, h4]],
              '9':[[v7, v8],[h2, h3]],
              '10':[[v6, v7],[h1, h2]],
              '11':[[v5, v6],[h2, h3]],
              '12':[[v4, v5],[h3, h4]],
              '13':[[v4, v5],[h1, h2]],
              '14':[[v3, v4],[h2, h3]],
              '15':[[v2, v3],[h1, h2]],
              '16':[[v1, v2],[h2, h3]],
              'A':[[v1, v2],[h5,h6]],
              'B':[[v3, v4],[h5,h6]],
              'C':[[v5, v6],[h5,h6]],
              'D':[[v7, v8],[h5,h6]],
              'E':[[v1, v2],[h3,h4]],
              'F':[[v3, v4],[h3,h4]],
              'G':[[v5, v6],[h3,h4]],
              'H':[[v7, v8],[h3,h4]],
              'I':[[v1, v2],[h1,h2]],
              'J':[[v3, v4],[h1,h2]],
              'K':[[v5, v6],[h1,h2]],
              'L':[[v7, v8],[h1,h2]]}

alleyInterType = {'A':['intersection'],
                  'B':['intersection'],
                  'C':['intersection'],
                  'D':['intersection'],
                  'E':['intersection'],
                  'F':['intersection'],
                  'G':['intersection'],
                  'H':['intersection'],
                  'I':['intersection'],
                  'J':['intersection'],
                  'K':['intersection'],
                  'L':['intersection'],
                  '0':['alley','horizontal'],
                  '1':['alley','horizontal'],
                  '4':['alley','horizontal'],
                  '6':['alley','horizontal'],
                  '8':['alley','horizontal'],
                  '10':['alley','horizontal'],
                  '12':['alley','horizontal'],
                  '13':['alley','horizontal'],
                  '15':['alley','horizontal'],
                  '2':['alley','vertical'],
                  '3':['alley','vertical'],
                  '5':['alley','vertical'],
                  '7':['alley','vertical'],
                  '9':['alley','vertical'],
                  '11':['alley','vertical'],
                  '14':['alley','vertical'],
                  '16':['alley','vertical']}

# midpoints of alleys and intersections
alleyInterSkeleton = {
        'sV1': v1+(v2-v1)/2,
        'sV2': v3+(v4-v3)/2,
        'sV3': v5+(v6-v5)/2,
        'sV4': v7+(v8-v7)/2,
        'sH1': h1+(h2-h1)/2,
        'sH2': h3+(h4-h3)/2,
        'sH3': h5+(h6-h5)/2,
        }

skelPts = []
for h in np.arange(alleyInterSkeleton['sH1'],alleyInterSkeleton['sH3']+1,5):
    skelPts.append([h, alleyInterSkeleton['sV1']])
    skelPts.append([h, alleyInterSkeleton['sV2']])
    skelPts.append([h, alleyInterSkeleton['sV3']])
    skelPts.append([h, alleyInterSkeleton['sV4']])
for v in np.arange(alleyInterSkeleton['sV1'],alleyInterSkeleton['sV4']+1,5):
    skelPts.append([alleyInterSkeleton['sH1'], v])
    skelPts.append([alleyInterSkeleton['sH2'], v])
    skelPts.append([alleyInterSkeleton['sH3'], v])
skelPts = np.array(skelPts)

# For skeleton smoothing
runningMeanN = 5

turnDirection = {}
# %%
#for a in alleyInterBounds:
#    
#    # For intersections (letters), skip (continue)
#    if a.isdigit() is False:
#        continue
#        
#    for b in alleyInterBounds:
#        
#        if a == b or b.isdigit() is False:
#            continue
#        
#        [[ax1, ax2], [ay1, ay2]] = alleyInterBounds[a]
#        [[bx1, bx2], [by1, by2]] = alleyInterBounds[b]
#        
#        # Test if alleys are adjacent; share y value and x values
##        equalArr = np.array([x==y for x in [ax1, ax2, ay1, ay2] for y in [bx1, bx2, by1, by2]]).reshape([4,4])
##        df = pd.DataFrame(data = equalArr, index = ['bx1', 'bx2', 'by1', 'by2'], columns = ['ax1', 'ax2', 'ay1', 'ay2'])
#        
#        
#        # Test if alleys are across from each other
#        fudgeFactor = 75 # bigger than width of one intersection, but smaller than width of intersection+alley
#
#        # share y values
#
#        shareY = np.any([ay1 == by1, ay2 == by1, ay2 == by2, ay1 == by2])
#        shareX = np.any([ax1 == bx1, ax2 == bx1, ax2 == bx2, ax1 == bx2])      
#        across = np.any([
#                        np.all([ay1 == by1, ay2 == by2, np.any([np.abs(ax2 - bx1) < fudgeFactor, np.abs(ax1 - bx2) < fudgeFactor]), ay2-ay1 < ax2-ax1, by2-by1 < bx2-bx1]), # horizontal
#                        np.all([ax1 == bx1, ax2 == bx2, np.any([np.abs(ay2 - by1) < fudgeFactor, np.abs(ay1 - by2) < fudgeFactor]), ay2-ay1 > ax2-ax1, by2-by1 > bx2-bx1])]) # vertical
#        
#                
#            
#        if np.all([shareY, shareX]) or across:
#            
#            # Find intersection
#            for c in alleyInterBounds:
#                
#                if a.isdigit() is True:
#                    continue
#                
#                [[cx1, cx2], [cy1, cy2]] = alleyInterBounds[c]
#                
#                xIn
#                yIn
#            pprint.pprint(a+' '+b+' share X and Y value')
#            
#        # share x values
##        if [ax1 == bx1 or ax2 == bx1 or ax2 == bx2 or ax1 == bx2]:
##            print(a+' '+b+' share X value')

        
# %%       

# alleybounds used so Will's code doesn't break. These values are wrong.       
alleyBounds = {0:[[165, 210],[345, 385]],
                   1:[[150, 215],[210, 250]],
                   2:[[105, 140],[272, 340]],
                   3:[[230, 270],[275, 340]],
                   4:[[295, 359],[352, 390]],
                   5:[[370, 410],[260, 327]],
                   6:[[430, 492],[347, 390]],
                   7:[[500, 540],[260, 325]],
                   8:[[420, 480],[210, 254]],
                   9:[[500, 550],[140, 205]],
                   10:[[415, 480],[75, 115]],
                   11:[[370, 410],[135, 205]],
                   12:[[280, 350],[210, 250]],
                   13:[[280, 350],[75, 115]],
                   14:[[230, 275],[132, 205]],
                   15:[[165, 230],[75, 115]],
                   16:[[105, 145],[115, 185]]}   
# %%

cameraFrameRate = 30# fps. Used to give accurate spiking rate.

## Algorithm for detecting place fields
min_PF_proportion = 0.02 #(approximately 22 bins under testing setup)
prctile = 95 # 
cutoffFromPeakValue = 0.2 # A field is determined by its peak and the surrounding area including bins that are at least 20% of the max firing rate
minPeakFiringRate = 2 # #Fields have to fire more than 2Hz at their peak
minFiringRateFactor = 0; # Removed at Jims request. Place field peaks must not be above a certain fraction of the global peak. Formerly 0.5.

## PF Splitting algorithm parameters
minPeakDistanceAwayPixels = 40 #bins #in pixel units 
minPeakDistanceAwayBins = int(np.round(minPeakDistanceAwayPixels/binWidth))
minValleyProportion = 0.75

# Find PF Firing Over Time
timeStampScale = 1000000 # neuralynx is in microseconds
minPFOccupancyTime = 0

# Visits to Place Field
minGapLength = 0.25 # units in s

## Calculated parameters

# For finding alley/intersection visits
minVisitSizeAI = 0.125