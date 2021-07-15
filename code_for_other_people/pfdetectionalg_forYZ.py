# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 12:01:44 2021

@author: whockei1

Place field detection program for WH repetition project
Algorithm written by Will Snider (rotation project)


Basic alg here for YZ 7/2021. Just the detection part, not the extra stuff
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: William Snider
"""

## import modules 
import numpy as np
from scipy.ndimage import label, center_of_mass, binary_dilation, generate_binary_structure
from scipy.ndimage.morphology import binary_fill_holes
from skimage.feature import peak_local_max
from skimage.morphology import watershed



# Parameters
# WH added params from ratterdam_defaults
ptsCm_macaulay = 4.72
ptsCm = ptsCm_macaulay

#old bins = 15,30
singleAlleyBins = [13,8]
wholeAlleyBins = [480/ptsCm,640/ptsCm]


## Binning values and size of arena
# Values are now hardcoded, commented out lines calculated more precisely for individual plots
binWidth = 10
xmax = 640 #np.max([np.max(unit.spikes[:,1]),np.max(unit.position[:,1])])
ymax = 480 #np.max([np.max(unit.spikes[:,2]),np.max(unit.position[:,2])])
xedges = np.arange(0, xmax+binWidth, binWidth) #np.arange(0, int(np.ceil(xmax/binWidth)*(binWidth+1)), binWidth)
yedges = np.arange(0, ymax+binWidth, binWidth) #np.arange(0, int(np.ceil(ymax/binWidth)*(binWidth+1)), binWidth)

## Gaussian Filter for 2D Rate map
sigma = 1.5

velocity_filter_thresh = 1 #TODO: play this with value #units are cm/s

cameraFrameRate = 30# fps. Used to give accurate spiking rate.

## Algorithm for detecting place fields
min_PF_proportion = 0.02 #(approximately 22 bins under testing setup)
prctile = 95 # 
cutoffFromPeakValue = 0.25 # A field is determined by its peak and the surrounding area including bins that are at least X% of the max firing rate
absFieldCutoff = 1.5# absolute FR below which an expanding field will be cutoff. cutoffFromPeakValue gives a realtive threshold of when to cutoff a field.
minPeakFiringRate = 1. # #Fields have to fire more than XHz at their peak
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



class RateMap():
    """
    Class that plots 2D ratemaps for ratterdam data analysis.
    """
    
    def __init__(self, position, ratemap):
        
        self.position = position
        self.rateMap2D = ratemap
        self.maxFiringRate = np.nanmax(self.rateMap2D)
        self.PF = []
        self.findPlaceFields()
        #self.findVisits()
        #self.PF_df = self.findPlaceFieldPatterns()

    
    def findPlaceFields(self):
        
        #######################################################################
        def checkField(self, PlaceFieldObject):
            
            """
            Function that tests if a place field meets certain requirements, like being above a min size and min peak firing rate. 
            """ 
            validField = True
            rateMap = self.rateMap2D 
            PF = PlaceFieldObject 
            
            # Peaks must be over a certain size
            PF_size = np.sum(PF.mask)
            binaryArena = np.where(np.isnan(rateMap), 0, 1)
            arenaSize = np.sum(binaryArena)
            PF_size_proportion = PF_size / arenaSize
            
            # Other tests
            _, numObjects = label(PF.mask)
            # Fields must be at least certain size
            if PF_size_proportion < min_PF_proportion:
                validField = False  # mark field object to be deleted
            
            # Peak Values must be over a certain firing rate
            elif PF.peakValue < minPeakFiringRate:
                validField = False  # mark field object to be deleted
            
            # Peaks must be at least a certain factor below the maxFiringRate of the rate map
            # Likely set to 0 but can check in williamDefaults
            elif PF.peakValue < minFiringRateFactor * self.maxFiringRate:
                validField = False  # mark field object to be deleted
            
            # Peaks must be contiguous (no split peaks)
            elif numObjects != 1:
                validField = False
                
            return validField

        def findUsefulInfo(self, PlaceFieldObject):
            """
            Function that finds useful info (e.g. center of mass) of a particular place field.
            """
            PF = PlaceFieldObject
            
            # Find field with values as opposed to the binary mask
            PF_array = np.where(PF.mask, rateMap, 0)
            field = PF_array
            PF.field = field
        
            # Find center of mass
            PF_com = center_of_mass(PF_array)
            com_list = (int(round(PF_com[0])), int(round(PF_com[1]))) #TODO: Remove this rounding as it loses data
            com = com_list
            PF.com = com

            # Find perimeter
            #ryl edit 1/2021: increase resolution of PF.mask before using binary_dilation to find the perimeter
            # increase resolution of PF.mask by inserting a column of zeros between every 2 columns in PF.mask and inserting a row of zeros between every 2 rows
            r, c = PF.mask.shape
            PFmaskExpanded = np.zeros((r*2,c*2), dtype=PF.mask.dtype)
            PFmaskExpanded[::2, ::2] = PF.mask
            
            # filling in the new pixels using values to the left and right
            for i in range(0,r*2,2):
                for j in range(1,c*2-1,2):
                    PFmaskExpanded[i,j] = PFmaskExpanded[i,j-1] and PFmaskExpanded[i,j+1]
                    
            # filling in the new pixels using values above and below
            for i in range(1,r*2-1,2):
                for j in range(c*2):
                    PFmaskExpanded[i,j] = PFmaskExpanded[i-1,j] and PFmaskExpanded[i+1,j]
                    
            # adding rows and columns of zeros at the ends of PFmaskExpanded
            PFmaskExpanded = np.vstack((np.zeros((1,PFmaskExpanded.shape[1])), PFmaskExpanded))
            PFmaskExpanded = np.hstack((np.zeros((PFmaskExpanded.shape[0],1)), PFmaskExpanded))
            
            PF_perimeter = np.where(binary_dilation(PFmaskExpanded, structure=generate_binary_structure(2, 2))!=binary_fill_holes(PFmaskExpanded)) #binary expanded - original
            #PF_perimeter = np.where(binary_dilation(PF.mask, structure=generate_binary_structure(2, 2))!=binary_fill_holes(PF.mask)) #binary expanded - original
            perimeter = PF_perimeter
            PF.perimeter = perimeter
            
            # Find peakValue
            peakValue = np.nanmax(PF.field)
            PF.peakValue = peakValue

            # Find peakXY
            #wh edit 2/2020: logic of finding peakXY to date is find loc of max
            # but if that isnt unique then you try to cast a tuple to int, not possible
            # chnage to a matrix based approach

            PF.peakXY = np.unravel_index(PF.field.argmax(), PF.field.shape)
            
            # Find Up/Down/Left/Right extrema, for analysis of navigation. 
            valsY, valsX = np.where(PF.mask)
            medianY = np.median(np.unique(valsY))
            minY = np.min(valsY)
            maxY = np.max(valsY)
            medianX = np.median(np.unique(valsX))
            minX = np.min(valsX)
            maxX = np.max(valsX)
            PF.navUp = [maxY, medianX]
            PF.navDown = [minY, medianX] 
            PF.navLeft = [medianY, minX] 
            PF.navRight = [medianY, maxX]
            
           
            
        #######################################################################

        # Algorithm for finding place fields:
        
        # Delete PFs if already found
        if any(self.PF[:]):
            print('Place Fields found previously. Deleted original and found again.')
            self.PF[:] = []
            
        # Find peaks (local extrema)
        rateMap = self.rateMap2D
        rateMap_no_nan = np.where(np.isnan(rateMap),0,rateMap)
        peaks = peak_local_max(rateMap_no_nan, indices=False) # indices -> returns boolean of local extrema
        peaks = np.where(peaks==True)
        
        # Store peak coordinates and firing rates in a 3xnumberofPFs array, and sort array in  descending max firing rate
        peakValuesList = rateMap[peaks].reshape([1,len(peaks[0])]) # column of peak firing values
        peaksList = np.concatenate((peaks,peakValuesList),axis=0)
        peaksList = peaksList[:,(-peaksList[2,:]).argsort()] #sorts according to 3rd column - firing rate, in descending order
        
        # Create peak objects to facilitate filtering further
        for idx, val in enumerate(peaksList[0]):
             self.PF.append(FieldObject())
             self.PF[idx].peakXY = (int(peaksList[0][idx]), int(peaksList[1][idx]))
             self.PF[idx].peakValue = rateMap[self.PF[idx].peakXY]
             
        
        # Expand peaks to include all pixels at least 20% of the max firing rate and touching
        # Crucially, the peaksList array above is sorted descendingly by peak firing rate. Therefore, when we remove avaiable bins for each PF, the next cannot contain any bins that the previous ones did. This ensures that there is no overlap of place fields.#PF with the highest firing rate wins out bins that are in competition 
        binsTaken = np.zeros(np.shape(self.rateMap2D), dtype='bool')
                
        for PF in self.PF[:]:
            
            # If bin of peak is already taken, delete the place field and continue
            if binsTaken[PF.peakXY] == True:
                self.PF.remove(PF)
                continue
            
            # Expand the peak out
            fieldThreshold = PF.peakValue * cutoffFromPeakValue
            aboveThreshold = np.where(rateMap >= max(fieldThreshold,0), True, False)
            peakMask = np.where(aboveThreshold, True, False) #np.where(aboveThreshold * belowPeak, True, False) # contains all points in rateMap above threshold; must choose single blob that contains peak
            peakMask = np.where(binsTaken, False, peakMask) #remove peaks that are already occupied
            labeledBlobs, numberOfBlobs = label(peakMask)
            
            blobList = range(1, numberOfBlobs+1)    
            for blob in blobList:
                blobMask = np.where(labeledBlobs==blob,True,False) #convert blob to binary     
                if blobMask[PF.peakXY] == False: 
                    labeledBlobs = np.where(blobMask, 0, labeledBlobs) # remove all blobs that do not contain the peak
            PF.mask = np.where(labeledBlobs,True,False) # contains single blob that has the peak
            binsTaken = np.where(PF.mask, True, binsTaken) #mark occupied bins
        
        
        # Find useful information about the peaks using findUsefulInfo and store
        for PF in self.PF[:]:
            findUsefulInfo(self,PF)
                
        # Filter peaks based on parameters using checkField 
        for PF in self.PF[:]: 
            validField = checkField(self, PF)
            if validField is False:
                self.PF.remove(PF) #delete the field object
                
        # Split fields that contain multiple peaks into smaller subfields
        for PF in self.PF[:]:
            field = PF.field
            mask = PF.mask
            peaks = peak_local_max(field, min_distance=minPeakDistanceAwayBins, exclude_border=False, indices=False) # peaks must be minPeakDistanceAwayBins # of bins away from another peak
            peaksLabeled, numberOfPeaks = label(peaks)
            
            # Skip place fields with only 1 local max
            if numberOfPeaks == 1:
                continue
             
            newFieldsOkay = False
            while newFieldsOkay is False:
                
            
                # Split field based on assigned peaks
                labels = watershed(-field, peaksLabeled, mask=mask)
                
                # Make list of new field labels. Remove 0 if it is there.
                newFieldList = np.unique(labels)[np.unique(labels)!=0]
                numberOfPeaks = len(newFieldList)
                
                # Store as FieldObjects to facilitate filtering
                newFields = []
                for fieldNum in newFieldList:
                    newFields.append(FieldObject())
                    newFields[-1].mask = np.where(labels==fieldNum,True,False)
                    findUsefulInfo(self, newFields[-1])
                                      
                # Test if each new field is okay according to checkField
                for NF in newFields:
                    validField = checkField(self, NF)
                    if validField is False:
                        peaksLabeled[NF.mask] = 0; # removes issue of discrepancy between peak and local max for a subfield
                        break
                    
                # If a field was invalid, reset the while loop to recalculate
                if validField is False:
                    continue
                
                # Recombine fields if the difference between them is small
                # Find border values of each newField, store in newField[-1].borderXY
                labelsWithBorder = watershed(-field, peaksLabeled, mask=mask, watershed_line=True)
                border = np.where(labels!=labelsWithBorder,True,False)
                borderLabeled, _ = label(border,  structure=generate_binary_structure(2, 2))
                borderList = np.unique(borderLabeled[borderLabeled>0])
                newFieldsSorted = []
                for NF in newFields:
                    NF.borderXY = [] #incase it already has contents
                    for bordNum in borderList:
                        bordXY = np.where(borderLabeled==bordNum)
                        if any(NF.mask[bordXY]):
                            NF.borderXY.append(bordXY)
                        elif any(binary_dilation(NF.mask[bordXY])):
                            NF.borderXY.append(bordXY)
                    
                    # Find number of borders for each new field
                    NF.numberOfBorders = len(NF.borderXY)
                    newFieldsSorted = newFieldsSorted + [NF] + [NF.numberOfBorders]
                
                # sort fields based on number
                newFieldsSorted = np.reshape(newFieldsSorted, [len(newFields),2])
                newFieldsSorted = newFieldsSorted[newFieldsSorted[:,1].argsort()] #sort by number of borders found

                # Find mean of points at border and compare if they're below threshold, start with new fields with only single border
                # Crucially, this starts with new fields with only one border. IF there is an issue, the program recalculates after dropping that new fields peak.
                for NF in newFieldsSorted[:,0]:
                    for bord in NF.borderXY: #some fields have two borders
                    # Find mean of border
                        avgBorderValue = np.nanmax(field[bord]) #mean here gets too weighetd down by edge values of border
                    # Compare mean to peak of new field
                        if avgBorderValue > minValleyProportion * NF.peakValue:
                            validField = False
                            peaksLabeled[NF.mask] = 0; # removes issue of discrepancy between peak and local max for a subfield
                            break
                        
                    if validField is False: # break to outer for loop
                        break
                    
                # If a field was invalid, reset the while loop to recalculate
                if validField is False:
                    #print(str(PF)+' had a field be recombined.')
                    continue
                
                # All fields passed the test, so escape while loop
                newFieldsOkay = True
            # Add newFields to current PF list and remove current PF that it was split from
            for newField in newFields:
                self.PF.append(newField) # add in the new sub fields
            self.PF.remove(PF)
 

        # Test if vertical, horizontal, or intersection    #patternDict = #dict contains: PF, location #PFdf[]

class FieldObject():
    
    def __init__(self):
        self.mask = None
        self.field = None
        self.com = None
        self.perimeter = None
        self.PF_num = None
        self.peakXY = None
        self.peakValue = None
        self.borderXY = [] # used to determine if fields should be split
        self.numberOfBorders = None
        self.visits = []
        self.allVisitsFR = None
        self.allVisitsRange = None
        self.allVisitsMidpoint = None
        self.navUp = None
        self.navDown = None
        self.navLeft = None
        self.navRight = None
        self.alleyInter = None
        self.df = None
        
class VisitObject():
    
    def __init__(self):
        self.path = None
        self.num = None
        self.pathExt = None # a few seconds before and after for plotting
        self.meanFR = None
        self.timeRange = None
        self.midpoint = None
        self.entry = None
        self.exit = None
        self.direction = None
        self.AISeries = None
        self.skelAISeries = None
        self.skelStartFrame = None
        self.skelEndFrame = None
        
