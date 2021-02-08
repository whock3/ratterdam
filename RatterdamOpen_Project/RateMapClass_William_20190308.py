#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: William Snider
"""

## import modules 
import numpy as np
import matplotlib.pyplot as plt
from williamDefaults import xedges, yedges, sigma, cameraFrameRate, min_PF_proportion, cutoffFromPeakValue, minPeakFiringRate, minFiringRateFactor, binWidth, minPeakDistanceAwayBins, minValleyProportion, timeStampScale, minPFOccupancyTime, minGapLength, alleyInterBounds, alleyInterType, skelPts, runningMeanN,minVisitSizeAI, absFieldCutoff
# import ratterdam_CoreDataStructures as core
from astropy.convolution import convolve, Gaussian2DKernel
from scipy.ndimage import label, center_of_mass, binary_dilation, generate_binary_structure
from scipy.ndimage.morphology import binary_fill_holes
from skimage.feature import peak_local_max
from skimage.morphology import watershed
import pdb
import pandas as pd
from scipy.spatial import KDTree
from collections import Counter
from itertools import groupby

## define path

## make class
class RateMap():
    """
    Class that plots 2D ratemaps for ratterdam data analysis.
    """
    
    def __init__(self, unit):
        
        self.unit = unit
        self.rateMap2D = self.makeRateMap2D()
        self.maxFiringRate = np.nanmax(self.rateMap2D)
        self.PF = []
        self.findPlaceFields()
        #self.findVisits()
        #self.PF_df = self.findPlaceFieldPatterns()

    
    def makeRateMap2D(self):
        
        # Make 2D histogram of spiking
        spikeHist, _ , _ = np.histogram2d(self.unit.spikes[:,1],self.unit.spikes[:,2], bins = (xedges, yedges))
        spikeHist = spikeHist.T
        
        # Make 2D histogram of occupancy
        occHist, _ , _ = np.histogram2d(self.unit.position[:,1], self.unit.position[:,2], bins = (xedges, yedges))
        occHist = occHist.T
        
        # Normalize firing rate by dividing spikeHist by occHist, suppress divide by zero warning
        with np.errstate(invalid='ignore'):
            spikeHistNorm = np.divide(spikeHist, occHist);
            
        #Adjust firing rate for fps of camera
        spikeHistNorm = spikeHistNorm*cameraFrameRate
        
        # Gaussian filter the data to smooth
        kernel = Gaussian2DKernel(sigma)
        spikeHistConv = convolve(spikeHistNorm, kernel, nan_treatment='fill', fill_value=0, preserve_nan=True)
        return spikeHistConv
    
    def plotRateMap2D(self):
        
        fig, ax = plt.subplots()
        ax.imshow(self.rateMap2D, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])                                      
        
    def makeNiceRateMap2D(self,ax):
        
        maxFiringRateRounded = np.str(round(self.maxFiringRate,2))
        rateMapTitle = ((self.unit.name + '   Peak FR (spikes/s): ' + maxFiringRateRounded))
        ax.set_title(rateMapTitle, fontsize = 3)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    
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
            
            # Find alley/intersection in which a place field's com is located
            for alleyInter in alleyInterBounds:
                [[x1, x2], [y1, y2]] = np.divide(alleyInterBounds[alleyInter], binWidth)
                (PFy, PFx) = PF.com
                if PFx >= x1 and PFx < x2 and PFy >= y1 and PFy < y2:
                    PF.alleyInter = alleyInter
                    break
            
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
 
    def findVisits(self, normalize=True):
        
        """
        For each place field, return the start time and mean firing rate for each visit.
        """
        #######################################################################
        def countSpikesWithinVisit(allSpikeTimes, visitTransitionsTimes):
            """
            This function accounts for gaps within a visit and does not count spikes occuring during a gap.
            """
            numSpikes = 0                
            # loop over just the subvisits and count spikes
            for idx,_ in enumerate(visitTransitionsTimes[1::2]):

                subvisitStart = visitTransitionsTimes[idx]
                subvisitEnd = visitTransitionsTimes[idx-1]
                
                spikesBefore = [allSpikeTimes >= subvisitStart]
                spikesAfter = [allSpikeTimes <  subvisitEnd]
                spikesWithinSubvisit = np.logical_and(spikesBefore,spikesAfter)
                numSpikes = numSpikes + np.sum(spikesWithinSubvisit)
            
            return numSpikes
        
        def findAI(positionsInBins, minVisitSizeAI):
            
            """Determine the alley or intersection for every position of the rat, as well as the preceding and succeeding alleys or intersections"""
            
            # Make array with Alley/Intersection graphed ('lookup table')
            AIBounds2D = np.empty([65,50],dtype='<U2')
            for AI in alleyInterBounds.keys():
                [[x1, x2], [y1, y2]] = np.floor(np.divide(alleyInterBounds[AI], binWidth)).astype(int)
                AIBounds2D[x1:x2,y1:y2] = AI
            
            allAI = np.empty([len(positionsInBins),1],dtype='<U2')
            allAI[:] = '-'
            for idx,pos in enumerate(positionsInBins):
                pos = tuple(pos[::-1])
                allAI[idx] = AIBounds2D[pos]
        
            # Condense allAI into non-repeating list of alley,interesction visits
            groupedAIVisits = groupby(allAI)
            AIVisitOrder = [[''.join(labelAI.tolist()), sum(1 for _ in group)] for labelAI, group in groupedAIVisits]

            # Remove short AI visits (and non-visits); expand to bigger neighboring visit
            minVisitAI = np.round(minVisitSizeAI*cameraFrameRate).astype(int)
            
            # Replace short visits with ''
            for i, [_ , count] in enumerate(AIVisitOrder):
                if count < minVisitAI:
                    AIVisitOrder[i][0] = ''
            
            # Expand and regroup
            expanded = [[AI]*int(count) for AI, count in np.array(AIVisitOrder)]
            expanded = [y for x in expanded for y in x] #flattens list
            groupedAIVisits = groupby(expanded)
            AIVisitOrder = [[labelAI, sum(1 for _ in group)] for labelAI, group in groupedAIVisits]

            # Assign largest neighbor to '' values. 
            for i, [alleyORintersection, count] in enumerate(AIVisitOrder):
                if alleyORintersection == '': #implies that neighbors cannot be '' themselves
                    
                    if i > 0: # handle edge case
                        prevCount = AIVisitOrder[i-1][1]
                        prevValue = AIVisitOrder[i-1][0]
                    else:
                        prevCount = 0
                        prevValue = ''
                    
                    # Find count of next AI, handle first element edge case 
                    try:
                        nextCount = AIVisitOrder[i+1][1]
                        nextValue = AIVisitOrder[i+1][0]
                    except:
                        nextCount = 0
                        nextValue = ''
                    
                    if prevCount >= nextCount:
                        AIVisitOrder[i][0] = prevValue
                    elif prevCount < nextCount:
                        AIVisitOrder[i][0] = nextValue
        
            # Expand and regroup again
            expanded = [[AI]*int(count) for AI, count in np.array(AIVisitOrder)]
            expanded = [y for x in expanded for y in x] #flattens list
            groupedAIVisits = groupby(expanded)
            AIVisitOrder = [[labelAI, sum(1 for _ in group)] for labelAI, group in groupedAIVisits]       
                        
            # Calculate previous and next alleys/intersections and add to df
            prevPrevAI = np.empty(np.shape(allAI),dtype='<U2')
            prevPrevAI[:] = '-'
            prevAI = np.empty(np.shape(allAI),dtype='<U2')
            prevAI[:] = '-'
            currAI = np.empty(np.shape(allAI),dtype='<U2')
            currAI[:] = '-'
            nextAI = np.empty(np.shape(allAI),dtype='<U2')
            nextAI[:] = '-'
            nextNextAI = np.empty(np.shape(allAI),dtype='<U2')
            nextNextAI[:] = '-'

            IDX = 0
            for i, [alleyORintersection, count] in enumerate(AIVisitOrder):
                                
                # Previous-Previous alley or intersection
                if i >= 2: 
                    prevPrevAI[IDX:IDX+count] = AIVisitOrder[i-2][0]

                # Previous alley or intersection
                if i >= 1: 
                    prevAI[IDX:IDX+count] = AIVisitOrder[i-1][0]

                # Current alley or intersection
                currAI[IDX:IDX+count] = AIVisitOrder[i][0]

                # Next alley or intersection
                if i <= len(AIVisitOrder)-2: 
                    nextAI[IDX:IDX+count] = AIVisitOrder[i+1][0]

                # Next Next alley or intersection
                if i <= len(AIVisitOrder)-3: 
                    nextNextAI[IDX:IDX+count] = AIVisitOrder[i+2][0]

                IDX = IDX + count
            
            # indices of starting frames for each AI
            prevPrevAIidxS = findTransitions(prevPrevAI, -2) 
            prevAIidxS = findTransitions(prevAI, -1)
            currAIidxS = findTransitions(currAI, 0) 
            nextAIidxS = findTransitions(nextAI, 1 ) 
            nextNextAIidxS = findTransitions(nextNextAI, 2) 
            
            # indices of ending frames for each AI
            prevPrevAIidxE = findTransitions(prevPrevAI, -1) 
            prevAIidxE = findTransitions(prevAI, 0)
            currAIidxE = findTransitions(currAI, 1) 
            nextAIidxE = findTransitions(nextAI, 2) 
            nextNextAIidxE = findTransitions(nextNextAI, 3) 
            
            # Concatenate
            AISequences = np.concatenate((prevPrevAI,prevAI,currAI,nextAI,nextNextAI), axis=1).tolist()           
            AIidxS = np.concatenate((prevPrevAIidxS, prevAIidxS, currAIidxS, nextAIidxS, nextNextAIidxS), axis=1).tolist()
            AIidxE = np.concatenate((prevPrevAIidxE, prevAIidxE, currAIidxE, nextAIidxE, nextNextAIidxE), axis=1).tolist()
            return AISequences, AIidxS, AIidxE
        
        def findTransitions(arr, adjust):
            """ Return the index of the alley or intersection as opposed to identity
                adjust is how much the transitions should shift. For example, 
                nextAI sequence should be +1. """
            
            transitions = [0,len(arr)]
            
            # Find transitions
            preVal = arr[0]
            for i, curVal in enumerate(arr):
                
                if curVal != preVal:
                    transitions.insert(-1,i)                    
                preVal = curVal
              
            # Find transition frames
            pret = 0
            transitionFrames = np.zeros(np.shape(arr),dtype='int32')
            for i, currt in enumerate(transitions):
                 
                
                if currt == 0: 
                    continue
                
                if adjust == -2 and i <= 2:
                    transitionFrames[pret:currt] = 0 #transitions[i-1+adjust]
                elif adjust == -1 and i <= 1:
                    transitionFrames[pret:currt] = 0 #transitions[i-1+adjust]
                elif adjust == 1 and i >= len(transitions) - 1:
                    transitionFrames[pret:currt] = transitions[-1]
                elif adjust == 2 and i >= len(transitions) - 2:
                    transitionFrames[pret:currt] = transitions[-1]
                elif adjust == 3 and i >= len(transitions) - 3:
                    transitionFrames[pret:currt] = transitions[-1]
                else:
                    transitionFrames[pret:currt] = transitions[i-1+adjust]
                pret = currt
            
            return transitionFrames
        
        def running_mean(y, N): # Copied from stack exchange, handles edges (not zero padded)
            y_padded = np.pad(y, (N//2, N-1-N//2), mode='edge')
            y_smooth = np.convolve(y_padded, np.ones((N,))/N, mode='valid') 
            return y_smooth
        #######################################################################
        
        # %% Find timestamps within this place field
        positions = self.unit.position[:,[2,1]]
        positionsInBins = np.floor(self.unit.position[:, [2, 1]]/binWidth).astype(int)
        allCameraTimes = self.unit.position[:,0]/timeStampScale
        allSpikeTimes = self.unit.spikes[:,0]/timeStampScale

        for PF in self.PF:            
            #Make vector of position values inside a place field
            inPF = np.zeros(len(positionsInBins),dtype='bool')
            for idx, pos in enumerate(positionsInBins):
                inPF[idx] = PF.mask[tuple(pos)] 
            
            # label nonVisits
            outPF = ~inPF
            gaps, _ = label(outPF)
            gapList = np.unique(gaps)[np.unique(gaps)!=0]
    
            for g in gapList:
                
                # Calculate gap length
                gapFrames = np.where(gaps == g)# test if shorter than minGapLength
                gapFrames = np.concatenate((gapFrames[0],np.array([gapFrames[0][-1]+1]))) # concatenate additional frame, necessary to get correct endtime
                gapStartTime = allCameraTimes[gapFrames[0]]
                try:
                    gapEndTime = allCameraTimes[gapFrames[-1]]
                except:
                    #print('Indexing error exception caught. This may  happen if the rat was not in a place field during the last recorded camera frame. Visit number: '+str(g))
                    gapEndTime = allCameraTimes[gapFrames[-2]]+1/cameraFrameRate   
                gapLength = gapEndTime - gapStartTime
                
                # Remove gap if length is too short
                if gapLength < minGapLength:
                    #print('Gap #'+str(g)+' was removed because its length was only: '+str(gapLength)+' s')
                    outPF[gapFrames[0:-1]] = False
            
            # Label visits to account for short gaps
            newInPF = ~outPF  # newInPF takes into account visits separated by only a very short gap
            visits, _ = label(newInPF)
            #visits = np.where(inPF, visits, 0)# Important! Remove all visit timepoints that were actually gaps 
            visitList = np.unique(visits)[np.unique(visits)!=0]
    
            # Create vector of transitions. [a, a, a, b, b, c] should be [1,0,0,1,0,1]
            transitions = np.zeros(np.shape(inPF),dtype='bool')
            transitions[0] = True # mark first frame as transition point
            transitions[-1] = True # mark last frame as transition point       
            for idx, _ in enumerate(transitions[0:-2]):
                t1 = inPF[idx]
                t2 = inPF[idx+1]
                
                if t1 != t2:
                    transitions[idx+1] = True
            transitionsFrames = np.where(transitions)
            
            # %% Find AI sequences
            AISequences, AIidxS, AIidxE = findAI(positionsInBins, minVisitSizeAI)
            
            # %% Find closest skeleton point of each position
            T = KDTree(skelPts)
            neighbor_dists, neighbor_indices = T.query(positions)
            skelPos = skelPts[neighbor_indices]
            
            # Smooth these points to avoid large jumps/ensure intersection is always visited          
            xSmooth = running_mean(skelPos[:,0], runningMeanN)
            ySmooth = running_mean(skelPos[:,1], runningMeanN)
            smoothSkelPos = np.array([xSmooth, ySmooth])
            
            # Again find closest skeleton point of the smoothed positions
            smoothSkelPos = [list(smoothSkelPos[:,x]) for x in np.arange(0,len(smoothSkelPos[0]))]
            neighbor_dists, neighbor_indices = T.query(smoothSkelPos)
            smoothSkelPos = skelPts[neighbor_indices]
            
            # Convert to bin coordinates
            smoothSkelPos = np.floor(np.divide(smoothSkelPos, binWidth)).astype(int)
            
            # Find smoothed AI sequences
            skelAISequences, skelAIidxS, skelAIidxE = findAI(smoothSkelPos, minVisitSizeAI)
            
            # %% Add time-relevant values to dataframe
            data = {'allCameraTimes': allCameraTimes, 'positionsInBinsX': positionsInBins[:,0],'positionsInBinsY': positionsInBins[:,1],'AISequences': AISequences, 'AIidxS':AIidxS,'AIidxE':AIidxE, 'skeletonPoints': smoothSkelPos.tolist(), 'skelAISequences' :skelAISequences, 'skelAIidxS': skelAIidxS,'skelAIidxE': skelAIidxE, 'inPF': inPF,'newInPF':newInPF, 'outPF': outPF, 'transitions': transitions}
            df = pd.DataFrame(data)
            PF.df = df
            
            # %% Find each visit
            
            # Vectors to store info relative to all 
            allVisitsFR = []
            allVisitsRange = []
            allVisitsMidpoint = []
            experimentStartTime = allCameraTimes[0]
            visitNumber = 0
            
            for vis in visitList:
                
                # Convoluted but more accurate than assuming each frame is exactly 1/cameraFrameRate
                # WithGaps means it contains frames that should be excluded because it is a gap
                visitFramesWithGaps = np.where(visits==vis)
                visitFramesWithGaps = np.concatenate((visitFramesWithGaps[0],np.array([visitFramesWithGaps[0][-1]+1]))) #actually one frame longer than the visit, necessary to get correct visitLength
                visitFramesLessGaps = np.intersect1d(np.where(inPF), visitFramesWithGaps)
                visitTransitionFrames = np.intersect1d(visitFramesWithGaps,transitionsFrames)
                visitTransitionsTimes = allCameraTimes[visitTransitionFrames]
                visitStartTime = allCameraTimes[visitFramesWithGaps[0]]
                try:
                    visitEndTime = allCameraTimes[visitFramesWithGaps[-1]]
                except:
                    #print('Indexing error exception caught. This should only happen if the rat was in a place field during the last recorded camera frame. Visit number: '+str(vis))
                    visitEndTime = allCameraTimes[visitFramesWithGaps[-2]]+1/cameraFrameRate   
                # Sum up every second timestamp--ending with final departure from PF and subtract the sum of every entrance to the PF. This simply subtracts out the gaps
                visitLength = np.sum(allCameraTimes[visitTransitionFrames[1::2]])-np.sum(allCameraTimes[visitTransitionFrames[0::2]])            
                
                # Same values as above but with time where 0 = start of experiment
                visitStartTimeNorm = round(visitStartTime - experimentStartTime,2)
                visitEndTimeNorm = round(visitEndTime - experimentStartTime,2)
                visitMidpointNorm = round((visitStartTime + visitLength/2)-experimentStartTime,2) #centered on visit and time starting from experiment start
                visitRange = (visitStartTimeNorm, visitEndTimeNorm) # Gap frames included here for visualization
    
                # Only include visits longer than min occupancy time
                if visitLength < minPFOccupancyTime:
                    continue
                
                # Count number of spikes in this time
                numSpikes = countSpikesWithinVisit(allSpikeTimes,visitTransitionsTimes)
                firingRate = numSpikes/visitLength 
                
                # If desired, normalize to spatial rate map
                if normalize is True:
                    visitBins = positionsInBins[visitFramesLessGaps]
                    visitBinsX, visitBinsY = zip(*visitBins) # unpack tuple to give two lists of x and y coordinates, which is necessary to index numpy array
                    expectedFR = np.mean(PF.field[visitBinsX,visitBinsY])
                    firingRate = firingRate/expectedFR # divide to normalize
        
                # Append allvisit info to PF
                allVisitsFR.append(firingRate)
                allVisitsRange.append(visitRange)
                allVisitsMidpoint.append(visitMidpointNorm)

                # Find positions of visit, also before and after in seconds
                visitPositions = positions[visitFramesWithGaps[:-1]] #last frame was added, not in PF
                timeBeforeAfter = 3 #units = s
                numFramesBA = timeBeforeAfter*cameraFrameRate
                visitPositionsBA = positions[int(visitFramesWithGaps[0]-numFramesBA):int(visitFramesWithGaps[-1]+numFramesBA+1)]
            
                # Create visit object
                PF.visits.append(VisitObject())

                # %% Find movement through PF direction 
                
                #Find which end of the PF the visit started 
                posEnter = positionsInBins[visitFramesWithGaps[0]]
                ptsEnter = np.array([posEnter,PF.navUp, PF.navDown, PF.navLeft, PF.navRight])
                T = KDTree(ptsEnter)
                nearest_dist, nearest_ind = T.query(ptsEnter, k=2)
                entry = nearest_ind[0,1]
                if entry == 1: PF.visits[-1].entry = 'Up'
                elif entry == 2: PF.visits[-1].entry = 'Down'
                elif entry == 3: PF.visits[-1].entry = 'Left'
                elif entry == 4: PF.visits[-1].entry = 'Right'
                else:
                    'Error in the visit entrypoint code'
                    
                # Find where the PF visit ended   
                posExit = positionsInBins[visitFramesWithGaps[-2]] # -1 was added, -2 is the last frame WITHIN the visit
                ptsExit = np.array([posExit,PF.navUp, PF.navDown, PF.navLeft, PF.navRight])
                T = KDTree(ptsExit)
                nearest_dist, nearest_ind = T.query(ptsExit, k=2)
                ext = nearest_ind[0,1]
                if ext == 1: PF.visits[-1].exit = 'Up'
                elif ext == 2: PF.visits[-1].exit = 'Down'
                elif ext == 3: PF.visits[-1].exit = 'Left'
                elif ext == 4: PF.visits[-1].exit = 'Right'
                else:
                    'Error in the visit exitpoint code'
                    
                # %% Find direction of turns before and after visit to a place field
                visitMidpoint = visitFramesWithGaps[0] # first point of entry #visitFramesWithGaps[int(len((visitFramesWithGaps)-1)/2)]
                PF.visits[-1].AISeries = df['AISequences'].iloc[visitMidpoint]
                PF.visits[-1].skelAISeries = df['skelAISequences'].iloc[visitMidpoint]
                                
                
                # %% Append visit info
                PF.visits[-1].path = visitPositions
                PF.visits[-1].pathExt = visitPositionsBA
                PF.visits[-1].num = visitNumber
                PF.visits[-1].meanFR = firingRate
                PF.visits[-1].timeRange = visitRange
                PF.visits[-1].midpoint = visitMidpointNorm
                PF.visits[-1].direction = [PF.visits[-1].entry+PF.visits[-1].exit]
                PF.visits[-1].skelStartFrame = PF.df['skelAIidxS'][visitFramesWithGaps[0]][2]
                PF.visits[-1].skelEndFrame = PF.df['skelAIidxE'][visitFramesWithGaps[0]][4]
    
                # Increment visitNumber for next visit
                visitNumber += 1
              
            # Assign all visit info to PF object   
            PF.allVisitsFR = allVisitsFR
            PF.allVisitsRange = allVisitsRange
            PF.allVisitsMidpoint = allVisitsMidpoint

    
    def findPlaceFieldPatterns(self):
        
        # Make dataframe with PFs and their location
        if len(self.PF) == 0: # removes error for no-PF units
            PF_df = pd.DataFrame([],columns=['PF','alleyInter','orientation'])
        else:
            alleyInterList = []
            for PF in self.PF:
                objType = ''.join(alleyInterType[PF.alleyInter])
                alleyInterList.append([PF, PF.alleyInter, objType])
            alleyInterList = np.array(alleyInterList)
            PF_df = pd.DataFrame(alleyInterList,columns=['PF','alleyInter','orientation'])
            PF_df.sort_values(by=['orientation'])
        
        return PF_df
#        # Count occurences of each orientation type
#        orientationCounts = Counter(PF_df['orientation'])
#        
#        # Make list of orientations and PFs corresponding, if over certain number
#        patternList = []
#        for orientations in orientationCounts:
#            
#            if orientationCounts[orientation] > 2:
#                
        

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
        
