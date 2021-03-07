# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 13:16:25 2021

@author: Ruo-Yah Lai
"""

import numpy as np
from random import randint, choice
from matplotlib import path
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation, generate_binary_structure, center_of_mass
from math import ceil
from williamDefaults import xmax, ymax, binWidth, xedges, yedges
from placeFieldBorders import reorderBorder
from newAlleyBounds import R859


def shuffleField(unit, subfield, rat, title="", graph=False):
    """
    Moves a field to a new location and assign firing rates to the new pixels
    using the firing rates in the original pixels
    
    rat: named tuple from newAlleyBounds
    """
    x = int(xmax/binWidth)
    y = int(ymax/binWidth)
    
    #find the field's area
    fieldPerim = path.Path(unit.perimeters[subfield])
    inField = fieldPerim.contains_points(coords)
    area = np.sum(inField)
    
    #find pixels with sampling, excluding pixels outside of the track
    binaryArena = np.where(np.isnan(unit.repUnit.rateMap2D), False, True) #rows=y, cols=x
    v1, h1 = rat.alleyLines[0,0]/binWidth
    v8, h6 = rat.alleyLines[-1,1]/binWidth
    onTrack = np.ones((y,x), dtype=bool)
    for row in range(0, ceil(h1)):
        onTrack[row,:] = False
    for row in range(ceil(h6), y):
        onTrack[row,:] = False
    for col in range(0, ceil(v1)):
        onTrack[:,col] = False
    for col in range(ceil(v8), x):
        onTrack[:,col] = False
    binaryArena = np.where(onTrack, binaryArena, False)
    
    #find new center of mass
    newCom = [0,0]
    while ~binaryArena[newCom[1], newCom[0]]:
        newCom = [randint(0,x-1), randint(0,y-1)] #[x,y]

    #dilate and mask until new area >= original area
    newFieldDilMask = np.zeros((y,x), dtype=bool) #rows=y, cols=x
    newFieldDilMask[newCom[1], newCom[0]] = True
    while np.sum(newFieldDilMask) < area:
        newFieldMask = newFieldDilMask
        dilated = binary_dilation(newFieldMask, structure=generate_binary_structure(2, 2)) #dilate in all directions, including diagonals
        newFieldDilMask = np.where(binaryArena, dilated, False)
    
    #remove pixels from the last dilation until new area = original area
    while np.sum(newFieldDilMask) > area:
        lastDil = np.where(newFieldDilMask != newFieldMask)
        lastDil = list(zip(lastDil[0], lastDil[1]))
        removePx = choice(lastDil)
        newFieldDilMask[removePx[0], removePx[1]] = False
    
    #find and sort firing rates in the original field
    field = np.where(unit.repUnit.PF[subfield].mask, unit.repUnit.rateMap2D, np.nan)
    fieldNonNaN = field[~np.isnan(field)]
    fieldRates = sorted(fieldNonNaN, reverse=True)
    
    #assign highest firing rate to the center of the new field
    newFieldDilMask2 = np.zeros((y,x), dtype=bool) #rows=y, cols=x
    newField = np.where(newFieldDilMask2, 0, np.nan) #with values, not just True or False
    newCenter = center_of_mass(newFieldDilMask)
    newCenter = (int(round(newCenter[0])), int(round(newCenter[1]))) #(y,x)
    newField[newCenter[0], newCenter[1]] = fieldRates[0]
    newFieldDilMask2[newCenter[0], newCenter[1]] = True
    fieldRates.pop(0)
    
    #assign firing rates to the pixels in the new field
    while np.sum(newFieldDilMask2) < area:
        #dilate from the pixels of the field with values and remove pixels outside of the field
        newFieldMask2 = newFieldDilMask2
        dilated = binary_dilation(newFieldMask2, structure=generate_binary_structure(2, 2))
        newFieldDilMask2 = np.where(newFieldDilMask, dilated, False)
        lastDil = np.where(newFieldDilMask2 != newFieldMask2)
        lastDil = list(zip(lastDil[0], lastDil[1]))
        
        #randomly assign rates to the pixels in the last dilation
        while len(lastDil) != 0:
            randomPx = choice(lastDil)
            lastDil.remove(randomPx)
            newField[randomPx[0], randomPx[1]] = fieldRates[0]
            fieldRates.pop(0)
        
    if graph:
        fig, ax = plt.subplots()
        ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                  cmap="jet", vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        ax.plot(unit.perimeters[subfield][:,0], unit.perimeters[subfield][:,1],color="k", label="original")
        ax.imshow(newField, origin='lower', aspect='auto', interpolation='None', 
                  cmap="jet", vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        newPerim = findPerim(newFieldDilMask)
        ax.plot(newPerim[:,0], newPerim[:,1], color="r", label="new")
        ax.set_title(title)
        ax.legend()
        ax.axis("equal")    
    

def findPerim(PFmask):
    """
    Finds the perimeter of a field given the boolean mask of the field
    Taken from RateMapClass_William
    """
    # increase resolution of PF.mask by inserting a column of zeros between every 2 columns in PF.mask and inserting a row of zeros between every 2 rows
    r, c = PFmask.shape
    PFmaskExpanded = np.zeros((r*2,c*2), dtype=PFmask.dtype)
    PFmaskExpanded[::2, ::2] = PFmask
            
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
            
    PF_perimeter = np.where(binary_dilation(PFmaskExpanded, structure=generate_binary_structure(2, 2))!=PFmaskExpanded) #binary expanded - original
    
    return reorderBorder(PF_perimeter,0)


Y, X = np.mgrid[0.5:ymax/binWidth+0.5, 0.5:xmax/binWidth+0.5]*binWidth
coords = np.array(list(zip(X.flatten(), Y.flatten())))