# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 19:31:38 2020

@author: Ruo-Yah Lai
"""
from matplotlib import path
import numpy as np
from RateMap import makeRM
from collections import namedtuple
from williamDefaults import alleyInterBounds, alleyInterType
import csv
from string import ascii_uppercase


def overlap(perim, alley):
    field = path.Path(perim)
    inField = field.contains_points(coords)
    alley = path.Path(alley)
    inAlley = alley.contains_points(coords)
    overlap = np.logical_and(inField, inAlley)
    return np.sum(overlap)


def mass(perim, alley, spikes, pos):
    """
    Finds the fraction of firing in an alley/intersection
    """
    hist = makeRM(spikes, pos, smoothing_2d_sigma=0)
    field = path.Path(perim)
    inField = field.contains_points(coordsBin)
    totalMass = np.nansum(hist.flatten()[inField])
    alley = path.Path(alley)
    inAlley = alley.contains_points(coordsBin)
    massInAlley = np.nansum(hist.flatten()[np.logical_and(inField,inAlley)]) #in the field and the alley/intersection
    return massInAlley / totalMass


def repeatingPF(perims, spikes, pos):
    """
    perims: subfields' perimeters
    """
    subfields = [[] for _ in range(len(perims))]
    overlaps = [[] for _ in range(len(perims))]
    allSubfields = []
    for i,perim in enumerate(perims):
        for j in range(17):
            alley = alleyInterBounds[str(j)]
            aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
            alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                   [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
            alleySize = (aRect.xmax-aRect.xmin) * (aRect.ymax-aRect.ymin)
            if overlap(perim, alleyPerim) > alleySize*0.75:
                subfields[i].append(alleyInterType[str(j)][1])
                overlaps[i].append(j)
        for j in ascii_uppercase[:12]:
            intersection = alleyInterBounds[j]
            iRect = Rectangle(intersection[0][0], intersection[1][0], 
                              intersection[0][1], intersection[1][1])
            iPerim = np.array([[iRect.xmin, iRect.ymin], [iRect.xmax, iRect.ymin],
                               [iRect.xmax, iRect.ymax], [iRect.xmin, iRect.ymax]])
            Size = (iRect.xmax-iRect.xmin) * (iRect.ymax-iRect.ymin)
            if overlap(perim, iPerim) > Size*0.75:
                subfields[i].append(alleyInterType[j][0])
                overlaps[i].append(j)
        
        print(subfields[i])
        if len(set(subfields[i])) > 1: #if the subfield is in multiple location types
            masses = np.empty(0) #mass in each location
            for k, j in enumerate(overlaps[i]):
                alley = alleyInterBounds[str(j)]
                aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
                alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                       [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
                masses = np.hstack((masses, mass(perim, alleyPerim, spikes, pos)))
            horiz = np.sum(masses[np.where(np.array(subfields[i]) == "horizontal")])
            vert = np.sum(masses[np.where(np.array(subfields[i]) == "vertical")])
            intxn = np.sum(masses[np.where(np.array(subfields[i]) == "intersection")])
            print(horiz,vert,intxn)
        
            subfields[i] = []
            if horiz > 0.25:
                subfields[i].append("horizontal")
            if vert > 0.25:
                subfields[i].append("vertical")
            if intxn > 0.25:
                subfields[i].append("intersection")
        
        for j in subfields[i]:
            allSubfields.append(j)
                
    print(allSubfields, set(allSubfields), subfields)
    return len(allSubfields) > len(set(allSubfields))


def repeatingPC(perims, spikes, pos, title, unitNames):
    a = np.empty(0, dtype=bool)
    for i,perim in enumerate(perims):
        a = np.hstack((a,repeatingPF(perim, spikes[i], pos)))
    with open(title, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        data = np.column_stack((unitNames, a))
        csvwriter.writerows(data)


Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")
Y, X = np.mgrid[0:450, 0:600]
coords = np.array(list(zip(X.flatten(), Y.flatten())))

rows = np.linspace(2.38, 477.62, round(480/4.72)-1) #matches the 2D histogram of spikes 
cols = np.linspace(2.37, 637.63, round(640/4.72)-1) #matches the 2D histogram of spikes
Xbin, Ybin = np.meshgrid(cols, rows)
coordsBin = np.array(list(zip(Xbin.flatten(), Ybin.flatten())))

#repeatingPC([borders1,borders2,borders4,borders5,borders6,borders7,borders8],[spikes1,spikes2,spikes4,spikes5,spikes6,spikes7,spikes8],pos,"20200802-134514",["unit1","unit2","unit4","unit5","unit6","unit7","unit8"])
