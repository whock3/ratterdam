# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 19:31:38 2020

@author: Ruo-Yah Lai
"""
from matplotlib import path
import numpy as np
from RateMap import makeRM
from collections import namedtuple
from newAlleyBounds import R781, R808, R859, alleyInterType
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


def PolyArea(x,y):
    """
    Found at https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    """
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))


def repeatingPF(perims, spikes, pos, rat):
    """
    perims: subfields' perimeters
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    """
    subfields = [[] for _ in range(len(perims))]
    overlaps = [[] for _ in range(len(perims))]
    allSubfields = np.empty(0)
    threshold = 0.6
    for i,perim in enumerate(perims):
        fieldSize = PolyArea(perim[:,0], perim[:,1])
        for j in range(17):
            alley = rat.alleyInterBounds[str(j)]
            aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
            alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                   [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
            alleySize = (aRect.xmax-aRect.xmin) * (aRect.ymax-aRect.ymin)
            overlap1 = overlap(perim, alleyPerim)
            if overlap1 > alleySize*threshold or overlap1 > fieldSize*threshold:
                subfields[i].append(alleyInterType[str(j)][1])
                overlaps[i].append(j)
        for j in ascii_uppercase[:12]:
            intersection = rat.alleyInterBounds[j]
            iRect = Rectangle(intersection[0][0], intersection[1][0], 
                              intersection[0][1], intersection[1][1])
            iPerim = np.array([[iRect.xmin, iRect.ymin], [iRect.xmax, iRect.ymin],
                               [iRect.xmax, iRect.ymax], [iRect.xmin, iRect.ymax]])
            Size = (iRect.xmax-iRect.xmin) * (iRect.ymax-iRect.ymin)
            overlap1 = overlap(perim, iPerim)
            if overlap1 > Size*threshold or overlap1 > fieldSize*threshold:
                subfields[i].append(alleyInterType[j][0])
                overlaps[i].append(j)
        
        print(subfields[i])
        if len(set(subfields[i])) > 1: #if the subfield is in multiple location types
            masses = np.empty(0) #mass in each location
            for k, j in enumerate(overlaps[i]):
                alley = rat.alleyInterBounds[str(j)]
                aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
                alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                       [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
                masses = np.hstack((masses, mass(perim, alleyPerim, spikes, pos)))
            horiz = np.sum(masses[np.where(np.array(subfields[i]) == "horizontal")])
            vert = np.sum(masses[np.where(np.array(subfields[i]) == "vertical")])
            intxn = np.sum(masses[np.where(np.array(subfields[i]) == "intersection")])
            print(horiz,vert,intxn)
        
            subfields[i] = []
            if horiz > 0.2:
                subfields[i].append("horizontal")
            if vert > 0.2:
                subfields[i].append("vertical")
            if intxn > 0.2:
                subfields[i].append("intersection")
        
        for j in subfields[i]:
            allSubfields = np.hstack((allSubfields, j))
    
    repeat = False
    PFtype = ""
    a = 0
    if len(np.where(allSubfields == "horizontal")[0]) > 1:
        repeat = True
        PFtype = PFtype + "horizontal"
        a += 1
    if len(np.where(allSubfields == "vertical")[0]) > 1:
        repeat = True
        PFtype = PFtype + "vertical"
        a += 1
    if len(np.where(allSubfields == "intersection")[0]) > 1:
        repeat = True
        PFtype = PFtype + "intersection"
        a += 1
    
    repeatType = ""
    if a > 1: #more than 1 type of location is repeated
        repeatType = "multiple"
        multLoc = [i for i in subfields if len(i)>1]
        for i,subfield1 in enumerate(multLoc):
            for j,subfield2 in enumerate(multLoc):
                if i == j:
                    continue
                if len(set(subfield1 + subfield2)) < len(subfield1 + subfield2)-1:
                    repeatType = "complex"
                    break
            else:
                continue
            break
    
    print(allSubfields, subfields)
    return repeat, PFtype, repeatType, overlaps


#units from UnitClass2
def repeatingPC(units, pos, title, unitNames, rat="", day="", tetrode=""):
    repeats = np.empty(0, dtype=bool)
    PFtypes = np.empty(0, dtype=str)
    repeatTypes = np.empty(0, dtype=str)
    for unit in units:
        repeat, PFtype, repeatType, _ = repeatingPF(unit.perimeters, unit.spikes, pos)
        repeats = np.hstack((repeats,repeat))
        PFtypes = np.hstack((PFtypes, PFtype))
        repeatTypes = np.hstack((repeatTypes, repeatType))
    rows = len(unitNames)
    with open(title+'.csv', "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        data = np.column_stack((np.full(rows,rat), np.full(rows,day), np.full(rows,tetrode), unitNames, repeats, PFtypes, repeatTypes))
        csvwriter.writerow(["Rat number", "Day", "Tetrode number", "Unit", "Repetition?", "Repeated location", "Repeat type"])
        csvwriter.writerows(data)


def fieldLocation(unit, title):
    """
    Prints the alleys/intersections that each field occupies into a csv
    """
    _,_,_,overlaps = repeatingPF(unit.perimeters,unit.spikes,unit.position)
    with open(title, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Subfield", "Alleys/intersections"])
        for i,a in enumerate(overlaps):
            a = [str(x) for x in a]
            csvwriter.writerow([i,a])


Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")
Y, X = np.mgrid[0:450, 0:600]
coords = np.array(list(zip(X.flatten(), Y.flatten())))

rows = np.linspace(2.38, 477.62, round(480/4.72)-1) #matches the 2D histogram of spikes 
cols = np.linspace(2.37, 637.63, round(640/4.72)-1) #matches the 2D histogram of spikes
Xbin, Ybin = np.meshgrid(cols, rows)
coordsBin = np.array(list(zip(Xbin.flatten(), Ybin.flatten())))

#repeatingPC([unit1,unit2,unit4,unit5,unit6,unit7,unit8],pos,genTimestamp(),["unit1","unit2","unit4","unit5","unit6","unit7","unit8"], "R859", "D3", "T6")
