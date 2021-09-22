# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 19:31:38 2020

@author: Ruo-Yah Lai
"""
from matplotlib import path
import numpy as np
from collections import namedtuple, Counter
from newAlleyBounds import alleyInterType
import utility_fx as util
from RateMap import makeRM
import ratterdam_RepetitionCoreFx as RepCore
import csv
from string import ascii_uppercase
import ratterdam_Defaults as Def

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
    hist = makeRM(spikes, pos, smoothing_2d_sigma=0) #different from util.makeRM
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


def repeatingPF(unit, rat):
    """
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    """
    subfields = [[] for _ in range(len(unit.perimeters))]
    overlaps = [[] for _ in range(len(unit.perimeters))]
    subfieldsAbbr = []
    threshold = Def.fieldOverlapThresh # expressed as decimal 0-1 of % region area field overlaps w 
    for i,perim in enumerate(unit.perimeters):
        fieldSize = PolyArea(perim[:,0], perim[:,1])
        for j in range(17):
            alley = rat.alleyInterBounds[str(j)]
            aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
            alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                   [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
            alleySize = (aRect.xmax-aRect.xmin) * (aRect.ymax-aRect.ymin)
            overlap1 = overlap(perim, alleyPerim)
            if overlap1 > alleySize*threshold: # or overlap1 > fieldSize*threshold: #WH commented out second conditional 7-14-21
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
            if overlap1 > Size*threshold: # or overlap1 > fieldSize*threshold: WH commented out 8-10-21. Should have comment out when line 65 (ie same conditional but for alleys) but didnt notice this
                subfields[i].append(alleyInterType[j][0])
                overlaps[i].append(j)
    
        
        #print(subfields[i])
        abbr = {"horizontal": "H", "vertical": "V", "intersection": "I"}
        if len(set(subfields[i])) == 1: #if the subfield is in one location type
            subfieldsAbbr.append(abbr[subfields[i][0]])
        elif len(set(subfields[i])) > 1: #if the subfield is in multiple location types
            masses = np.empty(0) #mass in each location
            for k, j in enumerate(overlaps[i]):
                alley = rat.alleyInterBounds[str(j)]
                aRect = Rectangle(alley[0][0], alley[1][0], alley[0][1], alley[1][1])
                alleyPerim = np.array([[aRect.xmin, aRect.ymin], [aRect.xmax, aRect.ymin],
                                       [aRect.xmax, aRect.ymax], [aRect.xmin, aRect.ymax]])
                masses = np.hstack((masses, mass(perim, alleyPerim, unit.spikes, unit.position)))
            horiz = np.sum(masses[np.where(np.array(subfields[i]) == "horizontal")])
            vert = np.sum(masses[np.where(np.array(subfields[i]) == "vertical")])
            intxn = np.sum(masses[np.where(np.array(subfields[i]) == "intersection")])
            #print(horiz,vert,intxn)
        
            #subfields[i] = []
            subfieldsAbbr.append("")
            if horiz > 0.2:
                #subfields[i].append("horizontal")
                subfieldsAbbr[-1] += "H"
            if vert > 0.2:
                #subfields[i].append("vertical")
                subfieldsAbbr[-1] += "V"
            if intxn > 0.2:
                #subfields[i].append("intersection")
                subfieldsAbbr[-1] += "I"
        
        #for j in subfields[i]:
        #    allSubfields = np.hstack((allSubfields, j))
    
    repeatTypeNumber = 0
    repeatType = ""
    repeat = False
    locCount = Counter(H=0, V=0, I=0, HV=0, HI=0, VI=0, HVI=0)
    locCount.update(subfieldsAbbr)
    for key, value in locCount.items():
        if value > 1:
            repeat = True
            repeatTypeNumber += 1
            if len(key) > 1:
                repeatType = "complex"
    if repeatTypeNumber > 1:
        repeatType += "multiple"
    
        
    #repeat = False
    #PFtype = ""
    #a = 0
    #if len(np.where(allSubfields == "horizontal")[0]) > 1:
    #    repeat = True
    #    PFtype = PFtype + "horizontal"
    #    a += 1
    #if len(np.where(allSubfields == "vertical")[0]) > 1:
    #    repeat = True
    #    PFtype = PFtype + "vertical"
    #    a += 1
    #if len(np.where(allSubfields == "intersection")[0]) > 1:
    #    repeat = True
    #    PFtype = PFtype + "intersection"
    #    a += 1
    
    #repeatType = ""
    #if a > 1: #more than 1 type of location is repeated
    #    repeatType = "multiple"
    #    multLoc = [i for i in subfields if len(i)>1]
    #    for i,subfield1 in enumerate(multLoc):
    #        for j,subfield2 in enumerate(multLoc):
    #            if i == j:
    #                continue
    #            if len(set(subfield1 + subfield2)) < len(subfield1 + subfield2)-1:
    #                repeatType = "complex"
    #                break
    #        else:
    #            continue
    #        break
    
    #print(subfieldsAbbr)
    return repeat, locCount, repeatType, overlaps


#units from UnitClass2
def repeatingPC(units, title, filepath, unitNames, rat):
    """
    units: list of units
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    
    Prints whether there is repetition, location of repetition, and type of
    repetition of each unit into a csv
    """
    with open(filepath+title+'.csv', "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        #data = np.column_stack((unitNames, repeats, subfieldsAbbrs, repeatTypes))
        csvwriter.writerow(["Unit", "Repetition?", "Repeated location", "Repeat type"])
        for i,unit in enumerate(units):
            repeat, locCount, repeatType, _ = repeatingPF(unit, rat)
            PFtype = ""
            for loc in locCount:
                if locCount[loc] > 1:
                    PFtype = PFtype + loc + ", "
            if len(PFtype) > 0:
                PFtype = PFtype[:-2]
            csvwriter.writerow([unitNames[i], repeat, PFtype, repeatType])
            

def repeatingPCAllLoc(units, title, filepath, unitNames, rat):
    """
    units: list of units
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    
    Prints whether there is repetition, lists all locations and location combinations
    , and type of repetition of each unit into a csv
    """
    with open(filepath+title+'.csv', "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        #data = np.column_stack((unitNames, repeats, subfieldsAbbrs, repeatTypes))
        csvwriter.writerow(["Unit", "Repetition?", "H","V","I","HV","HI","VI","HVI", "Repeat type"])
        for i,unit in enumerate(units):
            repeat, c, repeatType, _ = repeatingPF(unit, rat)
            csvwriter.writerow([unitNames[i], repeat, c["H"], c["V"], c["I"], c["HV"], c["HI"], c["VI"], c["HVI"], repeatType])


def fieldLocation(unit, title, rat, filepath):
    """
    Prints the alleys/intersections that each field occupies into a csv
    """
    _,_,_,overlaps = repeatingPF(unit, rat)
    with open(filepath+title, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Subfield", "Alleys/intersections"])
        for i,a in enumerate(overlaps):
            a = [str(x) for x in a]
            csvwriter.writerow([i,a])


def bulkFieldLocation(file, title, timestamp, rat, df, df2, df3):
    """
    file: file name of the file with repetition tabulations
    title: part of the title of the generated files (include rat and day)
    rat: named tuple R781, R808, or R859 from newAlleyBounds
    df: path to data
    df2: path to repetition tabulations
    df3: where to save the files
    """
    cells = util.readRepeatingCells(file, df2)
    for cell in cells:
        cellTitle = cell.replace("\\cl-maze", " ")
        unit = RepCore.loadRepeatingUnit(df, cell)
        fieldLocation(unit, f"{timestamp} - {title} {cellTitle} locations", rat, df3)


Rectangle = namedtuple("Rectangle", "xmin ymin xmax ymax")
Y, X = np.mgrid[0:450, 0:600]
coords = np.array(list(zip(X.flatten(), Y.flatten())))

rows = np.linspace(2.38, 477.62, round(480/4.72)-1) #matches the 2D histogram of spikes 
cols = np.linspace(2.37, 637.63, round(640/4.72)-1) #matches the 2D histogram of spikes
Xbin, Ybin = np.meshgrid(cols, rows)
coordsBin = np.array(list(zip(Xbin.flatten(), Ybin.flatten())))

#repeatingPC([unit1,unit2,unit4,unit5,unit6,unit7,unit8],pos,genTimestamp(),["unit1","unit2","unit4","unit5","unit6","unit7","unit8"], "R859", "D3", "T6")
