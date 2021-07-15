# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 23:33:13 2021

@author: whockei1

Scrap file
"""
#%% Functions to load
def checkCheating(unit, turn, turn_np1):
    """
    deprecated. solving problem a different way
    
    Check to see if any of the unit's fields overlap with either of the non-local
    regions involved in this turn (i.e. the pre and post alleys were decoding dir
     from) and only if not do we let this unit contribute data to this turn
    """
    threshold = 0.10
    anyExcessiveOverlap = False
    for perim in unit.perimeters:
        for a in [turn['Alley-'], turn_np1['Alley+']]:
            alley = ratborders.alleyInterBounds[str(a)]
            xmin,xmax = alley[0][0],alley[0][1]
            ymin,ymax = alley[1][0],alley[1][1]
            alleyPerim = np.array([[xmin, ymin], [xmax, ymin],
                                           [xmax, ymax], [xmin, ymax]])
            alleySize = (xmax-xmin) * (ymax-ymin)
            overlap1 = repPC.overlap(perim, alleyPerim)
            if overlap1 > alleySize*threshold:
                anyExcessiveOverlap = True
            
    return anyExcessiveOverlap

