# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 11:19:19 2021

@author: whockei1

Script to find repeated behavioral patterns in data, i.e. stereotyped behaviors
This is a qualitative analysis meant to uncover repeated patterns in assumption-less way

General approach is to take a sliding window and define behavioral snippet in
it as the template Slide across the session and compute similarity of new snippet
in each new window to original template. Similarity fx will be a string matching
like Hamming distance or a generalization like Damerau–Levenshtein distance

Each sliding window will serve as a template. So the 1st window is the first template,
slide it and compute similarity at each step. Then second window is template and repeat,
but starting at second window position. Add each as a new row in what will be a jagged,
increasingly short row-length array

Repeat with increasing window sizes, each window size will genrate a separate plot

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json, pickle
import ratterdam_Defaults as Def 
import ratterdam_RepetitionCoreFx as RepCore
from scipy.spatial.distance import hamming



rat, day = 'R859', 'D2'
turns, unit  = RepCore.loadTurns(rat, day)
alleys = list(turns["Alley+"])
headings = list(turns['Allo+'])
session_trajectory = pd.Series([f"{a}-{b}" for a,b in zip(alleys,headings)])



#%% define functions

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

#%% forward looking only 

window_length = 5

distances = []  

for i, template in enumerate(list(session_trajectory.rolling(window_length))[window_length-1:]):
    single_row = np.zeros((1,i))    
    for j, snippet in enumerate(list(session_trajectory.rolling(window_length))[window_length-1+i:]):
        dist = abs(window_length-levenshteinDistance(template, snippet))
        single_row = np.append(single_row, dist)
    distances.append(single_row)
    
distances = np.asarray(distances)


#%% forward and reverse looking  

window_length = 10

distances = []  

for i, template in enumerate(list(session_trajectory.rolling(window_length))[window_length-1:]):
    single_row = []   
    for j, snippet in enumerate(list(session_trajectory.rolling(window_length))[window_length-1:]):
        dist = abs(window_length-levenshteinDistance(template, snippet))
        single_row.append(dist)
    distances.append(single_row)
    
distances = np.asarray(distances)