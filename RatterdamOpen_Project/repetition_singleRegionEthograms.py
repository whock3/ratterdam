# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 14:15:26 2022

@author: whockei1

Create ethograms for individual regions: alleys and intersections

Behaviors are defined as three consecutive headings: previous,current,next dirs
Idea is to then correlate or regress this with a measure of neural variability
(specifically overdispersion).

"""

#%% Imports and Load Data

import numpy as np, matplotlib.pyplot as plt, pickle, os, pandas as pd

# Load csvs with data that is less filtered than what has been used in R GLMs
# and other recent analyses. Basically this csv filters individual passes through
# field to make sure the visit "goes in enough" to the field. The second stage,
# filters data at each alley to make sure there is enough sampling in each direction
# Below is just the first stage:
    
alley_data_path = 'E:\\Ratterdam\\R_data_repetition\\20220120-135920_superPopAlleyBehaviorResponse_1.5vfilt.csv'
inter_data_path = 'E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv'

with open(alley_data_path, "r") as afile:
    alleydf = pd.read_csv(afile)

with open(inter_data_path, "r") as ifile:
    interdf = pd.read_csv(ifile)


#%% 

rat, day = 'R781', 'D3'
rd_df = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]
