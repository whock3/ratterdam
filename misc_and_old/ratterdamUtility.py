
###########
# Imports #
###########

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys
sys.path.insert(0, 'C:\\Users\\whock\\Google Drive\\KnierimLab\\Ratterdam\\Code')
sys.path.insert(0, 'C:\\Users\\whock\\Google Drive\\Python_Code\\KLab\\mts_analysis')
import utility_fx as util

alleyBounds = {0:[[395, 495],[352, 390]],
               1:[[395, 495],[218, 255]],
               2:[[495, 535],[255, 352]],
               3:[[360, 395],[255, 352]],
               4:[[260, 360],[352, 390]],
               5:[[216, 260],[255, 352]],
               6:[[120, 216],[352, 390]],
               7:[[87, 120],[255, 352]],
               8:[[120, 216],[218, 255]],
               9:[[87, 120],[115, 218]],
               10:[[120, 216],[75, 115]],
               11:[[216, 260],[115, 218]],
               12:[[260, 360],[218, 255]],
               13:[[260, 360],[75, 115]],
               14:[[360, 395],[115, 218]],
               15:[[395, 495],[75, 115]],
               16:[[495, 535],[115, 218]]}



