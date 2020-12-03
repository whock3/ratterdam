# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 20:38:15 2020

@author: whockei1
"""

import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
from scipy.stats import sem
import utility_fx as util
import os
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
import ratterdam_RepetitionCoreFx as RepCore
import RateMapClass_William_20190308 as RateMapClass
import williamDefaults as wmDef
from matplotlib import path
from matplotlib.backends.backend_pdf import PdfPages
import more_itertools, itertools
from sklearn.metrics import auc
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from scipy.interpolate import splrep, splev
from scipy.spatial import ConvexHull
import scipy


rat = "R859"
day = "D2"
savepath = f'E:\\Ratterdam\\{rat}\\ratterdam_plots\\{day}\\decoding\\'
df = f'E:\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
clustList = util.getClustList(df)
population = {}
for clust in clustList:
    print(clust)
    unit = RepCore.loadRepeatingUnit(df, clust, smoothing=1)
    rm = util.makeRM(unit.spikes, unit.position)
    if np.nanpercentile(rm, 95) > 1.:
        population[clust] = unit
        print(f"{clust} included")
    else:
        print(f"{clust} is not included")