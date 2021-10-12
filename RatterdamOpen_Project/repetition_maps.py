# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:41:42 2021

@author: whockei1

Map switching across the population 
"""

import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl
import pickle

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)
    