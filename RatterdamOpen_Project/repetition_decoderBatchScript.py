# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 14:13:51 2021

@author: whockei1

Wrapper script to run classifiers (RF and naive) on data 
"""

import numpy as np, pandas as pd, copy
from sklearn.ensemble import RandomForestClassifier
import ratterdam_Defaults as Def
import ratterdam_RepetitionCoreFx as RepCore
import utility_fx as util
import json, os 
from matplotlib import pyplot as plt 
import repetition_decoderFunctions as Dfx

rats = ['R765','R781','R781','R808','R808','R859','R859','R886','R886']
days = ['RFD5','D3','D4','D6','D7','D1','D2','D1','D2']

verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]

population = Dfx.loadNeuralData(rats, days)
results = {}
for target in ['CurrDir']:
    print(target)
    for rat, day in zip(rats, days):
        print(rat, day)
        
        if rat not in results.keys():
            results[rat] = {}
        if day not in results[rat].keys():
            results[rat][day] = {}
        if target not in results[rat][day].keys():
            results[rat][day][target] = {}

        turns, refturns = population[rat][day]['turns'], population[rat][day]['refturns']
        labels = Dfx.createLabelArrayAlley(turns, refturns, target)
        Xrep = Dfx.createNeuralResponseArray(population[rat][day]['units'], labels, group='Repeating')
        Xnonrep = Dfx.createNeuralResponseArray(population[rat][day]['units'], labels, group='Non-Repeating')
        
        #some days may have no repeating neurons by chance. check nonrep size just to be safe though not necessary 
        if Xrep.shape[1] > 0 and Xnonrep.shape[1] > 0:
            results[rat][day][target]['naivePerf'] = {}
            results[rat][day][target]['naivePerf']['verticals'] = Dfx.naiveDirectionClassifier(turns,verticals)
            results[rat][day][target]['naivePerf']['horizontals'] = Dfx.naiveDirectionClassifier(turns,horizontals)
            
            for repType, _X in zip(['Repeating', 'Non-Repeating'],[Xrep,Xnonrep]):
                if repType not in results[rat][day][target].keys():
                    results[rat][day][target][repType] = {}
                for orienType, orienList in zip(['verticals', 'horizontals'],[verticals, horizontals]):
                    oobs = Dfx.runRandomForest(_X[labels['CurrAlley'].isin(orienList)],
                                               labels[labels['CurrAlley'].isin(orienList)]['Label']                                                
                                                )
                    results[rat][day][target][repType][orienType] = oobs