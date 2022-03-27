# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:41:20 2022

@author: whockei1


Quantifying time dynamics within CA1 neurons recorded on city-block maze (Ratterdam)

"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
import json, pickle
import ratterdam_RepetitionCoreFx as RepCore
import repetition_DfSamplingFiltering as DfFilt 
import repetition_manuscript_defaults as MDef 


df = pd.read_csv("E:\\Ratterdam\\R_data_repetition\\2022-03-23_AlleySuperpopDirVisitFiltered.csv")
filtDf = DfFilt.filterAlleyDatasets(df, passThresh=1)


#%% Panel A - example time signals
# taken from raw data, no analysis to do here

#%% Panel B - GLM results from R. repetition_timeModels.R

glmData = pd.read_csv("E:\\Ratterdam\\repetition_manuscript\\Figure6\\2022-03-23_timeGLMResults.csv")

fig, ax = plt.subplots()

ax.plot(glmData.rmse_base[glmData.sigs==0],glmData.rmse_alt[glmData.sigs==0],
        marker='o',
        markersize=10,
        linestyle='',
        color='black',
        label='Not Temporally Dynamic')
ax.plot(glmData.rmse_base[glmData.sigs==1],glmData.rmse_alt[glmData.sigs==1],
        marker='o',
        markersize=10,
        linestyle='',
        color='red',
        label='Temporally Dynamic')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(MDef.spine_width)
ax.spines['bottom'].set_linewidth(MDef.spine_width)
ax.set_ylabel("Base Model \n+ Time RMSE",fontsize=MDef.ylabelsize)
ax.set_xlabel("Base Model RMSE",fontsize=MDef.xlabelsize)
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', which='major', labelsize=MDef.ticksize)
lgnd = plt.legend(prop={'size':MDef.legend_size})
for lhand in lgnd.legendHandles:
    lhand._legmarker.set_markersize(MDef.legend_marker_size)
lgnd.get_frame().set_linewidth(MDef.legend_frame_width)



#%% Decode time using k-NN using raw time, i.e. time in us from start, not IFD vectors

from scipy.interpolate import PchipInterpolator as pchip
from bisect import bisect_left
from sklearn.neighbors import KNeighborsClassifier as KNN 
from sklearn.model_selection import train_test_split 
#%% 
#unit = RepCore.loadRepeatingUnit('R886','D1','TT5\\cl-maze1.3')
units, turns = RepCore.loadRecordingSessionData('R765','DFD4')
#%% 

# Parameters
shuffle = False
total_points = 500 
padThresh = 1/total_points # what we want is the time gap between start/stop of
                            # field and beginning/end session that is so small
                            # we cant even add pad with 1 point. 
nEpochs = 3 # number of epochs to divide the session into, these are the labels we are trying to decode
wnSize = 5*60*1e6 # time in minutes of window size 
wnStep = 2*60*1e6

nTechReps = 25
nShuffles = 1

shuffle_meanperfs = []


for s in range(nShuffles):
    pchipfields = []
    reps = []


    print(s)
    # Gather all fields, interpolate within using pchip, and pad to ends of session w closest value
    fieldArray = [u.fields for u in units.values()]
    fieldArray = [item for sublist in fieldArray for item in sublist]
    
    if shuffle:
        fieldArray = [np.column_stack((field[:,0], np.random.permutation(field[:,1]))) for field in fieldArray]
    
    fmax = int(np.ceil(max([max(field[:,0]) for field in fieldArray])))
    fmin = int(np.ceil(min([min(field[:,0]) for field in fieldArray])))
    total_length = (fmax-fmin)/1e6
                                                          
    session_time_vector = np.linspace(fmin,fmax,total_points)
    
    for field in fieldArray:
        
        beginningGap = ((field[0,0]-fmin)/1e6)/total_length
        endGap = ((fmax-field[-1,0])/1e6)/total_length
        
        if beginningGap > padThresh:
            npadbeginning = int(np.ceil(beginningGap*total_points))
            padbeginning = [field[0,1]]*npadbeginning
        else:
            npadbeginning = 0
            padbeginning = []
            
        if endGap > padThresh:
            npadend = int(np.ceil(endGap*total_points))
            padend = [field[-1,1]]*npadend
        else:
            npadend = 0
            padend = []
            
            
        nptsfield = total_points - npadbeginning - npadend
        
        time = np.linspace(int(field[0,0]), int(field[-1,0]),nptsfield)
        pchipfields.append(np.concatenate((padbeginning,
                                           pchip(field[:,0], field[:,1])(time),
                                           padend
                                           ))
                           )
                           
                    
    pchipfields = np.asarray(pchipfields) # n,p array n=number fields, p=total_points
    
    # plt.figure()                     
    # for rawfield,interpfield in zip(fieldArray, pchipfields):
    #     p = plt.plot(session_time_vector, interpfield,linestyle='-',marker='o',markersize=5)
    #     plt.plot(rawfield[:,0], rawfield[:,1], linestyle='-', color=p[0].get_color())
    
    # plt.vlines(session_time_vector,0,plt.ylim()[1],color='k',alpha=0.6)
    # Decode time using k-nn approach
    epoch_intervals = np.linspace(fmin-1, fmax, nEpochs+1) # fmin-1 because if the edge is equal to the value bisect_left doesnt work as I want it to, so make the interval ever so slightly larger towards the left
    
    wins = []
    winIndices = []
    labels = []
    begin = fmin
    stop = False
    
    while not stop:
        a,b = begin, begin + wnSize
        if b < np.ceil(fmax):
            wins.append((a,b))
            begin += wnStep
        else:
            stop = True
    
    popvecs = []
    for idx,win in enumerate(wins):
        winIndices.append(idx)
        labels.append(bisect_left(epoch_intervals, int(np.mean([win[0],win[1]]))))  
        finw=pchipfields[:,np.where((session_time_vector > win[0]) & (session_time_vector <= win[1]))[0]]
        popvec = np.mean(finw,axis=1)
        popvecs.append(popvec)
        
    popvecs = np.asarray(popvecs)
    
    for n in range(nTechReps):
        Xtrain, Xtest, Ytrain, Ytest = train_test_split(popvecs,labels)
        neigh = KNN(n_neighbors=3)
        neigh.fit(Xtrain,Ytrain)
        neigh.predict(Xtest)
        reps.append(sum(neigh.predict(Xtest)==Ytest)/len(Xtest))
    
    meanperf = np.mean(reps)
    if shuffle:
        shuffle_meanperfs.append(meanperf)
    