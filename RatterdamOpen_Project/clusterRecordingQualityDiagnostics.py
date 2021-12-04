# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 17:57:13 2021

@author: whockei1
"""
import numpy as np 
import utility_fx as util 
import ratterdam_RepetitionCoreFx as RepCore 
import ratterdam_Defaults as Def 
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages


numChannels = 4

parmsType = np.dtype([
        ('ts', '<u8'),
        ('data', '<f4', (numChannels*4+6))
    ])


nttType = np.dtype([
    ('ts', '<u8'),
    ('spikeEntryNumber', '<u4'),
    ('cellNumber', '<u4'),
    ('parameters', '<u4', (8)),
    ('samples', '<i2', (32,4))
])

def load_NTT(filename, count=-1, headerRaw=False):
    """
    Load Neuralynx NTT file
    Parameters
    ----------
    filename: str
        The full path of the file
    count: integer
        Number of instances to load
    Returns
    -------
    ntt : Numpy array
        Numpy array using the nvtType dtype
    """
    if headerRaw:
        nttFile = open(filename, 'rb')
        header = nttFile.read(16*1024)
        nttFile.close()
    else:
        nttFile = open(filename, 'r', encoding='mbcs')
        header = nttFile.read(16*1024)
        nttFile.close()
    nttFile = open(filename, 'rb')
    nttFile.seek(16*1024) #skip header
    return np.fromfile(nttFile, dtype=nttType, count=count), header

#%% Generate pdf for each recording day. Looking at events from a cluster
# plotting all four wires' peaks over time on one graph per cluster

rat, day = 'R886', 'D2'
datapath = f"E:\\Ratterdam\{rat}\\{rat}_RatterdamOpen_{day}\\"
#savepath = f"E:\Ratterdam\{rat}\clusterDriftDiagnostics\\{day}\\"
savepath = 'E:\\Ratterdam\\rawClustersOverTime\\'

if not os.path.exists(savepath):
    os.makedirs(savepath)
    
with open(datapath+"sessionEpochInfo.txt","r") as f:
    lines = f.readlines()
start, end = float(lines[0].split(',')[0]), float(lines[0].split(',')[1])

with PdfPages(savepath+f"{rat}{day}_clustParmsOverTime.pdf") as pdf:
    for subdir, dirs, fs in os.walk(datapath):
        clusters = {}
        tt = subdir.split('\\')[-1]

        for f in fs:
            if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f and "manuallyRedrawnFields" not in f:
                clusters[f"{tt}{f}"] = {}
                with open(subdir+"\\"+f) as clust:
                    raw_clust_data = clust.read().splitlines()
                ts = [float(i.split(',')[-1]) for i in raw_clust_data[13:]] 
                x = [int(i.split(',')[1]) for i in raw_clust_data[13:]] # data begins line 13. 
                y = [int(i.split(',')[2]) for i in raw_clust_data[13:]] # ts is 0 element, then x,y,a,b peaks
                a = [int(i.split(',')[3]) for i in raw_clust_data[13:]]
                b = [int(i.split(',')[4]) for i in raw_clust_data[13:]]
                
                clusters[f"{tt}{f}"]["ts"] = ts
                clusters[f"{tt}{f}"]["x"] = x
                clusters[f"{tt}{f}"]["y"] = y
                clusters[f"{tt}{f}"]["a"] = a
                clusters[f"{tt}{f}"]["b"] = b
                

                
        nclusts = len(clusters)
        if nclusts >= 1:
            fig, ax = plt.subplots(nclusts,1,figsize=(12,10))
            for i,clust in enumerate(clusters.keys()):
                ts = clusters[clust]['ts']
                for wire,color in zip(['x','y','a','b'],['b','g','r','k']):
                    fig.axes[i].scatter(ts,clusters[clust][wire],c=color,rasterized=True)
                fig.axes[i].set_xlim([start, end])
                fig.axes[i].set_title(clust)
            plt.suptitle(tt)
            
            pdf.savefig()
            plt.close()
                    
                                                                
                


#%% Generate pdf for each recording day, looking at all triggered events. 
# Each page is a tetrode. Each plot is a wire (X,Y,A,B) vs time

Def.velocity_filter_thresh = 0
rat, day = 'R859', 'D2'
datapath = f"E:\\Ratterdam\{rat}\\{rat}_RatterdamOpen_{day}\\"
savepath = f"E:\Ratterdam\{rat}\clusterDriftDiagnostics\\{day}\\"
if not os.path.exists(savepath):
    os.makedirs(savepath)
    
with open(datapath+"sessionEpochInfo.txt","r") as f:
    lines = f.readlines()
start, end = float(lines[0].split(',')[0]), float(lines[0].split(',')[1])

with PdfPages(savepath+f"{rat}{day}_nttParmsOverTime.pdf") as pdf:
    for subdir, dirs, fs in os.walk(datapath):
        cutCells = False
        for f in fs:
            if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f and "manuallyRedrawnFields" not in f:
                cutCells = True
            if 'ntt' in f and 'parms' not in f:
                ntt,header = load_NTT(subdir+'\\'+f)
                print('loading ntt')
                
                fig, ax = plt.subplots(4,1,figsize=(12,10))
                ts,maxx,maxy,maxa,maxb = [],[],[],[],[]
                for sample in ntt:
                    ts.append(sample[0])
                    x,y,a,b = sample[4].max(axis=0) # (32,4) array is the spike
                    maxx.append(x)
                    maxy.append(y)
                    maxa.append(a)
                    maxb.append(b)
                for i,wire in enumerate([maxx,maxy,maxa,maxb]):
                    fig.axes[i].scatter(ts,wire,c='k',s=3,rasterized=True)
                    fig.axes[i].set_xlim([start,end])
                    fig.axes[i].set_title(f"Wire {i}")
                tt = subdir.split('\\')[-1]
                plt.suptitle(f"Tetrode {tt}, cut cells = {cutCells}")
                
                pdf.savefig()
                plt.close()
                
#%% Generate plots of x,y position (separate subplots) over time for each neuron
# one tt per page, one neuron per row
Def.velocity_filter_thresh = 0
rat, day = 'R886', 'D1'
datapath = f"E:\\Ratterdam\{rat}\\{rat}_RatterdamOpen_{day}\\"
savepath = f"E:\Ratterdam\{rat}\clusterDriftDiagnostics\\{day}\\"
if not os.path.exists(savepath):
    os.makedirs(savepath)
    
with open(datapath+"sessionEpochInfo.txt","r") as f:
    lines = f.readlines()
start, end = float(lines[0].split(',')[0]), float(lines[0].split(',')[1])

clustList, clustQuals = util.getClustList(datapath, quals=False)

tts = np.unique([i.split("\\")[0] for i in clustList])

recordingData = {i:[] for i in tts}

for clust in clustList:
    print(clust)
    tt = clust.split("\\")[0]
    unit = RepCore.loadRepeatingUnit(rat, day, clust, smoothing=1)                                   
    recordingData[tt].append(unit)
    

with PdfPages(savepath+f"{rat}{day}_cutClustersDrift.pdf") as pdf:

    for tt in tts:
        nclusts = len(recordingData[tt])
        fig, ax  = plt.subplots(nclusts,2,figsize=(15,12))
        for i,unit in enumerate(recordingData[tt]):
            
            ax[i,0].plot(unit.spikes[:,0], unit.spikes[:,1],color='k',marker='.',linestyle='')
            ax[i,0].set_xlim([start, end])
            
            ax[i,1].plot(unit.spikes[:,0], unit.spikes[:,2],color='k',marker='.',linestyle='')
            ax[i,1].set_xlim([start, end])
            
            ax[i,0].set_title(f"{unit.name} x position")
            ax[i,1].set_title(f"{unit.name} y position")
            
        plt.suptitle(f"Tetrode {tt}",fontsize=20)
        
        pdf.savefig()
        plt.close()

    
