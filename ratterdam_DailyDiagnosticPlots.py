# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:39:06 2019

@author: whockei1
"""
import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
import utility_fx as util
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis


datafile = "E:\\Ratterdam\\R808\\R808_Beltway_D6\\2019-08-20_18-54-16\\"
expCode = "BRD6"

figpath = f"E:\\Ratterdam\\R808\\beltway_plots\\{expCode}\\"
beltwayAlleys = [16,17,3,1,5,7,8,10,11]

ratplot = Vis.BasicRateMaps()

alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)

stimfiles = Parse.getStimFileName(datafile)
stimData = Parse.loadBeltwayData(datafile,stimfiles,expCode)

cmap = util.makeCustomColormap()


alley_name_conversion = {16:1, 17:2, 3:3, 1:4, 5:5, 7:6, 8:7, 10:8, 11:9}



with PdfPages(figpath+f"{expCode}_Diagnostics.pdf") as pdf:
    
    
    fig, ax  =plt.subplots(10,6,figsize=(10,10))
    for lap in range(60):
        for i,a in enumerate(beltwayAlleys):
            a = a-1
            begin, end = alleyVisits[a][lap][0], alleyVisits[a][lap][1]
            p = p_sess[(p_sess[:,0]>begin)&(p_sess[:,0]<=end)]
            fig.axes[lap].scatter(p[:,1],p[:,2],c=np.random.rand(3,1))
            fig.axes[lap].set_yticklabels([])
            fig.axes[lap].set_xticklabels([])
    for l in range(60):
        for i in range(1,18):
            ratplot.drawAlleyBounds(fig.axes[l],i)
            
    pdf.savefig()
    plt.close()
    
    
