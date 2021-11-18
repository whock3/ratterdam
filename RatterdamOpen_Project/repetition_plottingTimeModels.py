# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:38:21 2021

@author: whockei1

Load glm data from R relating to time encoding
and plot. 

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ratterdam_Defaults as Def


databasepath = 'E:\\Ratterdam\\temp\TimeCD_Model\\'
fname = "alleyVisit_timeData.csv"

df = pd.read_csv(databasepath+fname)


#%% plot
fig, ax = plt.subplots()

ax.bar([0,0.5],
       [df[(df.sigs==True)&(df.repOrNot==True)].shape[0]/df[df.repOrNot==True].shape[0],
       df[(df.sigs==True)&(df.repOrNot==False)].shape[0]/df[df.repOrNot==False].shape[0]
       ],
       width=0.3
       )
       