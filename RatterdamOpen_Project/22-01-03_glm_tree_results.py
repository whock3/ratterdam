# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 21:40:36 2022

@author: whockei1
"""

from matplotlib import pyplot as plt
import numpy as np 
from matplotlib.ticker import MaxNLocator


model_labels = ["CD",
                "PD",
                "ND", 
                "CD versus CD+PD",
                "CD versus CD+ND",
                "PD versus PD+CD",
                "PD versus PD+ND",
                "ND versus ND+CD",
                "ND versus ND+PD",
                "CD+PD versus CD+PD+ND",
                "CD+ND versus CD+ND+PD",
                "PD+CD versus PD+CD+ND",
                "PD+ND versus PD+ND+CD",
                "ND+CD versus ND+CD+PD",
                "ND+PD versus ND+PD+CD"          
                ]


# model results are average number of field chunks which have a signifiant LRT
# of the indicated comparison. Not mutually exclusive - fields can have multiple
# significant test results.

repeating = [0.256,
             0.200,
             0.233,
             0.108,
             0.067,
             0.167,
             0.150,
             0.058,
             0.141,
             0.050,
             0.108,
             0.050,
             0.033,
             0.108,
             0.033
            ]

#
singlefield = [0.222,
               0.111,
               0.185,
               0.037,
               0.111,
               0.222,
               0.222,
               0.111,
               0.037,
               0.074,
               0.037,
               0.074,
               0.074,
               0.037,
               0.074
                ]


x = np.arange(len(model_labels))
width = 0.35

fig, ax = plt.subplots()
rect1 = ax.bar(x-width/2,repeating,width,label='Repeating Field Pieces',color='firebrick')
rect2 = ax.bar(x+width/2,singlefield,width,label='Single-fields',color='dodgerblue')
ax.set_ylabel("Proportion of Field Pieces With \n Significant LRT Comparison between Models",fontsize=30)
ax.hlines(0.05,-0.5,14,linestyle='--',color='k',linewidth=4)
ax.set_title("GLM LRT Analysis 22-01-04",fontsize=30)
ax.tick_params(axis='y', which='major', labelsize=18)

ax.set_xticks(range(len(model_labels)))
ax.set_xticklabels(model_labels,rotation=45,fontsize=18,ha='right')

ax.legend(fontsize=20)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.subplots_adjust(bottom=0.3)
plt.show()