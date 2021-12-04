# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 20:20:19 2021

@author: whockei1
"""
import numpy as np, matplotlib.pyplot as plt

x1 = np.random.normal(51,3,100)
x2 = np.random.normal(55,3,100)
x3 = np.random.normal(60,3,100)
x = np.concatenate((x1,x2,x3))
fig, ax = plt.subplots()
ax.plot(x,color='k')
ax.set_xlabel("Time (s)", fontsize=45)
ax.set_ylabel("Pip Frequency (Hz)", fontsize=45)
ax.hlines(50,0,300,color='r',alpha=0.6,linestyle='--')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=30)
