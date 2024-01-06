# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 14:02:39 2021

@author: whockei1

Default parameters for repetition manuscript paper
These mostly relate to plotting params, like tick sizes and so forth
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

ylabelsize = 24
xlabelsize = 24
spine_width = 1

ticksize =  16 # this is label not tick itself

titlesize = 30

legend_size = 44
legend_marker_size = 30
legend_frame_width = 2

scatter_size = 250

# set tick width
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.sans-serif'] = 'Calibri'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.labelsize'] = ticksize # in case there are figures that dont call rcparams themselves
mpl.rcParams['ytick.labelsize'] = ticksize 