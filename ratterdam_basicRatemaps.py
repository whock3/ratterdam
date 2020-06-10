# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:32:48 2018

@author: whockei1

Script to create multipage PDF for a unit with the following plots:
    Pg 1 - overall ratemap
    (Possible pg 2) - overall rate maps stiched together from alleys
    Pg 2 - Track ratemap of alleys when txt A
    pg 3 - Same, B
    pg 4 - Same C
    pg 5 - Lin Ratemaps Common Scale Y-Lim
    pg 6 - Lin Ratemaps Alleys Scaled Y-lim independently
    pg 7 - permutation test results
"""

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket, os
from scipy.stats import sem
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as GS
from matplotlib.backends.backend_pdf import PdfPages
from importlib import reload

if socket.gethostname() == 'Tolman':
    codeDirBase = 'C:\\Users\\whockei1\\Google Drive'
elif socket.gethostname() == 'DESKTOP-BECTOJ9':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'
    
sys.path.insert(0, codeDirBase + '\\KnierimLab\\Ratterdam\\Code')
sys.path.insert(0, codeDirBase + '\\Python_Code\\KLab\\mts_analysis')

import utility_fx as util
import ratterdam_ParseBehavior as pBehav
import ratterdam_CoreDataStructures as core
import ratterdam_visBasic as Vis
import ratterdam_PermutationTests as Perm
from ratterdam_Defaults import *


exp = "RFD7"
dayCode = f"R765{exp}\\"
figpath = f"C:\\Users\\whockei1\\Google Drive\\KnierimLab\\Ratterdam\\Figures\\R765{exp}\\velocityFiltered_{velocity_filter_thresh}\\"
if not os.path.isdir(figpath):
    os.mkdir(figpath)
datafile = f'E:\\Ratterdam\\R765\\R765{exp}\\'
behav = core.BehavioralData(datafile, exp, velocity_filter_thresh)
ts, position, alleyTracking, alleyVisits,  txtVisits = behav.loadData()
plot = Vis.BasicRateMaps()
gs = GS.GridSpec(5,7)
gs.update(wspace=0.025, hspace=0.05)

print(dayCode)

#################################
# Plotting Defaults and Lookups #
#################################

### Begin Contiguous Alley plotting defaults and fx
"""
I drew ratterdam on paper s.t. the alleys were 2 units short axis and
4 units long axis (the exact num is somewhat arb but ratio is important)
Then you scale these by the total length along each dim so they equal 1
(B/c were plotting subplots in data frame coords [0,1] each dim)
Then you have some lookup dicts and fx so you can go from alley ID (e.g. alley 1)
to what row, col entry it will be when plotting (here 2, 0). NB the vert and
horiz are separate because it was easier to deal with their perpendicularity
in frame coords by considering them separately (i.e. which axis is 2 or 4 is flipped)
Code in this script will loop over the alleys and plot the horiz and vert
alleys to scale using this info.

"""
hoffset = 2/22 # real sum is 20 units wide, 14 tall. compress a bit to fit colorbar
voffset = 2/16
hh, hw, vh, vw = 2/16, 4/22, 4/16, 2/22
alleyVorH = {'vertical':[3, 17, 4, 15, 6, 12, 8, 10], 'horizontal':[1, 5, 7, 2, 13, 9, 16, 14, 11]}
manualLookup = {1:[2,0],
                2:[1,0],
                3:[1,0],
                4:[1,1],
                5:[2,1],
                6:[1,2],
                7:[2,2],
                8:[1,3],
                9:[1,2],
                10:[0,3],
                11:[0,2],
                12:[0,2],
                13:[1,1],
                14:[0,1],
                15:[0,1],
                16:[0,0],
                17:[0,0]}

## End Contiguous alley plotting defaults and fx

def alleyLookup(r,c,orien):
    """
    Based on the row,col ID and whether
    we're doing a vert or horiz alley,
    find alley that should be plotted 
    """
    possibleAlleys = alleyVorH[orien]
    for alley, placement in manualLookup.items():
        if placement == [r,c] and alley in possibleAlleys:
            return alley
        


for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            print(clustname)
            unit = core.UnitData(clustname, datafile, "", alleyBounds, alleyVisits, txtVisits, position, ts)
            unit.loadData_raw()
            cn = clustname.split("\\")[0] + clustname.split("\\")[1]
            
            with PdfPages(figpath+cn+".pdf") as pdf:
                
                ###########################################
                # Whole track ratemap, no alley divisions #
                ###########################################
                
                plot.plot_2dWholeTrack(position, unit)
                try:
                    pdf.savefig()
                    print(f"Finished plot of 2d overall data")
                except:
                    print(f"Error for overall")
                plt.close()

                
                ##############################################################################
                # Stimulus Conditions A,B,C, (add overall to get overall, collapsed alleys) #
                #############################################################################
                
                mymax_txts = plot.getMaxArrays(unit, "textures") 
                mymax_overall = plot.getMaxArrays(unit, "overall")
                    
                for stim in ["A", "B", "C"]:
                    
                    if stim == "overall":
                        which_max = mymax_overall
                    else:
                        which_max = mymax_txts
                    
                    haveIM = False #Toggle for whether we've found a valid imshow we can use for the colorbar
                    fig = plt.figure()
                    for row in [0, 1]:
                        for col in [0, 1, 2, 3]:
                            left =  (hw + vw)*col
                            alley = alleyLookup(row, col, "vertical")
                            bottom = voffset + (voffset*row) + (vh*(row))
                            ax = fig.add_axes([left, bottom, vw, vh])
                            
                            if not haveIM:
                                im = plot.plot_2dAlleys(ax, unit, alley, stim, **{"max":which_max, "return":True})
                                if im:
                                    haveIM = True
                            else:
                                plot.plot_2dAlleys(ax, unit, alley, stim, **{"max":which_max})
                            ax.set_xticks([])
                            ax.set_yticks([])
                            ax.set_autoscalex_on(False)
                            ax.set_autoscaley_on(False)
                            ax.set_xmargin(0.0)
                            ax.set_ymargin(0.0)

                    for row in [0, 1, 2]:
                        for col in [0, 1, 2]:
                            left = hoffset + (hoffset*col) + (hw*(col))
                            bottom = (hh + vh)*row
                            alley = alleyLookup(row, col, "horizontal")
                            ax = fig.add_axes([left, bottom, hw, hh])
                            if not haveIM:
                                im = plot.plot_2dAlleys(ax, unit, alley, stim, **{"max":which_max, "return":True})
                                if im:
                                    haveIM = True
                            else:
                                plot.plot_2dAlleys(ax, unit, alley, stim, **{"max":which_max})
                            ax.set_xticks([])
                            ax.set_yticks([])
                            ax.set_autoscalex_on(False)
                            ax.set_autoscaley_on(False)
                            ax.set_xmargin(0.0)
                            ax.set_ymargin(0.0)

                            
                    if im:
                        plot.add_colorbar(fig, im)
                    plt.suptitle(f"{stim}", fontsize=17)
                    try:
                        pdf.savefig(dpi=200)
                        print(f"Finished plot of 2d {stim} data")
                    except:
                        print(f"Error for {stim}")
                    plt.close()
                    
                    
                ###########################################
                # Linearized Rate Maps - Common Y Scale  #
                ##########################################
                
                ylims = []
                fig, axes = plt.subplots(5,4, figsize=(10,10))
                for i in range(17):
                    ax = fig.axes[i]
                    plot.plot_linAlley(ax, unit,i+1)
                    ylims.append(ax.get_ylim()[1]) # get all max ylims of plots
                    
                # we don't know what all the ylims are until we calc and plot
                # so loop back over and adjust everything
                ymax = max(ylims)
                for i in range(17):
                    ax = fig.axes[i]
                    ax.set_ylim([0, ymax])

                pdf.savefig()
                plt.close()
                print("Finished linear RM subplots (common scale)")


                #####################################################
                # Linearized Rate Maps - Alleys Y-lim Independent   #
                #####################################################
                
                fig, axes = plt.subplots(5,4, figsize=(10,10))
                for i in range(17):
                    ax = fig.axes[i]
                    plot.plot_linAlley(ax, unit,i+1)
     
                pdf.savefig()
                plt.close()
                print("Finished linear RM subplots (independent scales)")
                
                
                ##############################
                #  Permutation Test Results  #
                ##############################
                
                Perm.unitPermutationTest_AllPairsAllAlleys(unit, 500,'', logger=True, plot='addFile')
                pdf.savefig()
                plt.close()
                print("Finished Permutation Tests")
                
                
                
                




