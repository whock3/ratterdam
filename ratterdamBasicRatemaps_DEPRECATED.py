import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket, os
from scipy.stats import sem
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys
from sys import argv
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
from ratterdam_Defaults import *
import ratterdam_visBasic as Vis

if argv and len(argv) > 1:
    exp = argv[1]
else:
    exp = "DFD4"
    
dayCode = f"R765{exp}\\"
figpath = f"C:\\Users\\whockei1\\Google Drive\\KnierimLab\\Ratterdam\\Figures\\R765{exp}\\velocityFiltered_{velocity_filter_thresh}\\"
if not os.path.isdir(figpath):
    os.mkdir(figpath)
datafile = f'E:\\Ratterdam\\R765\\R765{exp}\\'
behav = core.BehavioralData(datafile, exp, velocity_filter_thresh)
ts, position, alleyTracking, alleyVisits,  txtVisits = behav.loadData()
plot = Vis.BasicRateMaps()
gs = GS.GridSpec(5,7)
for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            print(clustname)
            unit = core.UnitData(clustname, datafile, "", alleyBounds, alleyVisits, txtVisits, position, ts)
            unit.loadData_raw()
            cn = clustname.split("\\")[0] + clustname.split("\\")[1]
            with PdfPages(figpath+cn+".pdf") as pdf:
                
                mymax = plot.getMaxArrays(unit)                    
                haveIM = False #Toggle for whether we've found a valid imshow we can use for the colorbar

                for stim in ["overall", "A", "B", "C"]:
                    fig = plt.figure(figsize=(10,10))
                    for i in range(1,18):
                        r,c = plot.gs_lookup[i]
                        ax = plt.subplot(gs[r,c])
                        if not haveIM:
                            im =  plot.plot_2dAlleys(ax, unit, i, stim, **{'max':mymax, 'return':True})

                            if im is not False:
                                haveIM = True  #So once we get a valid imshow, stop looking.
                        else:
                            plot.plot_2dAlleys(ax, unit, i, stim, **{'max':mymax})
                        ax.set_title(f"{i}")
                    if im:
                        plot.add_colorbar(fig, im)
                    plt.suptitle(f"{stim}", fontsize=17)
                    pdf.savefig()
                    plt.close()
                    print(f"Finished plot of 2d {stim} data")


                fig, axes = plt.subplots(5,4, figsize=(10,10))
                for i in range(17):
                    ax = fig.axes[i]
                    plot.plot_linAlley(ax, unit,i+1)

                pdf.savefig()
                plt.close()
                print("Finished linear RM subplots")
                