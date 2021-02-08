# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 13:40:53 2019

@author: whockei1

Ratterdam Beltway
Visualization Script #1
Plot linear ratemaps as matrix of heatmaps, divided by txt
"""
import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import numpy as np
from scipy.stats import sem
import utility_fx as util
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import ratterdam_visBasic as Vis
from scipy.interpolate import interp1d


rat = "R808"
expCode = "BRD7"


datafile = f"E:\\Ratterdam\\{rat}\\{rat}{expCode}\\"
figpath = f"E:\\Ratterdam\\{rat}\\beltway_plots\\{expCode}\\"
beltwayAlleys = [16,17,3,1,5,7,8,10,11]
date = util.genTimestamp()

ratplot = Vis.BasicRateMaps()

alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)

stimfiles = Parse.getStimFileName(datafile)
stimData = Parse.loadBeltwayData(datafile,stimfiles,expCode)

cmap = util.makeCustomColormap()

page_layout = gridspec.GridSpec(4,6)

alley_name_conversion = {16:1, 17:2, 3:3, 1:4, 5:5, 7:6, 8:7, 10:8, 11:9}

for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            print(clustname)
            if True:  # this is here so you can replace w a given tt and not have to change indents 
                unit = Core.UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
                unit.loadData_raw()
                cn = clustname.split("\\")[0] + clustname.split("\\")[1]
                
                with PdfPages(figpath+date+"_"+unit.name+f"_smooth{Def.smoothing_1d_sigma}_{Def.singleAlleyBins[0]-1}bins_{Def.velocity_filter_thresh}vfilt_R{Def.includeRewards}"+".pdf") as pdf:
    
                    for alley in beltwayAlleys:
                        
                        fig = plt.figure(figsize=(14,16), dpi=300)
                        
                        #overall ratemap
                        ax1 = fig.add_subplot(page_layout[0,:2])
                        
                        # linear plots by txts
                        ax2 = fig.add_subplot(page_layout[0,2:])
                        
                        # stack of linear ratemaps
                        ax3 = fig.add_subplot(page_layout[1,:2])
                        ax4 = fig.add_subplot(page_layout[1,2:4])
                        ax5 = fig.add_subplot(page_layout[1,4:])
                        
                         # rasters space
                        ax6 = fig.add_subplot(page_layout[2,:2])
                        ax7 = fig.add_subplot(page_layout[2,2:4])
                        ax8 = fig.add_subplot(page_layout[2,4:])
                        
                        # rasters time
                        # Am dividing each texture's raster plot into two to create
                        # a broken x axis to show short time responses (variable maybe to say 3s)
                        # and then break then everything after. Otherwise early stuff is often squished on long trials
                        ax9 = fig.add_subplot(page_layout[3,0])  # A short time response
                        ax10 = fig.add_subplot(page_layout[3,1]) # A longer time
                        ax11 = fig.add_subplot(page_layout[3,2]) # B short
                        ax12 = fig.add_subplot(page_layout[3,3]) # B long
                        ax13 = fig.add_subplot(page_layout[3,4]) # C short 
                        ax14 = fig.add_subplot(page_layout[3,5]) # C long
                        axTimeRasterA = [ax9, ax10]
                        axTimeRasterB = [ax11, ax12]
                        axTimeRasterC = [ax13, ax14]

                        #################################
                        #    Linear average ratemaps    #
                        #################################
                    
                        visits = {txt:np.empty((0, Def.singleAlleyBins[0]-1)) for txt in ['A','B','C']}
                        visitsSpk = {txt:[] for txt in ['A','B','C']}
                        visitIdx = {txt:[] for txt in ['A','B','C']}
                        visitTickColor = {txt:[] for txt in ['A','B','C']}
    
                        
                        for i,visit in enumerate(unit.alleys[alley]):
                            txt  = visit['metadata']['stimulus']
                            visits[txt] = np.vstack((visits[txt], visit['ratemap1d']))
                            visitsSpk[txt].append(visit['spikes'])
                            if Def.includeRewards == 2:
                                visitIdx[txt].append(i)
                            elif Def.includeRewards == 0 or Def.includeRewards == 1:
                                visitIdx[txt].append(visit['metadata']['nrtrialnum'])
                            
                            
                        
                        ###############################
                        #   Linear Ratemap Stacks     #
                        ###############################
                        
                        #max_for_imshow = np.percentile([np.nanmax(visits[txt]) for txt in ['A','B','C']],pct)
    
                        flatalldata = np.ndarray.flatten(np.vstack((visits['A'],visits['B'],visits['C'])))
                        max_for_imshow = np.nanpercentile(flatalldata, Def.singlealley_imshow_pct_cutoff)
                        
                        for ax,txt in zip([ax3, ax4, ax5], ['A','B','C']):
                            ax.imshow(visits[txt], cmap=cmap,vmax=max_for_imshow,interpolation='None',aspect='auto',origin='lower')
                            ax.set_yticks(range(len(visitIdx[txt])))
                            ax.set_yticklabels(visitIdx[txt])
                            
                            ax2.tick_params(axis='x', labelsize=14)
                            ax2.tick_params(axis='y', labelsize=14)
                            
                            
                            #for ticklabel, tickcolor in zip(ax.get_yticklabels(), visitTickColor[txt]):
                             #   ticklabel.set_color(tickcolor)
                            
                            if txt == 'A':
                                ax.set_title(f"{txt}, color saturates at {Def.singlealley_imshow_pct_cutoff}% {round(max_for_imshow,2)}Hz", fontsize=16)
                            else:
                                ax.set_title(txt, fontsize=16)
                            
                        for c,data in zip(['r','b','g'],[visits['A'],visits['B'],visits['C']]):
                            
                            mask = np.ma.masked_invalid(data)
                            avg = mask.mean(axis=0) # ignores inf and nan
                            err = np.std(mask,axis=0)/np.sqrt(np.sum(~mask.mask,axis=0))
                            
                            ax2.plot(avg,color=c)
                            ax2.fill_between(range(len(avg)), avg+err, avg-err, color=c,alpha=0.5)
                            
                            ax2.tick_params(axis='x', labelsize=14)
                            ax2.tick_params(axis='y', labelsize=14)
                            
                            
                        ##################################
                        # Overall ratemap - whole track  #
                        ##################################  
                        
                        rbins, cbins = Def.wholeAlleyBins
                        rows, cols = np.linspace(0, 480, rbins), np.linspace(0,640, cbins)
                        hs,xe,ye = np.histogram2d(unit.spikes[:,2],unit.spikes[:,1],bins=[rows, cols])
                        ho = np.histogram2d(unit.position[:,2],unit.position[:,1],bins=[rows, cols])[0]
                        n = (hs*np.reciprocal(ho))*30
                        n[np.where(ho==0)] = np.nan
                        n = util.weird_smooth(n,Def.smoothing_2d_sigma)
                        n[np.where(ho==0)] = np.nan
                        
                        #todo - make sep rms for start box and track
                        wholermmax = np.nanpercentile(n,Def.wholetrack_imshow_pct_cutoff)
                        ax1.imshow(n, aspect='auto', interpolation='none', cmap=cmap, origin='lower',
                                   extent=[xe[0],xe[-1],ye[0],ye[-1]], vmax = wholermmax)
                        ax1.set_title(f"{Def.wholetrack_imshow_pct_cutoff}% Percentile FR: {round(wholermmax,2)}")
                        
    #                    for i in beltwayAlleys:
    #                        a = Def.alleyBounds[i-1]
    #                        x1,y1,x2,y2 = a[0][0], a[1][0], a[0][1], a[1][1]
    #                        for x,y in zip([[x1, x1], [x1, x2], [x2, x2], [x1, x2]], [[y1, y2], [y2, y2], [y1, y2], [y1, y1]]):
    #                            ax1.plot(x, y, 'k')
                        
                        
                        #############################
                        #      Rasters  (Spatial)   #
                        #############################
                            
                        numP = 1000
                        ul, ll, ur, lr = util.extractCorners(Def.alleyBounds[alley-1])
                        if alley in [1,5,7,16,11]:
                            mpA = util.midpoint(ll,ul)
                            mpB = util.midpoint(lr,ur)
                            segmentX = np.linspace(mpA[0], mpB[0], numP)
                            segmentY = [mpA[1]]*numP
                            segment = list(zip(segmentX,segmentY))
                    
                        elif alley in [3,17,8,10]:
                            mpA = util.midpoint(ll,lr)
                            mpB = util.midpoint(ul,ur)
                            segmentX = [mpA[0]]*numP
                            segmentY = np.linspace(mpA[1], mpB[1], numP)
                        
                        segment = list(zip(segmentX,segmentY))
                            
                        
                        
                        projP = {txt:[] for txt in ['A', 'B', 'C']}
                        for txt in ['A', 'B', 'C']:
                            for i in range(len(visitsSpk[txt])):
                                if visitsSpk[txt][i].shape[0] > 0:
                                    projP[txt].append(util.project_points(segment, visitsSpk[txt][i][:,1:]))
                                else:
                                    projP[txt].append([])
                                    
                        for ax,txt in zip([ax6, ax7, ax8],['A','B','C']):
                            for i in range(len(projP[txt])):
                                if type(projP[txt][i]) == np.ndarray and projP[txt][i].shape[0] >0:
                                    ax.scatter(projP[txt][i], [i]*projP[txt][i].shape[0], s=6, c='k',alpha=0.5)
                                    ax.set_xlim([0,numP])
                            ax.set_yticks(range(len(visitIdx[txt])))
                            ax.set_yticklabels(visitIdx[txt])
                            
                        
                        ########################
                        # Rasters (Temporal)   #
                        ########################
#                        def findTrialStart(realIdx):
#                            for trial in unit.alleys[alley]:
#                                if trial['metadata']['nrtrialnum'] == realIdx:
#                                    return trial['occs'][0,0]
#                        
#                        timeCutoff = 2*1e3 # cutoff time in trial to divide spikes between subrasters (long, short timescale responses)
#                        
#                        for txt, ax in zip(['A','B','C'], [axTimeRasterA, axTimeRasterB, axTimeRasterC]):
#                            for i,(trial,idx) in enumerate(zip(visitsSpk[txt],visitIdx[txt])):
#                                trialStart= findTrialStart(idx)
#                                spikesZeroed = (trial[:,0] - trialStart)/1e3 # zero spike array 'trial' to trial start and convert from us to ms
#                                ax[0].plot(spikesZeroed, [i]*spikesZeroed.shape[0], marker='|', color='k',linestyle='')
#                                ax[1].plot(spikesZeroed, [i]*spikesZeroed.shape[0], marker='|', color='k',linestyle='')
#                            ax[0].set_xlim(left=0,right=timeCutoff)
#                            ax[1].set_xlim(left=timeCutoff)
#                            ax9.set_title("Time Raster Zeroed to trial start (ms)",fontsize=10)
#                            ax[0].set_yticks(range(len(visitIdx[txt])))
#                            ax[0].set_yticklabels(visitIdx[txt])
#                            
#                            ax[0].spines['right'].set_visible(False)
#                            ax[1].spines['left'].set_visible(False)
#                            ax[1].set_yticklabels([])
#                            ax[1].set_yticks([])
                                                    


                            
                        plt.suptitle(f"Alley {alley_name_conversion[alley]}")
                        try:
                            pdf.savefig()
                        except:
                            print(f"Error with {unit.name} alley {alley_name_conversion[alley]}")
                        plt.close()
                    