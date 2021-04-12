# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 14:28:25 2018

@author: whockei1
"""

import numpy as np, matplotlib.pyplot as plt, random, json, pickle, datetime, copy, socket
from scipy.stats import sem
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter as gauss # for smoothing ratemaps
import sys
from bisect import bisect_left

import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_DataFiltering as Filt
import ratterdam_Defaults as Def



class UnitData():
    ''' 
    Class to store data for a single unit. 
    
    Data divided by alley in a dict.
    
    Data is either loaded from raw cl-files using behavioral
    parsing functions or it loads saved data from jsons depending
    on loadType value.
    
    Needed in namespace: 
    - alleyTracking, txtVisits, alleyVisits (if loading raw)
    - datafile
    - position
    
    '''
    def __init__(self, unitName, datafile, experimentCode, alleyBounds, alleyVisits, txtVisits, position, ts, loadType='raw'):
        self.name = unitName.split("\\")[0]+unitName.split("\\")[1]
        self.namePath = unitName
        self.datafile = datafile
        self.experimentCode = experimentCode
        self.alleys = {i:[] for i in range(1,18)}
        self.alleyBounds = alleyBounds
        self.alleyVisits = alleyVisits
        self.txtVisits = txtVisits
        self.position = position
        self.ts = ts
        self.alleyTracking = Parse.compute_alleyTracking(position)
        
        
    def unitVelocityFilter(self, clust):
        '''
        Remove spikes emitted when the animal's
        velocity was below threshold. 
        
        This class does not explicitly store
        the thresh, or even whether a velocity
        filter was applied. That occurs in the BehavioralData
        class and is reflected implicitly in the position array
        having been filtered. the ts dict is unfiltered and is used
        as a reference.
        '''
        
        allSpikeTs = np.asarray([util.takeClosest(self.ts, i) for i in clust])
        filtTs = clust[np.isin(allSpikeTs, self.position[:,0])]
        return filtTs
        
    def loadData_raw(self):
        '''Loads from parseBehavioral helper fx.
        Outputs of which must be in namespace'''
        self.computeAlleyBins()
        clust = util.read_clust(self.datafile+self.namePath)
        clust = np.asarray(clust)
        clust = self.unitVelocityFilter(clust)
        spikexy = util.getPosFromTs(clust,self.position)
        spikes = np.column_stack((clust,spikexy))
        #spikes = spikes[np.where(spikes[:,2] > 50)]
        self.spikes = spikes
        self.linRMS = {i:{'A':None, 'B':None, 'C':None} for i in range(1,18)}
        self.alleyRMS = {i:{'A':None, 'B':None, 'C':None} for i in range(1,18)} # ie. 2d RM
        
        stimfiles = Parse.getStimFileName(self.datafile)
        stimData = Parse.loadBeltwayData(self.datafile, stimfiles, self.experimentCode)
        
        
        for alley in Def.beltwayAlleys:
            for visit in range(len(self.alleyVisits[alley-1])):
                #visitsOcc = util.getVisitPos(alley-1, visit, self.ts, self.alleyVisits, self.position)
                #visitsSpk = util.getVisitPos(alley-1, visit, self.spikes[:,0], self.alleyVisits, self.position)
                
                visitsOcc = self.position[(self.position[:,0]>self.alleyVisits[alley-1][visit][0])&(self.position[:,0]<=self.alleyVisits[alley-1][visit][1])]
                visitsSpk = self.spikes[(self.spikes[:,0]>self.alleyVisits[alley-1][visit][0])&(self.spikes[:,0]<=self.alleyVisits[alley-1][visit][1])]
                
                #The above gets all spikes and occs within a lap. WH EDit 5/29/20 and see Analysis&Code notebook I for more info
                # Now I will use the txt file with lap starts I create to set lap time bounds and get EVERYTHING in that for all alleys
                # Here below you apply alleybounds a a cookie cutter to punch out only whats in the alley
                # Previously I had, for legacy reasons, more complicated approach that took those lap times then
                # found all spikes/occs in alleytracking that were in that window and took alley start and stop to be first and last occs there
                # but it turned out to be not perfect
                visitsOcc = util.checkInAlley(visitsOcc, alley)
                visitsSpk = util.checkInAlley(visitsSpk, alley)
                
                # If flag includeRewards==False, check if trial is one that was
                # rewarded and exclude it. This is to prevent reward-related
                # firing rate confounds. 
                nrtrialnum = None
           
                interval = Parse.checkInterval(visitsOcc[0,0],[i[0] for i in self.alleyVisits[alley-1]]+[self.alleyVisits[alley-1][-1][1]])
                reward = stimData['rewards'][alley-1][interval]

                if Def.includeRewards == 0:
                    if reward == 1:
                        include = False
                    else:
                        include = True
                        nrtrialnum = visit
                elif Def.includeRewards == 1:
                    if reward == 1:
                        include = True
                        nrtrialnum = visit
                    else:
                        include = False
                    
                elif Def.includeRewards == 2:
                    include = True
                    nrtrialnum = visit

                                  
                if include == True:
                    self.alleys[alley].append({'spikes': visitsSpk, 
                                           'occs': visitsOcc, 
                                           'ratemap2d': self.computeSingleRM(visitsSpk, visitsOcc, alley, dim=2),
                                           'ratemap1d': self.computeSingleRM(visitsSpk, visitsOcc, alley, dim=1),
                                           'metadata': self.generateVisitMetadata(alley, visit, nrtrialnum, reward)
                                            })
                
        allS, allO = np.empty((0,3)), np.empty((0,3)) # for collecting all spikes/occs to make overall RM
        for alley in range(1,18):
            for txt in ['A', 'B', 'C']:
                s,p = self.collapseVisits(alley, txt)
                if len(s) > 0 and len(p) > 0:
                    s,p = util.list_to_arrays(s),util.list_to_arrays(p)
                    allS, allO = np.vstack((allS, s)), np.vstack((allO, p))
                    collapsed1DRateMap = self.computeSingleRM(s, p, alley,dim=1)
                    collapsed2DRateMap = self.computeSingleRM(s, p, alley,dim=2)
                    self.linRMS[alley][txt] = collapsed1DRateMap
                    self.alleyRMS[alley][txt] = collapsed2DRateMap
            self.alleyRMS[alley]['overall'] = self.computeSingleRM(allS, allO, alley,dim=2)
                
    def computeAlleyBins(self):
        longDimBins, shortDimBins = Def.singleAlleyBins
        self.alleyBins = {i:{'rows':None,'cols':None} for i in range(17)}
        for i,v in enumerate(self.alleyBounds.values()):
            x,y = v
            if i in [i-1 for i in [1,5,7,2,13,9,16,14,11]]:
                bins = [shortDimBins, longDimBins] #again, np.2dhist takes [x,y] which means [c, r]
            elif i in [i-1 for i in [3,4,6,8,17,15,12,10]]:
                bins = [longDimBins, shortDimBins]
            else:
                print("error")
            self.alleyBins[i]['rows'] = np.linspace(self.alleyBounds[i][1][0], self.alleyBounds[i][1][1],num=bins[0])
            self.alleyBins[i]['cols'] = np.linspace(self.alleyBounds[i][0][0], self.alleyBounds[i][0][1],num=bins[1])
            
    
                
    def computeSingleRM(self, spikes,position,alley,dim):
        '''Compute a single rate map from an array of spikes
        and pos of array form (ts,x,y). Not using track coords, just 
        0-max bin size in either array. Chooses which dim is long
        depending on which alley (lookup in fx)'''
        #So np.hist2d takes bins in [x,y] which is [c,r].
        # this is opposite of the [r,c] convntion thts more common
        rbins,cbins = self.alleyBins[alley-1]['rows'], self.alleyBins[alley-1]['cols']
        hs = np.histogram2d(spikes[:,2],spikes[:,1],bins=[rbins, cbins])[0]
        ho = np.histogram2d(position[:,2],position[:,1],bins=[rbins, cbins])[0]
        if dim == 2:
            n = (hs*np.reciprocal(ho))*30
            n[np.where(ho==0)] = np.nan
            n = util.weird_smooth(n,Def.smoothing_2d_sigma)
            n[np.where(ho==0)] = np.nan
        elif dim == 1:
            ls,lo  = np.sum(hs,axis=util.getAxType(ho)), np.sum(ho,axis=util.getAxType(ho))
            n = (ls* np.reciprocal(lo)) * 30
            if np.count_nonzero(~np.isnan(n))>1:
                n = util.stepsmooth(n,Def.smoothing_1d_sigma)
            n[np.where(lo==0)] = np.nan
        return n
    
    def collapseVisits(self, alley, txt):
        visitSpikes, visitPos = [], []
        for visit in self.alleys[alley]:
            if visit['metadata']['stimulus'] == txt:
                visitSpikes.append(visit['spikes'])
                visitPos.append(visit['occs'])
                
        return visitSpikes, visitPos
    
    def generateVisitMetadata(self, alley, visit, nrtrialnum, reward):
        """
        Generates metadata for each trial. Metadata is a dict so can
        easily be extended as more info needs to be included
        
        txt - stimulus present on that trial at that alley
        nrtrialnum - the overall trial number. Used when excluding reward so one has a record of what 
                    the actual trial num is as opposed to just numbering the surviving trials in their relative order
        reward - 0,1 if reward was not given, or given respectively on trial. This added 10/10/20. 
        """
        txt = self.txtVisits[alley-1][visit]
        return {"stimulus":txt, "nrtrialnum":nrtrialnum, "reward":reward}
        
    
    
    
class BehavioralData():
    """
    Class to read in position data from pos.p
    Also reads in behavior about alley visits and stimuli present.
    Does not store, just returns. Wrapper fx loadData.
    """
    def __init__(self,  datafile, experimentCode, vfilt):
        self.datafile = datafile
        self.experimentCode = experimentCode
        self.vfilt = vfilt # value of velocity filter threshold
                            # if 0, don't filter. Else filter using thresh given
        
    def loadData(self):
        """
        Returns ts, position, alleyTracking, alleyVisits,  txtVisits
        """
        pos = util.read_pos(self.datafile)
        ts = np.asarray(sorted(list(pos.keys())))
        posx = [640 - pos[i][0] for i in ts]
        posy = [pos[i][1] for i in ts]
        for t in ts:  
            pos[t][0] = 640 - pos[t][0]
        position = np.column_stack((ts,posx,posy))
        position = position[np.where(position[:,2] > 50)]
        if self.vfilt > 0:
            position = Filt.velocity_filtering(position)
        alleyTracking, alleyVisits,  txtVisits = Parse.getDaysBehavioralData(self.datafile, self.experimentCode)
        return ts, position, alleyTracking, alleyVisits,  txtVisits
    
    
if __name__ == '__main__':
    
    rat = "R859"
    expCode = "BRD3"
    datafile = f"E:\\Ratterdam\\{rat}\\{rat}{expCode}\\"
    clustname = 'TT8\\cl-maze1.4'
      
    print(f"Beginning {rat} {expCode} {clustname}")
    alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
    unit = UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
    unit.loadData_raw()
    print(f"{rat} {expCode} {clustname} loaded")
