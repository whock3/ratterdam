# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:39:14 2021

@author: whockei1

Testing independence of directionality within multiple repeating fields
Primarily using two-way ANOVA from statsmodels
"""
import numpy as np, matplotlib.pyplot as plt, statsmodels.api as sm, pandas as pd
from statsmodels.formula.api import ols
import utility_fx as util
import ratterdam_RepetitionCoreFx as RepCore
import williamDefaults as wmDef 
import ratterdam_Defaults as Def
import matplotlib as mpl

savepath = 'E:\\Ratterdam\\repetition_models\\21-09-29_anovas\\'

passThresh = 3
cmap = util.makeCustomColormap()

#%% Alleys
datapath = "E:\\Ratterdam\\R_data_repetition\\20210924-145911_superPopAlleyBehaviorResponse_1.5vfilt.csv"
df = pd.read_csv(datapath)

totalAnovas = 0
fdirIntsAlley = []
fdirAlley = []
for cellid in df['CellID'].unique():
    u = df[df['CellID']==cellid]

    for orientation, directions in zip(['V', 'H'],[['N','S'],['E','W']]):
        uOriented = u[u['Orientation']==orientation]
        # get num passes by dir per field
        fielddirs = {}
        for fn, fg in uOriented.groupby('FieldNum'):
            fielddirs[fn] = {}
            for dn, dg in fg.groupby('CurrDir'):
                fielddirs[fn][dn] = dg.shape[0] 
        includedFields = []
        for fnum, dircounts in fielddirs.items():
            if directions[0] in dircounts.keys() and directions[1] in dircounts.keys():
                if dircounts[directions[0]] >= passThresh and dircounts[directions[1]] >= passThresh:
                    includedFields.append(fnum)
        if len(includedFields)>=2:
            totalAnovas += 1
            mod = ols('Rate ~ C(FieldNum)*C(CurrDir)',data=uOriented).fit()
            aov = sm.stats.anova_lm(mod,typ=2)
            print("-----------------")
            print(cellid)
            print(orientation)
            print(aov)
            if aov['PR(>F)']['C(FieldNum):C(CurrDir)'] < 0.025:
                fdirIntsAlley.append((cellid, orientation))
            if aov['PR(>F)']['C(CurrDir)'] < 0.025:
                fdirAlley.append((cellid, orientation))
                    
#Results 9/29/21 1054am
# [(2, 'H'),
#  (10, 'H'),
#  (32, 'H'),
#  (45, 'H'),
#  (82, 'H'),
#  (85, 'H'),
#  (110, 'H'),
#  (127, 'H')]

#%% Alleys: plot ratemap and field/dir responses for each unit assigned significance
allocolors = {'N':'red','S':'blue','E':'green','W':'black'}
egocolors = {'S':'red','R':'blue','B':'green','L':'black'}

for cellid in fdirIntsAlley:
    cdf = df[(df['CellID']==cellid[0])&(df['Orientation']==cellid[1])]
    fieldgroups = cdf.groupby("FieldNum")
    ncols=2
    fig, ax_ = plt.subplots(int(np.ceil(len(fieldgroups)/ncols)+1),ncols, figsize=(8,8))
    for i, (fname, fg) in enumerate(fieldgroups):
         i=i+1 # shift plots, want rm in first 
         dgroups = fg.groupby("CurrDir")['Rate'].mean()
         dgroups.plot(kind='bar',
                      yerr=fg.groupby("CurrDir")['Rate'].sem(),
                      ax=fig.axes[i],
                      ylabel="Mean Rate +/- SEM",
                      title=f"Field {fname}",
                      fontsize=8)
    plt.suptitle(f"{np.unique(cdf['CellName'])[0]} {cellid[1]} Alleys (Current Dir:Field)")
    

    clustName=np.unique(cdf['CellName'])[0]
    rat = np.unique(cdf['Rat'])[0]
    day = np.unique(cdf['Day'])[0]
        
    unit = RepCore.loadRepeatingUnit(rat, day, clustName)
    ax = fig.axes[0]
    ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                  cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
           extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
    ax.set_title(f"{unit.name}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
    ax.axis('equal')
    ax.set_ylim([0,480])
    ax.set_xlim([0, 640])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for i, field in enumerate(unit.smoothedFields):       
        ax.plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i], label=f"{i}")
    ax.legend()
    
    clustName = clustName.replace("\\","_")
    
    #rescale ylim of bar graphs, cant share y bc ratemap is on different scale
    ymax = max([ax.get_ylim() for ax in fig.axes][1:])[1]
    for ax in fig.axes[1:]:
        ax.set_ylim([0,ymax])
    
    plt.savefig(savepath+f"{rat}{day}_{clustName}_{cellid[1]}"+"AlleyFieldDirAnovas.png",dpi=300)    
    plt.close()
#%%  Intersections
datapath = "E:\\Ratterdam\\R_data_repetition\\20210922-184818_superPopInterBehaviorResponse_1.5vfilt.csv"
df = pd.read_csv(datapath)
totalAnovas = 0
fdirIntsInter = []
fdirInter = []
for cellid in df['CellID'].unique():
    u = df[df['CellID']==cellid]

    
    # get num passes by dir per field
    fielddirs = {}
    for fn, fg in u.groupby('FieldNum'):
        fielddirs[fn] = {}
        for dn, dg in fg.groupby('CurrentEgo'):
            fielddirs[fn][dn] = dg.shape[0] 
    includedFields = []
    for fnum, dircounts in fielddirs.items():
        # for intersection fields i am not currently pairing off directions (N-S,E-W)
        # and consider all possibilities valid, e.g. N above E. So just look across
        # all dirs to see if two or more have enough passes
        acceptableDirs = 0
        for _,v in dircounts.items():
            if v >= passThresh:
                acceptableDirs +=1
        if acceptableDirs >= 2:
            includedFields.append(fnum)
        
    if len(includedFields)>=2:
        totalAnovas += 1
        mod = ols('Rate ~ C(FieldNum)*C(CurrentEgo)',data=u).fit()
        aov = sm.stats.anova_lm(mod,typ=2)
        print("-----------------")
        print(cellid)
        print(aov)
        if aov['PR(>F)']['C(FieldNum):C(CurrentEgo)'] < 0.05:
            fdirIntsInter.append(cellid)
        if aov['PR(>F)']['C(CurrentEgo)'] < 0.05:
            fdirInter.append(cellid)
            

#%% Intersections: Plot ratemap and field/dir responses for each unit assigned significance

allocolors = {'N':'red','S':'blue','E':'green','W':'black'}
egocolors = {'S':'red','R':'blue','B':'green','L':'black'}

for cellid in fdirIntsInter:
    cdf = df[df['CellID']==cellid]
    fieldgroups = cdf.groupby("FieldNum")
    ncols=2
    fig, ax_ = plt.subplots(int(np.ceil(len(fieldgroups)/ncols)+1),ncols, figsize=(8,8))
    for i, (fname, fg) in enumerate(fieldgroups):
         i=i+1 # shift plots, want rm in first 
         dgroups = fg.groupby("CurrentEgo")['Rate'].mean()
         dgroups.plot(kind='bar',
                      yerr=fg.groupby("CurrentEgo")['Rate'].sem(),
                      ax=fig.axes[i],
                      ylabel="Mean Rate +/- SEM",
                      title=f"Field {fname}",
                      fontsize=8)
    plt.suptitle(f"{np.unique(cdf['CellName'])[0]} Intersections (Current Ego:Field)")
    
    #this block is bc tt7 on one day is called 7_0001 bc i changed rec params
    # fix this in raw data dir to be tt7
    rat,day,clustName=np.unique(cdf['CellName'])[0].split("_")
        
    unit = RepCore.loadRepeatingUnit(rat, day, clustName)
    ax = fig.axes[0]
    ax.imshow(unit.repUnit.rateMap2D, origin='lower', aspect='auto', interpolation='None', 
                  cmap=cmap, vmax=np.nanpercentile(unit.repUnit.rateMap2D, 98),
           extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]])
    ax.set_title(f"{unit.name}, cutoff = {round(np.nanpercentile(unit.repUnit.rateMap2D, 98),2)}Hz", fontsize=14)
    ax.axis('equal')
    ax.set_ylim([0,480])
    ax.set_xlim([0, 640])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for i, field in enumerate(unit.smoothedFields):       
        ax.plot(unit.perimeters[i][:,0], unit.perimeters[i][:,1],color=unit.colors[i], label=f"{i}")
    ax.legend()
    
    cn = clustName.replace("\\","_")
    
    #rescale ylim of bar graphs, cant share y bc ratemap is on different scale
    ymax = max([ax.get_ylim() for ax in fig.axes][1:])[1]
    for ax in fig.axes[1:]:
        ax.set_ylim([0,ymax])
    
    plt.savefig(savepath+f"{rat}{day}_{cn}"+"InterFieldDirAnovas.png",dpi=300)    
    plt.close()
    
    
#%% Making directional ratemaps using net travel (ie. schematic track) method (not pt-by-pt dir)
bins=[50,70]
savepath = 'E:\\Ratterdam\\repetition_models\\netTravelDirectionalityRMs\\'
for rat in superpop.keys():
    for day in superpop[rat].keys():
        
        refturns = superpop[rat][day]['refturns']
        ballisticTurnIdx = []
        for i,row in refturns.iterrows():
            if i < refturns.shape[0]-1:
                inter = row['Inter']
                # 9/30/21 changing turnaround exclusion for this analysis and 
                #plan to do similar for others that rely on alley+ to determine
                # the current ROI. only exclude turnaround on alley+ itself not
                # if there was a turnaround on alley-. this should still be done
                # for analysis that looks at whole turn, eg a traj decoding
                if row['Ego'] != '3' and refturns.iloc[i+1].Inter != inter:
                    ballisticTurnIdx.append(i)
                 
        turns = refturns.iloc[np.asarray(ballisticTurnIdx)]
        for unitname,unit in superpop[rat][day]['units'].items():
            print(unitname)
            
            # data for ballistic travel
            directionalRawDataB = {d:{'spikes':np.empty((0,3)), 'occs':np.empty((0,3))} for d in ['N','E','S','W']}
            directionalRMB = {d:None for d in ['N','E','S','W']}
            for direction in ['N','E','S','W']:
                for t, turn in turns.iterrows():
                    if Def.allocodedict[turn['Allo+']] == direction:
                        start, stop = float(turn['Ts exit']), float(refturns.iloc[t+1]['Ts entry'])
                        spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)]
                        occs = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=stop)]
                        
                        directionalRawDataB[direction]['spikes'] = np.vstack((directionalRawDataB[direction]['spikes'], spikes))
                        directionalRawDataB[direction]['occs'] = np.vstack((directionalRawDataB[direction]['occs'], occs))
                               
                directionalRMB[direction] = util.makeRM(directionalRawDataB[direction]['spikes'],
                                                       directionalRawDataB[direction]['occs'],bins=bins)
                
                
            # data for turnarounds
            directionalRawDataT = {d:{'spikes':np.empty((0,3)), 'occs':np.empty((0,3))} for d in ['N','E','S','W']}
            directionalRMT = {d:None for d in ['N','E','S','W']}
            for direction in ['N','E','S','W']:
                for t, turn in refturns.iterrows():
                    if Def.allocodedict[turn['Allo+']] == direction and t not in ballisticTurnIdx and t < refturns.shape[0]-1:
                        start, stop = float(turn['Ts exit']), float(refturns.iloc[t+1]['Ts entry'])
                        spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=stop)]
                        occs = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=stop)]
                        
                        directionalRawDataT[direction]['spikes'] = np.vstack((directionalRawDataT[direction]['spikes'], spikes))
                        directionalRawDataT[direction]['occs'] = np.vstack((directionalRawDataT[direction]['occs'], occs))
                               
                directionalRMT[direction] = util.makeRM(directionalRawDataT[direction]['spikes'],
                                                       directionalRawDataT[direction]['occs'],bins=bins)
            
            # this will be kept for turnaround analysis below too 
            mymax = max([np.nanpercentile(r.flatten(),98) for r in directionalRMB.values()]+[np.nanpercentile(r.flatten(),98) for r in directionalRMT.values()])
            norm = mpl.colors.Normalize(vmin=0,vmax=mymax) # for colorbar
            fig, ax = plt.subplots(4,2,figsize=(12,15))
            i=0
            for directionalRM, title in zip([directionalRMB, directionalRMT],["Ballistic","Turnaround"]):
                for _,d in enumerate(['N','E','S','W']):
                    rm = directionalRM[d]
                    fig.axes[i].imshow(rm, aspect='auto',
                                       interpolation='None',
                                       origin='lower',
                                       cmap=cmap, 
                                       vmax=mymax,
                                       extent=[wmDef.xedges[0], wmDef.xedges[-1], wmDef.yedges[0], wmDef.yedges[-1]],zorder=70)
                    fig.axes[i].set_title(f"{title} {d}")
                    util.drawTrack(rat,day,ax=fig.axes[i])
                    fig.axes[i].axis('off')
                    
                    for f, field in enumerate(unit.smoothedFields):       
                        fig.axes[i].plot(unit.perimeters[f][:,0], unit.perimeters[f][:,1],color=unit.colors[f], label=f"{f}",zorder=99)
                    i += 1
            cax = fig.add_axes([0.9,0.2,0.01,0.5])
            cbar = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
            cbar.set_label("Hz")
            cn = unit.name.replace("\\","_")
            
            
            plt.suptitle(f"{rat}{day} {unit.name} {title}")
            plt.savefig(savepath+f"{rat}{day}_{cn}_netTravelDirectionality.png",dpi=300)
            plt.close()
        
        