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
from scipy.stats import sem

savepath = 'E:\\Ratterdam\\repetition_models\\21-09-29_anovas\\'

passThresh = 3
cmap = util.makeCustomColormap()

#%% Alleys
datapath = "E:\\Ratterdam\\R_data_repetition\\20211003-201105_superPopAlleyBehaviorResponse_1.5vfilt.csv"
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
            uValid = uOriented[uOriented['FieldNum'].isin(includedFields)]
            totalAnovas += 1
            mod = ols('Rate ~ C(FieldNum)*C(CurrDir)',data=uValid).fit()
            aov = sm.stats.anova_lm(mod,typ=2)
            print("-----------------")
            print(cellid)
            print(orientation)
            print(aov)
            if aov['PR(>F)']['C(FieldNum):C(CurrDir)'] < 0.025:
                fdirIntsAlley.append((cellid, 
                                      orientation, 
                                      includedFields, 
                                      aov['PR(>F)']['C(FieldNum):C(CurrDir)'],
                                      aov['PR(>F)']['C(CurrDir)']
                                      ))
                
            if aov['PR(>F)']['C(CurrDir)'] < 0.025:
                fdirAlley.append((cellid, orientation))
                    

#%% Alleys: plot ratemap and field/dir responses for each unit assigned significance
allocolors = {'N':'red','S':'blue','E':'green','W':'black'}
egocolors = {'S':'red','R':'blue','B':'green','L':'black'}

for celldata in fdirIntsAlley:
    cdf = df[(df['CellID']==celldata[0])&(df['Orientation']==celldata[1])&(df['FieldNum'].isin(celldata[2]))]
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
    plt.suptitle(f"{np.unique(cdf['CellName'])[0]} {celldata[1]} Alleys (Current Dir*Field), pdir = {round(celldata[4],4)}, pint = {round(celldata[3],4)}")
    

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
    
    plt.savefig(savepath+f"{rat}{day}_{clustName}_{celldata[1]}"+"AlleyFieldDirAnovas.png",dpi=300)    
    plt.close()
#%%  Intersections
datapath = "E:\\Ratterdam\\R_data_repetition\\20211005-170829_superPopInterBehaviorResponse_1.5vfilt.csv"
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
            if dg.shape[0] >= passThresh:
                fielddirs[fn][dn] = dg.shape[0] 
    
    
    # we want to only look at the directions that are common to all fields
    common_dirs = [list(fielddirs[d].keys()) for d in fielddirs.keys()]
    common_dirs = set.intersection(*[set(i) for i in common_dirs])
    commonfielddirs = {i:{} for i in fielddirs.keys()}
    for c in common_dirs:
        for d in commonfielddirs.keys():
            commonfielddirs[d][c] = fielddirs[d][c]
            
    # for intersection fields i am not currently pairing off directions (N-S,E-W)
    # and consider all possibilities valid, e.g. N above E. So just look across
    # all dirs to see if two or more have enough passes
    includedFields = []
    for fn,fd in commonfielddirs.items():
        if len(fd) >= 2: # direction threshold - number of dirs w sufficient passes (given by passThresh) to be included
            includedFields.append(fn)
    
    #so the logic is the common_dirs part looks for directions all fields
    # share in common. And the includedFields part makes sure each field
    # has at least 2 directions. So we are looking at the common directions shared
    # across fields that have at least two of the common directions. losing data but safer/valid anova 
    
    if len(includedFields)>=2:
        uValid = u[(u['FieldNum'].isin(includedFields))&(u['CurrentEgo'].isin(common_dirs))]
        totalAnovas += 1
        mod = ols('Rate ~ C(FieldNum)*C(CurrentEgo)',data=uValid).fit()
        aov = sm.stats.anova_lm(mod,typ=2)
        print("-----------------")
        print(cellid)
        print(aov)
        if aov['PR(>F)']['C(FieldNum):C(CurrentEgo)'] < 0.05:
            fdirIntsInter.append([cellid, includedFields, common_dirs, aov['PR(>F)']['C(FieldNum):C(CurrentEgo)'], aov['PR(>F)']['C(CurrentEgo)']])
        if aov['PR(>F)']['C(CurrentEgo)'] < 0.05:
            fdirInter.append(cellid)
            

#%% Intersections: Plot ratemap and field/dir responses for each unit assigned significance

allocolors = {'N':'red','S':'blue','E':'green','W':'black'}
egocolors = {'S':'red','R':'blue','B':'green','L':'black'}

for celldata in fdirIntsInter:
    cdf = df[(df['CellID']==celldata[0])&(df['FieldNum'].isin(celldata[1]))&(df['CurrentEgo'].isin(celldata[2]))]
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
    plt.suptitle(f"{np.unique(cdf['CellName'])[0]} Intersections (Current Ego*Field), pdir={round(celldata[4])}, pint={round(celldata[3])}")
    
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
    
    clustName=np.unique(cdf['CellName'])[0]
    cn = clustName.replace("\\","_")
    
    #rescale ylim of bar graphs, cant share y bc ratemap is on different scale
    ymax = max([ax.get_ylim() for ax in fig.axes][1:])[1]
    for ax in fig.axes[1:]:
        ax.set_ylim([0,ymax])
        
    
    clustName = clustName.split("\\")[0]+"_"+clustName.split("\\")[1]
    plt.savefig(savepath+f"{cn}_"+"InterFieldDirAnovas.png",dpi=300)    
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
        

#%% Turnarounds vs Traversals
passthresh = 3
data = []
errs = []
reps = []
for cellid, cellgroup in df.groupby('CellID'):
    rep = cellgroup['Repeating'].unique()
    for fnum, fieldgroup in cellgroup.groupby("FieldID"):
        trav = fieldgroup[fieldgroup["Traversal"]==True]
        turnaround = fieldgroup[fieldgroup["Traversal"]==False]
        if trav.shape[0] >= passthresh and turnaround.shape[0] > passthresh:
            data.append([np.mean(trav['Rate']), np.mean(turnaround['Rate'])])
            errs.append([sem(trav['Rate']), sem(turnaround['Rate'])])
            reps.append(rep)
            
data = np.asarray(data)
errs = np.asarray(errs)
reps = np.asarray(['r' if i==True else 'b' for i in reps])

#%% Plot double-error scatter for turnarounds vs traversals repeating vs non
import matplotlib.patches as mpatches

plt.figure()
plt.errorbar(data[:,0], data[:,1], 
             xerr=errs[:,0], 
             yerr=errs[:,1],
             ls='',
             color=reps,
             alpha=0.5)
plt.scatter(data[:,0], data[:,1],
            c=reps,)
plt.plot([0,25],[0,25],color='k')
plt.xlabel("Mean Traversal Rate +/- SEM (Hz)",fontsize=16)
plt.ylabel("Mean Turnaround Rate +/- SEM (Hz)",fontsize=16)
plt.title("Average Rate Per Field, Turnarounds vs Traversals",fontsize=20)
repPatch = mpatches.Patch(color='red',label='Repeating Fields')
nonrepPatch = mpatches.Patch(color='blue',label='Non-repeating Fields')
plt.legend(handles=[repPatch, nonrepPatch])




#%% shuffle turnaround vs traversal
shuffdiff = []
for i in range(1000):
    print(i)
    data = []
    for cellid, cellgroup in df.groupby('CellID'):
        rep = cellgroup['Repeating'].unique()
        for fnum, fieldgroup in cellgroup.groupby("FieldID"):
            visitTypes = np.random.permutation(fieldgroup['Traversal'])
            trav = fieldgroup[visitTypes==True]
            turnaround = fieldgroup[visitTypes==False]
            if trav.shape[0] >= passthresh and turnaround.shape[0] > passthresh:
                data.append([np.mean(trav['Rate']), np.mean(turnaround['Rate'])])
               
    data = np.asarray(data)
    shuffdiff.append(np.mean(abs(data[:,0]-data[:,1])))
shuffdiff = np.asarray(shuffdiff)


#%% ANVOA turnaround vs traversal

vdf = df[df['Orientation']=='V']
mod = ols('Rate ~ Traversal*CurrDir',data=vdf).fit()
aov = sm.stats.anova_lm(mod,typ=2)
print(aov)

hdf = df[df['Orientation']=='H']
mod = ols('Rate ~ Traversal*CurrDir',data=hdf).fit()
aov = sm.stats.anova_lm(mod,typ=2)
print(aov)


#%% Want to look at COM shift for each filt. Separated into V and H 
import pickle 
from matplotlib import path
import repeatingPC as repPC
import newAlleyBounds as nab
from scipy.ndimage import center_of_mass

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
     
    
verticals = [str(i) for i in [2,3,5,7,16,14,11,9]]
horizontals = [str(i) for i in [0,4,6,1,12,8,15,13,10]]

bins=[50,70]
savepath = 'E:\\Ratterdam\\temp\\'
comshifts = []
dirs = []
for rat in superpop.keys():
    for day in superpop[rat].keys():
        ratborders = nab.loadAlleyBounds(rat,day)
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
            repeat, locCount, repeatType, overlaps = repPC.repeatingPF(unit,ratborders)

            for perim,foverlap in zip(unit.perimeters, overlaps):
                contour = path.Path(perim)
                field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                
                
                if any(str(i) in verticals for i in foverlap):
                    coms = [] # quickanddirty way is N/E is first, S/W is second, keep ordering bc want to know if COM shift (if exists) is with or against animal travel
                    for direction in ['N','S']:
                        dirspikes, dirpos = np.empty((0,3)), np.empty((0,3))

                        for t, turn in turns.iterrows():
                            if Def.allocodedict[turn['Allo+']] == direction:
                                start, stop = float(turn['Ts exit']), float(refturns.iloc[t+1]['Ts entry'])
                                spikes = field_spikes[(field_spikes[:,0]>start)&(field_spikes[:,0]<=stop)]
                                occs = field_pos[(field_pos[:,0]>start)&(field_pos[:,0]<=stop)]
                                
                                dirspikes = np.vstack((dirspikes, spikes))
                                dirpos = np.vstack((dirpos, occs))
                                
                        rm = util.makeRM(dirspikes, dirpos, bins=bins)
                        coms.append(np.asarray(center_of_mass(np.nan_to_num(rm))))
                        
                    comshifts.append(np.linalg.norm(coms[0]-coms[1]))
                    dirs.append(direction)
                    
                    
                    
                if any(str(i) in horizontals for i in foverlap):
                    coms = [] # quickanddirty way is N/E is first, S/W is second, keep ordering bc want to know if COM shift (if exists) is with or against animal travel
                    for direction in ['E','W']:
                        dirspikes, dirpos = np.empty((0,3)), np.empty((0,3))

                        for t, turn in turns.iterrows():
                            if Def.allocodedict[turn['Allo+']] == direction:
                                start, stop = float(turn['Ts exit']), float(refturns.iloc[t+1]['Ts entry'])
                                spikes = field_spikes[(field_spikes[:,0]>start)&(field_spikes[:,0]<=stop)]
                                occs = field_pos[(field_pos[:,0]>start)&(field_pos[:,0]<=stop)]
                                
                                dirspikes = np.vstack((dirspikes, spikes))
                                dirpos = np.vstack((dirpos, occs))
                                
                        rm = util.makeRM(dirspikes, dirpos, bins=bins)
                        coms.append(np.asarray(center_of_mass(np.nan_to_num(rm))))
                        
                    comshifts.append(np.linalg.norm(coms[0]-coms[1]))
                    dirs.append(direction)
                    
                    
#%% Looking at directionality on loops
# looking at directionality as a fx of loop fq, on loop vs same dir off loop
from scipy.stats import sem 

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
    
rat, day, = 'R859', 'D2'

population, refturns = superpop[rat][day]['units'], superpop[rat][day]['refturns']

onloop, offloop, opposite = [],[],[]

loop = np.array(['4','6','7','9','10','13','14','3'])
loopDirs = np.array(['E','E','S','S','W','W','N','N'])
l=0
for t, turn in refturns.iterrows():
    try:
        if turn['Alley+'] in loop:
            # get location of alley in loop
            loopLoc = np.where(loop==turn['Alley+'])[0][0]
            
            start, end  = float(turn['Ts entry']), float(refturns.iloc[t+1]['Ts exit'])
            
            #check if rat did loop. by grabbing the right number of alleys and checking if they match the loop pattern
            ratpath = refturns.iloc[t-loopLoc:t+(loop.shape[0]-loopLoc)]
            if all(ratpath['Alley+']==loop):    
                print(ratpath['Alley+'])
                for _, unit in population.items():
                    #if the cell has a field here 
                    if any([int(turn['Alley+']) in overlap for overlap in unit.overlaps]):
                        spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=end)]
                        pos = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=end)]
                        rate = spikes.shape[0]/((end-start)/1e6)
                        onloop.append(rate)
    
            # direction is same as loop but weve dropped down here so we're not on loop
            elif Def.allocodedict[turn['Allo+']] == loopDirs[loopLoc]:
                for _, unit in population.items():
                    #if the cell has a field here 
                    if any([int(turn['Alley+']) in overlap for overlap in unit.overlaps]):
                        spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=end)]
                        pos = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=end)]
                        rate = spikes.shape[0]/((end-start)/1e6)
                        offloop.append(rate)
            
            # must be other direction
            else:
                for _, unit in population.items():
                    #if the cell has a field here 
                    if any([int(turn['Alley+']) in overlap for overlap in unit.overlaps]):
                        spikes = unit.spikes[(unit.spikes[:,0]>start)&(unit.spikes[:,0]<=end)]
                        pos = unit.position[(unit.position[:,0]>start)&(unit.position[:,0]<=end)]
                        rate = spikes.shape[0]/((end-start)/1e6)
                        opposite.append(rate)
           
    except:
        pass
            
plt.bar(range(3),[np.mean(i) for i in [onloop, offloop, opposite]],yerr=[sem(i) for i in [onloop, offloop, opposite]])
        
        
#%% Popvector correlations at each alley
import numpy.ma as ma
import pickle, numpy as np, matplotlib.pyplot as plt
import newAlleyBounds as nab
import repeatingPC as repPC
from matplotlib import path

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
 
rat, day = 'R859', 'D2'
population, refturns = superpop[rat][day]['units'], superpop[rat][day]['refturns']

population, turns, refturns = superpop[rat][day]['units'], superpop[rat][day]['turns'], superpop[rat][day]['refturns']
df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
ratborders = nab.loadAlleyBounds(rat, day)
with open(df+"sessionEpochInfo.txt","r") as f:
    data = f.readlines()
session_start, session_end = [float(i) for i in data[0].split(',')]

# create time windows
seconds = (session_end - session_start)/1e6
wnSize= 2*60*1e6 # set num mins you want in a window
wnStep = 2*60*1e6
time_windows = []
stop = False
begin = session_start
while not stop:
    a,b = begin, begin + wnSize
    if b <= session_end:
        time_windows.append((a,b))
        begin += wnStep
    else:
        stop = True
time_windows = np.asarray(time_windows)

pop_alleys = {i:None for i in range(17)}
for alley in range(17):
    pop_rates = np.empty((0, time_windows.shape[0]))
    repeating = []
    includeTurnaround = []
    repeating = []
    for _,unit in population.items():
        for perim,foverlap in zip(unit.perimeters, unit.overlaps):
            if alley in foverlap:
                contour = path.Path(perim)
                field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                unit_rates = []
                for win in time_windows:
                    winStart, winEnd = win[0], win[1]
                    winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                    winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                    if winPos.shape[0] > 0:
                        winRate = winSpikes.shape[0]/((winEnd-winStart)/1e6)
                    else:
                        winRate = np.nan
                    unit_rates.append(winRate)
                
                pop_rates = np.vstack((pop_rates, unit_rates))
                corrs = np.empty((time_windows.shape[0], time_windows.shape[0]))
                # at least one day has no repeating cells
                
                for i, pi in enumerate(pop_rates.T):
                    for j,pj in enumerate(pop_rates.T):
                        corr = ma.corrcoef(ma.masked_invalid(pi), ma.masked_invalid(pj))[0,1]
                        corrs[i,j] = corr  
    pop_alleys[alley] = corrs


fig, ax = plt.subplots(3,6,figsize=(18,12))

for i,(alley, mat) in enumerate(pop_alleys.items()):
    
    if alley in verticals:
        dirA, dirB = '1', '3'
    elif alley in horizontals:
        dirA, dirB = '2', '4'
    
    dirbias = [0]
    dirtime = [time_windows[0][0]]
    for t, turn in refturns.iterrows():
        if turn['Alley+']==f"{alley}":
            if turn['Allo+'] == dirA:
                dirbias.append(dirbias[-1]+1)
                dirtime.append(float(turn['Ts entry']))

            elif turn['Allo+'] == dirB:
                dirbias.append(dirbias[-1]-1)
                dirtime.append(float(turn['Ts entry']))
    dirtime = np.asarray(dirtime)
    dirbias = np.asarray(dirbias)
    
    
    fig.axes[i].imshow(mat, aspect='auto',interpolation='None',origin='lower',extent=[time_windows[0][0], time_windows[-1][1], 0, mat.shape[0]])
    ax2 = fig.axes[i].twinx()
    ax2.plot(dirtime,dirbias,color='k')
    fig.axes[i].set_title(alley)
plt.suptitle(f"{rat}{day}, {(wnSize/1e6)/60}min windows, {(wnStep/1e6)/60}min step size")
    
#%% Pop vector corr between blocks

import numpy.ma as ma, pickle 
import newAlleyBounds as nab
import repeatingPC as repPC
import matplotlib.path as path

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]
blocks = [[2,0,3,1],
          [3,4,5,12],
          [5,6,7,8],
          [16,1,14,15],
          [14,12,11,13],
          [11,8,9,10]]

with open("E:\\Ratterdam\\R_data_repetition\\superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
    
rat, day = 'R886', 'D2'
population, turns, refturns = superpop[rat][day]['units'], superpop[rat][day]['turns'], superpop[rat][day]['refturns']
df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
ratborders = nab.loadAlleyBounds(rat, day)
with open(df+"sessionEpochInfo.txt","r") as f:
    data = f.readlines()
session_start, session_end = [float(i) for i in data[0].split(',')]

# create time windows
seconds = (session_end - session_start)/1e6
wnSize= 10*60*1e6 # set num mins you want in a window
wnStep = 1*60*1e6
time_windows = []
stop = False
begin = session_start
while not stop:
    a,b = begin, begin + wnSize
    if b <= session_end:
        time_windows.append((a,b))
        begin += wnStep
    else:
        stop = True
time_windows = np.asarray(time_windows)

pop_alleys = {i:None for i in range(6)}

for orientation, alleys in zip(range(6), blocks):
    pop_rates = np.empty((0, time_windows.shape[0]))
    repeating = []
    includeTurnaround = []
    repeating = []
    for _,unit in population.items():
        for perim,foverlap in zip(unit.perimeters, unit.overlaps):
            if any(a in alleys for a in foverlap):
                repeating.append(unit.repeating)
                contour = path.Path(perim)
                field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                unit_rates = []
                for win in time_windows:
                    winStart, winEnd = win[0], win[1]
                    winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                    winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                    if winPos.shape[0] > 0:
                        winRate = winSpikes.shape[0]/((winEnd-winStart)/1e6)
                    else:
                        winRate = np.nan
                    unit_rates.append(winRate)
                
                pop_rates = np.vstack((pop_rates, unit_rates))
                corrs = np.empty((time_windows.shape[0], time_windows.shape[0]))
                # at least one day has no repeating cells
                
                for i, pi in enumerate(pop_rates.T):
                    for j,pj in enumerate(pop_rates.T):
                        corr = ma.corrcoef(ma.masked_invalid(pi), ma.masked_invalid(pj))[0,1]
                        corrs[i,j] = corr  

    pop_alleys[orientation] = corrs
    
fig, ax = plt.subplots(2,3, figsize=(15,12))
for i, block in enumerate(blocks):
    fig.axes[i].set_title(f"Block {i}")
    fig.axes[i].imshow(pop_alleys[i], aspect='auto', interpolation='None')
plt.suptitle(f"{rat}{day} Track block-based population time correlation ({(wnSize/1e6)/60}min, {(wnStep/1e6)/60}step)")


#%% Pop vector corr between V H,

import numpy.ma as ma, pickle 
import newAlleyBounds as nab
import repeatingPC as repPC
import matplotlib.path as path

verticals = [2,3,5,7,16,14,11,9]
horizontals = [0,4,6,1,12,8,15,13,10]
perimeter = [0,4,6,7,9,15,13,10,2,16]
interior = [3,14,5,11,12,1,8]
blocks = [[2,0,3,1],
          [3,4,5,12],
          [5,6,7,8],
          [16,1,14,15],
          [14,12,11,13],
          [11,8,9,10]]

with open("E:\\Ratterdam\\R_data_repetition\\21-10-19_superPopulationRepetition.pickle","rb") as f:
    superpop = pickle.load(f)      
    
rat, day = 'R886', 'D2'
population, turns, refturns = superpop[rat][day]['units'], superpop[rat][day]['turns'], superpop[rat][day]['refturns']
df = f'E:\\Ratterdam\\{rat}\\{rat}_RatterdamOpen_{day}\\'
ratborders = nab.loadAlleyBounds(rat, day)
with open(df+"sessionEpochInfo.txt","r") as f:
    data = f.readlines()
session_start, session_end = [float(i) for i in data[0].split(',')]

# create time windows
seconds = (session_end - session_start)/1e6
wnSize= 5*60*1e6 # set num mins you want in a window
wnStep = 1*60*1e6
time_windows = []
stop = False
begin = session_start
while not stop:
    a,b = begin, begin + wnSize
    if b <= session_end:
        time_windows.append((a,b))
        begin += wnStep
    else:
        stop = True
time_windows = np.asarray(time_windows)

pop_alleys = {i:None for i in range(6)}

for orientation, alleys in zip(range(2), [verticals, horizontals]):
    pop_rates = np.empty((0, time_windows.shape[0]))
    repeating = []
    includeTurnaround = []
    repeating = []
    for _,unit in population.items():
        for perim,foverlap in zip(unit.perimeters, unit.overlaps):
            if any(a in alleys for a in foverlap):
                repeating.append(unit.repeating)
                contour = path.Path(perim)
                field_pos = unit.position[contour.contains_points(unit.position[:,1:])]
                field_spikes = unit.spikes[contour.contains_points(unit.spikes[:,1:])] 
                unit_rates = []
                for win in time_windows:
                    winStart, winEnd = win[0], win[1]
                    winSpikes = field_spikes[(field_spikes[:,0]>winStart)&(field_spikes[:,0]<=winEnd)]
                    winPos = field_pos[(field_pos[:,0]>winStart)&(field_pos[:,0]<=winEnd)]
                    if winPos.shape[0] > 0:
                        winRate = winSpikes.shape[0]/((winEnd-winStart)/1e6)
                    else:
                        winRate = np.nan
                    unit_rates.append(winRate)
                
                pop_rates = np.vstack((pop_rates, unit_rates))
                corrs = np.empty((time_windows.shape[0], time_windows.shape[0]))
                # at least one day has no repeating cells
                
                for i, pi in enumerate(pop_rates.T):
                    for j,pj in enumerate(pop_rates.T):
                        corr = ma.corrcoef(ma.masked_invalid(pi), ma.masked_invalid(pj))[0,1]
                        corrs[i,j] = corr  

    pop_alleys[orientation] = corrs
    
fig, ax = plt.subplots(2,1, figsize=(15,12))
for i, orien in enumerate(['Verticals','Horizontals']):
    fig.axes[i].set_title(orien)
    fig.axes[i].imshow(pop_alleys[i], aspect='auto', interpolation='None')
plt.suptitle(f"{rat}{day} Verticals vs Horizontals time correlation ({(wnSize/1e6)/60}min, {(wnStep/1e6)/60}step)")