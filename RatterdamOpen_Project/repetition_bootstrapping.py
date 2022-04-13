# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 11:01:09 2022

@author: whockei1


Bootstrapping analysis to test putative result in "number paths x directionality" analysis
That analysis shows an inverted-U shape in directionality strength as a function of
number of paths realized through each alley. Want to know if that shape is real.
"""

import pandas as pd, matplotlib.pyplot as plt, numpy as np
import utility_fx as util 

rat_list = ['R765',
            'R765',
            'R781', 
            'R781', 
            'R808', 
            'R808', 
            'R859', 
            'R859', 
            'R886', 
            'R886']

day_list = ['RFD5',
            'DFD4',
            'D3', 
            'D4',
            'D6',
            'D7',
            'D1',
            'D2',
            'D1',
            'D2']

# rat_list = ['R781']
# day_list = ['D3']

alleydatapath = "E:\\Ratterdam\\R_data_repetition\\20220404-210901_superPopAlleyBehaviorResponse_1.5vfilt_FieldNormedTrue.csv"
alleydf = pd.read_csv(alleydatapath)

if 'Code' not in alleydf.columns:
    codes = []
    for r, row in alleydf.iterrows():
       code = f'{row["PrevDir"]}{row["CurrDir"]}{row["NextDir"]}'
       codes.append(code)
    alleydf = alleydf.assign(Code=codes)

for rat,day in zip(rat_list, day_list):
    
    rdf = alleydf[(alleydf.Rat==rat)&(alleydf.Day==day)]
 
    # interdatapath = "E:\\Ratterdam\\R_data_repetition\\20220120-164311_superPopInterBehaviorResponse_1.5vfilt.csv"
    # interdf = pd.read_csv(interdatapath)
    # alleydf = alleydf[(alleydf.Rat==r)&(alleydf.Day==d)]
    # interdf = interdf[(interdf.Rat==r)&(interdf.Day==d)]
    
    #% Calculate empirical relationship between num paths and directionality
    
    # ntrajs = []
    # diffs = []
      
    # ###### No time windowing 
    
    # traj_sample_thresh = 0
    # included_field_chunks = [] # list of tuples: [(fid,alley)]
    
    # for rname, rat in alleydf.groupby("Rat"):
    #     for daynum, day in rat.groupby("Day"):
    
    #         for alleyNum, alley in day.groupby("Alleys"):
                
    #             alley_trajs = []
                
    #             for codename, code in alley.groupby("Code"):
    #                 alley_trajs.append(code.shape[0])
                    
    #             #alley_trajs = [i for i in alley_trajs if i > traj_sample_thresh]
                
    #             alleyNumTrajs = len(alley_trajs)
    #             alleySpreadTrajs = (max(alley_trajs)-min(alley_trajs))/len(alley_trajs)
                
    #             for fid, field in alley.groupby("FieldID"):
                    
    #                 dirs  = np.unique(field.CurrDir)
                    
    #                 if len(dirs)>1:
    #                     included_field_chunks.append((np.unique(field.FieldID)[0], np.unique(field.Alleys)[0]))
    #                     diff = abs(field[field.CurrDir==dirs[0]].Rate.mean() - field[field.CurrDir==dirs[1]].Rate.mean())          
    #                     diffs.append(diff)
    #                     ntrajs.append(alleyNumTrajs)
    
    # included_field_chunks = np.asarray(included_field_chunks)      
    
    
    
    #### Time windowing 
    
    included_field_chunks = [] # list of tuples: [(fid,alley)]
    
    rdf.StartTimes = (rdf.StartTimes - rdf.StartTimes.min())/1e6 # s, ref'd to 1st sample
    
    window = 10*60
    offset = 2*60
    wins = []
    begin = 0
    stop = False
    
    while not stop:
        a,b = begin, begin + window
        if b < np.ceil(rdf.StartTimes.max()):
            wins.append((a,b))
            begin += offset
        else:
            stop = True
            
            
    numpaths_time = []
    directionality_time = []
    bias_time = []
    
    for win in wins:
        
        numpaths_win = []
        directionality_win = []
        bias_win = []
        
        windf = rdf[(rdf.StartTimes>win[0])&(rdf.StartTimes<=win[1])]    
    
        # for aNum, alley in windf.groupby("Alleys"):
        #     pathcount = 0
        #     for c, code in alley.groupby("Code"):
        #         pathcount += 1
        #     numpaths_win.append(pathcount)
            
        for fid, field in windf.groupby("FieldID"):
            
            for o, ofield in field.groupby("Orientation"):
    
                dirs = np.unique(ofield.CurrDir)
                if len(dirs)>1:
                    included_field_chunks.append((np.unique(ofield.FieldID)[0], np.unique(ofield.Alleys)[0]))
                    diff = abs(ofield[ofield.CurrDir==dirs[0]].Rate.mean()-ofield[ofield.CurrDir==dirs[1]].Rate.mean())
                    directionality_win.append(diff)
                    pathcount = 0
                    for cname, code in ofield.groupby("Code"):
                        pathcount +=1
                    numpaths_win.append(pathcount)
                    ndira = ofield[ofield.CurrDir==dirs[0]].shape[0]
                    ndirb = ofield[ofield.CurrDir==dirs[1]].shape[0]
                    bias_win.append(max(ndira,ndirb)/ofield.shape[0])
                
                    
                
        numpaths_time.append(numpaths_win)
        directionality_time.append(directionality_win)
        bias_time.append(bias_win)
    
    
    diffs  = np.asarray([i for j in directionality_time for i in j])
    ntrajs = np.asarray([i for j in numpaths_time for i in j])
    included_field_chunks = np.asarray(included_field_chunks)      
    
    #% Attempt based on fitting a quadratic to the data and bootstrapping to get many null curves
    #, the distribution of which can be compared to the empirical curve. 
    
    # first calculate # samples in each group
    ntrajdict = {i:None for i in np.unique(ntrajs)}
    
    for nt in np.unique(ntrajs):
        
        ntrajdict[nt] = len([diffs[i] for i in range(len(diffs)) if ntrajs[i]==nt])
        
    #% calculate real quadratic fit 
    npts = 50
    z =np.polyfit(ntrajs, diffs, 2)
    f = np.poly1d(z)
    xf = np.linspace(min(ntrajs),max(ntrajs),npts)
    yreal = f(xf)
     
    #% run bootstrap. code it procedurely, wrap it up if it's something I pursue
    
    nboots = 1000
    
    boot_curves = np.empty((0,npts))
    
    for nb in range(nboots):
        print(nb)
        
        shuff_diffs_sample = []
        shuff_ntrajs_sample = [] # the number of samples of each value of the trajectory 
                          # num is the same, but associated with different fields
                          
        reference_indices = np.asarray(range(included_field_chunks.shape[0])) # size of field chunks remove from this i.e w/o replacement
        
        for nt in np.unique(ntrajs):
            nsamples = ntrajdict[nt]
            idx = np.random.choice(reference_indices,nsamples,replace=False)
            
            #there has to be a better way
            for i in idx:
                argi = np.where(reference_indices==i)
                reference_indices = np.delete(reference_indices, argi)
            
            boot_fields = included_field_chunks[idx]
            
            for bf in boot_fields:
                alley = bf[1]
                fid = bf[0]
                bootdf = rdf[(rdf.Alleys==alley)&(rdf.FieldID==fid)]
                #dont need to check if num dirs > 1 bc its only included if it does
                dirs = np.unique(bootdf.CurrDir)
                bootdiff = abs(bootdf[bootdf.CurrDir==dirs[0]].Rate.mean() - bootdf[bootdf.CurrDir==dirs[1]].Rate.mean())
                shuff_diffs_sample.append(bootdiff)
                shuff_ntrajs_sample.append(nt)
                
        z = np.polyfit(shuff_ntrajs_sample, shuff_diffs_sample,2)
        f = np.poly1d(z)
        yf = f(xf)
        boot_curves = np.vstack((boot_curves, yf))

    
    #% visualize 
    
    upper = np.mean(boot_curves,axis=0)+np.std(boot_curves,axis=0)*1.96
    lower = np.mean(boot_curves,axis=0)-np.std(boot_curves,axis=0)*1.96
    
    fig, ax = plt.subplots(figsize=(17,12))
    ax.plot(xf,boot_curves.T,color='k',alpha=0.1)
    plt.plot(xf,boot_curves[0,:].T,color='k',alpha=1,label="Bootstrap Sample Quadratic") # replot one for legend labeling purposes
    ax.fill_between(xf,upper,lower,color='blue',alpha=0.7,label="97.5th Percentile Band")
    ax.plot(xf,upper,color='blue',alpha=0.7)
    ax.plot(xf,lower,color='blue',alpha=0.7)
    ax.plot(xf,yreal,color='r',linewidth=4,label="Real Quadratic",zorder=99)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel("Quadratic Fit Value of Unsigned Firing Rate Directionality",fontsize=20)
    ax.set_xlabel("Number of Realized Paths", fontsize=20)
    plt.title(f"{rat}{day} Real fitted quadratic versus 1000 bootstrap samples \n \
              {window/60}min windows, {offset/60}min offset",fontsize=25)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    lgnd = plt.legend(prop={'size':25}) 
    
    ts = util.genTimestamp()
    
    plt.savefig(f"E:\\Ratterdam\\temp\\PathTuningByLocation\\singleRecordingDays_NumPathsxDirectionality\\daybootstrap_timewindowing\\{ts}_{rat}{day}_bootstrapwindowed.png",dpi=300)
    plt.close()
