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
    totalpaths_time = []
    
    for win in wins:
        
        numpaths_win = []
        directionality_win = []
        bias_win = []
        totalpaths_win = []

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
                    numpaths_win.append(ofield.shape[0]) # numpaths is the behavioral variable we're using. in other script it stores pathcount. here its # samples (ofield.shape)
                    totalpaths_win.append(ofield.shape[0]) # unused
                    ndira = ofield[ofield.CurrDir==dirs[0]].shape[0]
                    ndirb = ofield[ofield.CurrDir==dirs[1]].shape[0]
                    bias_win.append(max(ndira,ndirb)/ofield.shape[0])
                
                
        numpaths_time.append(numpaths_win)
        totalpaths_time.append(totalpaths_win)
        directionality_time.append(directionality_win)
        bias_time.append(bias_win)
    
    # So the data above were calculated within time windows and now 
    # we are pooling all that data together. 
    diffs  = np.asarray([i for j in directionality_time for i in j])
    ntrajs = np.asarray([i for j in numpaths_time for i in j])
    totaltrajs = np.asarray([i for j in totalpaths_time for i in j])
    included_field_chunks = np.asarray(included_field_chunks)      
    
    #% Attempt based on fitting a quadratic to the data and bootstrapping to get many null curves
    #, the distribution of which can be compared to the empirical curve. 
    
    # first calculate # samples in each group
    ntrajdict = {i:None for i in np.unique(totaltrajs)}
    
    for nt in np.unique(ntrajs):
        
        ntrajdict[nt] = len([diffs[i] for i in range(len(diffs)) if ntrajs[i]==nt])
        
    #% calculate real quadratic fit 
    npts = 50
    z =np.polyfit(totaltrajs, diffs, 2)
    f = np.poly1d(z)
    xf = np.linspace(min(totaltrajs),max(totaltrajs),npts)
    yreal = f(xf)
     
    #% run bootstrap. code it procedurely, wrap it up if it's something I pursue
    
    nboots = 1000
    
    boot_curves = np.empty((0,npts))
    
    for nb in range(nboots):
        print(nb)
        
        shuff_diffs_sample = []
        shuff_ntrajs_sample = [] # the number of samples of each value of the trajectory 
                          # num is the same, but associated with different fields
                          
        reference_indices = np.asarray(range(totaltrajs.shape[0])) # size of field chunks remove from this i.e w/o replacement
        
        for nt in np.unique(totaltrajs):
            nsamples = ntrajdict[nt]
            idx = np.random.choice(reference_indices,nsamples,replace=False)
            
            #there has to be a better way
            for i in idx:
                argi = np.where(reference_indices==i)
                reference_indices = np.delete(reference_indices, argi)
            
            boot_fields = diffs[idx]
            
            for bf in boot_fields:
                shuff_diffs_sample.append(bf)
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
    ax.set_xlabel("Number of Straight-Through Passes Through Field", fontsize=20)
    plt.title(f"{rat}{day} Real fitted quadratic versus 1000 bootstrap samples \n \
              {window/60}min windows, {offset/60}min offset",fontsize=25)
    
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    lgnd = plt.legend(prop={'size':25}) 
    
    ts = util.genTimestamp()
    
    plt.savefig(f"E:\\Ratterdam\\temp\\PathTuningByLocation\\reassignment_behaviorVersusDirectionality\\{ts}_{rat}{day}_totalPathsWindowedReassignment.png",dpi=300)
    plt.close()
