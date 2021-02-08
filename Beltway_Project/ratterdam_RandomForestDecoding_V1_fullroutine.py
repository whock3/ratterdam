# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:30:14 2019

@author: whockei1
"""

############
# Imports  #
############

import sklearn as skl
import matplotlib.pyplot as plt
from sklearn import svm, preprocessing, metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

import numpy as np, socket, os, sys,json
from datetime import datetime
from scipy.stats import sem
import matplotlib.colors as colors
from importlib import reload

if socket.gethostname() == 'Tolman':
    codeDirBase = 'C:\\Users\\whockei1\\Google Drive'
elif socket.gethostname() == 'DESKTOP-BECTOJ9':
    codeDirBase = 'C:\\Users\\whock\\Google Drive'
    
sys.path.insert(0, codeDirBase + '\\KnierimLab\\Ratterdam\\Code')
import utility_fx as util
import ratterdam_ParseBehavior as Parse
import ratterdam_CoreDataStructures as Core
import ratterdam_PermutationTests as Perm
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt

from itertools import cycle


#####################
# Helper Functions  #
#####################
def compute_epoch(val,size):
    """hardcode that session is divided
    into thirds. find in which third the trial is in"""
    propthrusess = val/size
    if propthrusess < 0.25:
        epoch = '0'
    elif propthrusess <= 0.5:
        epoch = '1'
    elif propthrusess <= 0.75:
        epoch = '2'
    elif propthrusess < 1.:
        epoch = '3'
    return epoch

def checkRM(ratemap):
    """Utility function to take a 1-d linear ratemap
    and see if it is valid.
    May 2019: it's not empty ie. there's data
    and the nanmax of that data exceeds firing rate
    thresh defined locally"""
    if type(ratemap) == np.ndarray and np.nanmax(ratemap) >= frThresh:
        return True
    else:
        return False
    
def generateLabel(target, alley, stimulus, epoch):
    if target == 'Alley':
        label = str(alley)
    elif target == 'Stimulus':
        label = stimulus
    elif target == 'Epoch':
        label = epoch
    elif target == 'AlleyXStimulus':
        label = f"{alley}{stimulus}"
    elif target == 'AlleyXEpoch':
        label  = f"{alley}{epoch}"
    elif target == 'StimulusXEpoch':
        label = f"{stimulus}{epoch}"
    elif target == 'AlleyXStimulusXEpoch':
        label = f"{alley}{stimulus}{epoch}"
    
    return label


target_classes = {'Alley':[str(i) for i in [16,17,3,1,5,7,8,10,11]],
                  'Stimulus': ['A','B','C'],
                  'AlleyXStimulus': ['10A', '10B', '10C', '11A', '11B', '11C', '16A', '16B', '16C',
                                    '17A', '17B', '17C', '1A', '1B', '1C', '3A', '3B', '3C', '5A', '5B',
                                   '5C', '7A', '7B', '7C', '8A', '8B', '8C']
                 }


###############
# Load Data  #
###############

datafile = "E:\\Ratterdam\\R808\\R808_Beltway_D7\\"
rat = 'R808'
expCode = "BRD7"
figurePath = 'E:\\Ratterdam\\multidayFigures\\randomForest\\'
alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess = Parse.getDaysBehavioralData(datafile, expCode)
population = {}
for subdir, dirs, fs in os.walk(datafile):
    for f in fs:
        if 'cl-maze1' in f and 'OLD' not in f and 'Undefined' not in f:
            clustname = subdir[subdir.index("TT"):] + "\\" + f
            unit = Core.UnitData(clustname, datafile, expCode, Def.alleyBounds, alleyVisits, txtVisits, p_sess, ts_sess)
            unit.loadData_raw()
            if Filt.checkMiniumUnitActivity(unit, alleyVisits, threshold=50): # threshold refers to number of spikes emitted on track outside of start box
                print(clustname)
                population[unit.name] = unit
                
                
###################
# Set parameters  #
###################
                
frThresh = 0 #measured in Hz. Pick something close to 0, or 0 itself. This applies to individual visits to alleys, linearized ratemap
#target = 'AlleyXStimulus' #choices are Alley, Texture, Epoch, or some 2- or 3-member combination of these
targetList = ['Alley', 'Stimulus', 'AlleyXStimulus']
beltwayAlleys = [16, 17, 3, 1, 5, 7, 8, 10, 11] # beltway alley IDs in terms of their full track, 17-alley ID
nbins = Def.singleAlleyBins[0]-1
avgType = 'macro' # for signal detection / performance metrics which are not inherently multiclass (e.g. all but accuracy), pick how to aggregate individual class results
nRuns = 500
#shuffleOrNot = True
split_size = 0.75 # defined in terms of train size, proportion 0-1

fig, axes = plt.subplots(2,3,figsize=(8,6))

for ax, (shuffleOrNot, target) in enumerate(zip([False, False, False, True, True, True], targetList*2)):
    
#########################
# Create data matrices  #
#########################

    X = np.empty((0, nbins))
    Y = []
    
    for alley in beltwayAlleys:
        visitSize = len(population[list(population.keys())[0]].alleys[alley]) # all units have same behavioral data obviously so use first unit by default to get num viists to alley.
        for visitNum in range(round(visitSize/2), visitSize): 
            #dataRow = np.empty((0))
            epoch = compute_epoch(visitNum, visitSize)
            stimulus = population[list(population.keys())[0]].alleys[alley][visitNum]['metadata']['stimulus'] #again, stims are same for all units so use first unit to grab it
            label =  generateLabel(target, alley, stimulus, epoch)
            
            invalidRow = 0 # initialize to valid, set to invalid upon finding an invalid rm
            for unitname, Unit in population.items():
                dataRow = np.empty((0))
                try:
                    rm = Unit.alleys[alley][visitNum]['ratemap1d']
                    if checkRM(rm)== True:
                        dataRow = np.concatenate((dataRow, rm))
                        Y.append(label)
                        X = np.vstack((X, dataRow))
    
                    else:
                        dataRow = np.concatenate((dataRow, np.zeros((nbins)) ))
                except:
                    pass
                    
            # Y.append(label)
            # X = np.vstack((X, dataRow))
    
    X[np.where(~np.isfinite(X))] = 0
    X = preprocessing.StandardScaler().fit_transform(X)
    
    
    if shuffleOrNot is True:
        if target == 'AlleyXStimulus':
            txt = [i[-1] for i in Y]
            np.random.shuffle(txt)
            a = [i[:-1] for i in Y]
            Y = [f"{x}{y}" for x,y in zip(a,txt)]
        elif target == 'Stimulus':
            np.random.shuffle(Y)
        elif target == 'Alley':
            np.random.shuffle(Y)
            
            
    ##################################
    # Run random forest classifier   #
    ##################################
            
    oobs, precisions, recalls, f1s, accuracies = [], [], [], [],[]
    
    print(f"{expCode}, {target}")
    for i in range(nRuns):
        print(i)
        clf = RandomForestClassifier(n_estimators=700, oob_score=True)
        Xtrain, Xtest, ytrain, ytest = train_test_split(X,Y,shuffle=True,random_state=0)
        clf.fit(Xtrain,ytrain)
        oobs.append(clf.oob_score_)
        yfit = clf.predict(Xtest)
        p = precision_score(ytest, yfit, average=avgType)
        r = recall_score(ytest, yfit, average=avgType)
        f1 = f1_score(ytest, yfit, average=avgType)
        acc = accuracy_score(ytest,yfit)
        precisions.append(p)
        recalls.append(r)
        f1s.append(f1)
        accuracies.append(acc)
        
    data = {'oobs':oobs, 'precisions':precisions, 'accuracies':accuracies, 'recalls':recalls, 'f1s':f1s}
    ts =  "{:%Y%m%d_%H%M}".format(datetime.now())
    with open(figurePath+f"RFdecoding_{target}_{rat}{expCode}_{ts}.json", 'w') as outfile:
        json.dump(data, outfile)
    
    
    
    ################################
    # Plot Data and save as a pdf  #
    ################################
        
    fig.axes[ax].hist(precisions, color='b', alpha=0.5)
    fig.axes[ax].hist(recalls, color='r', alpha=0.5)
    fig.axes[ax].hist(f1s, color='g', alpha=0.5)
    fig.axes[ax].hist(accuracies, color='k', alpha=0.5)
    fig.axes[ax].hist(oobs,color='purple',alpha=0.5)
    #fig.axes[ax].legend(["Precision", "Recall", "F1 Score", "Accuracy", "OOB"])
    #plt.vlines(weightedChance(.56,.3),0,plt.ylim()[1])
    fig.axes[ax].set_ylabel("Frquency")
    fig.axes[ax].set_xlabel("Performance")
    fig.axes[ax].set_title(f"{target}, {(lambda x: 'Real' if x == False else 'Shuffle')(shuffleOrNot)}",fontsize=12)
# loop to adjust xlims for each pair of graphs for a condition, real and shuffle
for ax_ in [[0, 3],[1,4],[2,5]]:
    mymin = min([fig.axes[ax_[0]].get_xlim()[0], fig.axes[ax_[1]].get_xlim()[0]])
    mymax = max([fig.axes[ax_[0]].get_xlim()[1], fig.axes[ax_[1]].get_xlim()[1]])
    fig.axes[ax_[0]].set_xlim([mymin, mymax])
    fig.axes[ax_[1]].set_xlim([mymin, mymax])


plt.suptitle(f"Random Forest Decoding for {rat} {expCode}. {nRuns} Runs")
plt.tight_layout()
plt.savefig(figurePath+f"{rat}{expCode}_{ts}.png", format='png', dpi=300)
