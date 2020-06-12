"""
Script to pre-process behavior. Defines alley visits and the texture present during those visits.
First, aligns timestamps from arduino/python (ratterdam itself) and the corresponding Neuralynx ts
"""
import numpy as np
import csv
from matplotlib import pyplot as plt
import ratterdam_Defaults as Def
import  datetime, os, pickle
import sys

sys.path.insert(0, 'E:\\UserData\\Documents\\GitHub\\ratterdam\\')
import statistics as stats
from difflib import SequenceMatcher
import utility_fx as util
from bisect import bisect_left
import matplotlib.path as mplPath
import ratterdam_Directionality as Dir
import ratterdam_Defaults as Def
import ratterdam_DataFiltering as Filt

def parseEvents(p, stimData, alleySwapTS, expCode):
    alleyTracking = compute_alleyTracking(p)
    alleyVisits = {i:[] for i in range(17)} # each list will contain ts corresponding to time of alley visit beginning
    txtVisits = {i:[] for i in range(17)}
    for alley,timestamps in alleyTracking.items():
        if ("BR" in expCode or "BS" in expCode) and alley not in [i-1 for i in [1,3,5,7,8,10,11,16,17]]:
            pass
        else:
            vts = getVisits(timestamps[:,0],Def.default_visit_gap)
            vts = [i for i in vts if len(i)>1]
            vts = [i for i in vts if (i[-1]-i[0])>Def.default_visit_duration]
            visitTs = [(i[0], i[-1]) for i in vts] # vts returns list of list, each of which is all ts in a group. So first ts is the first ts of a visit, grab it
            alleyVisits[alley] = visitTs
            for visit in visitTs:
                intervalIdx = checkInterval(visit[0], alleySwapTS[alley])
                txtVisits[alley].append(stimData['stimuli'][alley][intervalIdx])
    return alleyTracking, alleyVisits, txtVisits


def parseEvents_Beltway_oldCopy(p, stimData, alleySwapTS, expCode):
    alleyTracking = compute_alleyTracking(p)
    alleyVisits = {j:[] for j in range(17)} # each list will contain ts corresponding to time of alley visit beginning
    txtVisits = {j:[] for j in range(17)}
    for lap in range(len(alleySwapTS[0])-1):
        begin, end = alleySwapTS[0][lap], alleySwapTS[0][lap+1]
        for alley, timestamps in alleyTracking.items():
            if alley in [i-1 for i in [1,3,5,7,8,10,11,16,17]]:
                posInAlley = timestamps[np.logical_and(timestamps[:,0] >= begin,timestamps[:,0] <= end)]
                if posInAlley.shape[0]>1:
                    alleyVisits[alley].append((posInAlley[0,0], posInAlley[-1,0]))
                    txtVisits[alley].append(stimData['stimuli'][alley][lap])
                else:
                    print(f"alley {alley}, {lap} lap")
    return alleyTracking, alleyVisits, txtVisits


def parseEvents_Beltway(p, stimData, alleySwapTS, expCode):
    alleyTracking = compute_alleyTracking(p)
    alleyVisits = {j:[] for j in range(17)} # each list will contain ts corresponding to time of alley visit beginning
    txtVisits = {j:[] for j in range(17)}
    for lap in range(len(alleySwapTS[0])-1):
        begin, end = alleySwapTS[0][lap], alleySwapTS[0][lap+1]
        for alley in Def.beltwayAlleys:
            alleyVisits[alley-1].append((begin,end))
            txtVisits[alley-1].append(stimData['stimuli'][alley-1][lap])
    return alleyTracking, alleyVisits, txtVisits
        
    


def parseEvents_new(position, stimData, alleySwapTS, epsilon=2):
    """
    Create data structures with visit position data and matching txt info
    Epsilon is distance in camera coordinate values. I.e. within epsilon you're
    close enough to lickport/alley border to qualify as a visit.
    """
    alleyTracking = compute_alleyTracking(position)
    alleyVisits = {i:[] for i in range(17)} # each list will contain ts corresponding to time of alley visit beginning
    txtVisits = {i:[] for i in range(17)}
    for alley,timestamps in alleyTracking.items():
        visitIdx,visits = getVisits(timestamps[:,0],Def.default_visit_gap)
        bounds = Dir.extractCorners(Def.alleyBounds[alley])
        ul,ll,ur,lr = bounds
        validVisits = []
        for entry,data in zip(visitIdx, visits):
            potVisit = alley[entry]
            mins = [abs(mp-i) for i in v[:,1]]
            if min(mins) < epsilon:
                validVisits.append(data)
        visitTs = [(i[0], i[1]) for i in validVisits]
        for visit in visitTs:
            intervalIdx = checkInterval(visit[0]/1e6, alleySwapTS[alley])
            txtVisits[alley].append(stimData['stimuli'][alley][intervalIdx])
    return alleyTracking, alleyVisits, txtVisits


def compute_alleyTracking(pos, alleyShape='rectangle'):
    alleyTracking = {i:np.empty((0,3)) for i in range(17)}
    
    if alleyShape == 'rectangle':   
        for entry in pos:
            for a in Def.alleyBounds.keys():
                if entry[1] > Def.alleyBounds[a][0][0] and entry[1] < Def.alleyBounds[a][0][1] and entry[2] > Def.alleyBounds[a][1][0] and entry[2] < Def.alleyBounds[a][1][1]:
                    alleyTracking[a] = np.vstack((alleyTracking[a], entry))
                    
    elif alleyShape == 'polygon':
        # mpl has a function to define a polygon as a series
        # of vertices. Then the path has a method to check if
        # given points are inside. Gives bool mask. So filter position by it
        
        for alley, vertices in Def.alleyBounds.items():
            polygon = mplPath.Path(vertices)
            pointsInAlley = polygon.contains_points(position[:,1:])
            alleyTracking[alley] = Def.position[pointsInAlley]
            
    return alleyTracking
            
        

def checkInterval(x, intervals):
    if type(intervals) == list:
        size = len(intervals)-1
    elif type(intervals) == np.ndarray:
        size = intervals.shape[0]-1
    target = None
    for i in range(size):
        if intervals[i] <= x <= intervals[i+1]:
            target = i
    return target

def getNLSwapTimes_RF(rec):
    swapTimes = []
    for i in range(len(rec)-4):
        if rec[i][5] == 64 and rec[i+1][5] == 0 and rec[i+2][5] == 64 and rec[i+3][5] == 0:
            swapTimes.append(rec[i][3])
    swapTimes = [i/1e6 for i in swapTimes]
    swapclean = []
    for i in range(len(swapTimes)-1):
        if swapTimes[i+1] - swapTimes[i] < 1:
            swapclean.append(swapTimes[i])
    return swapclean


def getNLSwapTimes_DF(rec):
    swapTimes = []
    for i in range(len(rec)-4):
        if rec[i][5] == 1 and rec[i+1][5] == 0:
            swapTimes.append(rec[i][3])
    swapTimes = [i/1e6 for i in swapTimes]
    return swapTimes


def getTTLCode(rec):
    ttls = []
    for entry in rec:
      if b'TTL' in entry[10]:
        ttls.append(str(entry[10]).split(" ")[-1][5:7])
    ttls = [int(i) for i in ttls]
    ttls = [i for i in ttls if i !=0]
    x = f'0x00{stats.mode(ttls)}'
    return bytes(x, "UTF-8")

def getNLSwapTimes_Beltway(rec, expCode):

    swapTimes = []
    if expCode == 'BRD6':
        swapTimes.append(2484730000)

    ttlCode = getTTLCode(rec)

    if "BR" in expCode:
        targetBegin = b'begin session again clockwise'
        targetEnd = b'end session'
    elif "BS" in expCode:
        targetBegin = b'begin session counterclockwise'
        targetEnd = b'end session counterclockwise'


    for i,entry in enumerate(rec):
        if similar(entry[10], targetBegin) > 0.95:
            nlBegin = i
        elif similar(entry[10], targetEnd) > 0.95:
            nlEnd = i

    recChunk = [i for i in rec[nlBegin+1:nlEnd]]

    for i,entry in enumerate(recChunk):
        if entry[5] == int(ttlCode, 16) and recChunk[i+1][5] == int(b'0x0000',16):
            swapTimes.append(entry[3])
            
    swapTimes.append(rec[nlEnd][3])

    #swapTimes = [i/1e6 for i in swapTimes]

    return  swapTimes


def beltway_getLapTimes(datafile, stimData, expCode):
    hdr, rec = util.readNEV(datafile)
    swapTimes = getNLSwapTimes_Beltway(rec, expCode)
    alleySwapTS = convertPytoNLts_Beltway(swapTimes, stimData)
    return alleySwapTS


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def convertPytoNLts_OLDDFKEEPTOCHECKDIFFs(rec, nltimes, pydata):
    '''nl times cleaned first''' 
    pySwapTimes = [i[0] for i in pydata['order'][17:]]
    adjPy = [i - pySwapTimes[0] for i in pySwapTimes]
    projs = [min(enumerate(nltimes), key=lambda x: abs(x[1]-i)) for i in [j+nltimes[0] for j in adjPy]]
    swapTS = [nltimes[i[0]] for i in projs]
    
    
    beginTS = [i[3] for i in rec if b'begin session' in i[10]][0]/1e6
    endTS = [i[3] for i in rec if b'end session' in i[10]][0]/1e6
    
    
    alleySwapTS = {i:[] for i in range(17)}
    for i in range(len(swapTS)):
        alleySwapTS[pydata['order'][i+17][1]].append(swapTS[i])
    for i in range(17):
        alleySwapTS[i] =  [beginTS] + alleySwapTS[i] + [endTS]
    return alleySwapTS

def convertPytoNLts(rec, nltimes, pydata):
    '''nl times cleaned first''' 
    pySwapTimes = [i[0] for i in pydata['order'][17:]]
    adjPy = [i - pySwapTimes[0] for i in pySwapTimes]
    projs = [min(enumerate(nltimes), key=lambda x: abs(x[1]-i)) for i in [j+nltimes[0] for j in adjPy]]
    swapTS = [nltimes[i[0]] for i in projs]
    
    beginTS = [i[3] for i in rec if similar(b'begin session',i[10]) > 0.85][0]
    endTS = [i[3] for i in rec if similar(b'end session', i[10]) > 0.85][0]
    
    
    
    #beginTS = [i[3] for i in rec if b'begin session' in  i[10]][0]/1e6
    #endTS = [i[3] for i in rec if b'end session' in i[10]][0]/1e6
    
    
    alleySwapTS = {i:[] for i in range(17)}
    for i in range(len(swapTS)):
        alleySwapTS[pydata['order'][i+17][1]].append(swapTS[i])
    for i in range(17):
        alleySwapTS[i] = [beginTS] + alleySwapTS[i] + [endTS]
    return alleySwapTS

def convertPytoNLts_Beltway(swapTimes, stimData):
    """
    DEPRECATED FOR BELTWAY. DOESNT ALWAYS WORK
    """
#    swapTimes = [i/1e6 for i in swapTimes]
#    diffs = [i - stimData['lapStarts'][0] for i in stimData['lapStarts']]
#    projs = [min(enumerate(swapTimes), key=lambda x: abs(x[1]-i)) for i in [j+swapTimes[0] for j in diffs]]
#    swapTS = [swapTimes[i[0]]*1e6 for i in projs]
    
    swapTS = []
    
    
    alleySwapTS = {i:[] for i in range(17)}
    
    for i, swap  in enumerate(swapTS):
        for alley in alleySwapTS.keys():
            alleySwapTS[alley].append(swap)
            
        
    return alleySwapTS

def getBeltwayLapTimes(datafile):
    """
    This version is current as of 9/19. All other fxs related
    to syncing NL and PY ts for beltway are considered deprecated
    
    Read in txt file called "sessionLapInfo.txt" Must be this name.
    
    Reads in a txt file that contains all laptimes
    The file should include a ts that marks the end of the last lap. I.e. if
    there are 60 trials there should be 61 entries. This last entry can be 
    the end session ts or the ts of another lap youve queued up (but obviously did not run)
    for the purposes of using its start ts
    """
    
    with open(datafile + 'sessionLapInfo.txt','r') as file:
        lines = file.readlines()
        
    swapTS = [int(i.rstrip().rstrip(',')) for i in lines if '#' not in i]
    
    alleySwapTS = {i:[] for i in range(17)}
    
    for i, swap  in enumerate(swapTS):
        for alley in alleySwapTS.keys():
            alleySwapTS[alley].append(swap)
            
        
    return alleySwapTS
    
    
      

def truncateWithinSession(position,ts, alleySwapTS):
    """this is for beltway so assuming all alleys have same ts bc its
    synchronous updating"""
    beginTS, endTS = alleySwapTS[0][0], alleySwapTS[0][-1]
    position = position[np.logical_and(position[:,0] >= beginTS, position[:,0] < endTS)]
    ts = ts[np.logical_and(ts >= beginTS, ts < endTS)]
    return position, ts
      

def getVisits(data, maxgap=3*1e6):
    """
    """
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


#def getVisits(data, maxgap=3*1e6):
#    data.sort()
#    groups = [[data[0]]]
#    idx = [[0]]
#    for i,x in enumerate(data[1:]):
#        if abs(x - groups[-1][-1]) <= maxgap:
#            groups[-1].append(x)
#            idx[-1].append(i)
#        else:
#            groups.append([x])
#            idx.append([i])
#            
#    for v in range(len(idx)):
#        idx[v] = np.asarray(idx[v])
#            
#    return idx,groups



def getStimFileName(datafile):
    fs = []
    for f in os.listdir(datafile):
        # as of 3/19 _Stimuli.p is a pickled dict with fields for
        # stimuli order and alley swap order. _BehavioralDaya is for
        # beltway and has rows for alley rewards [0,1,0,0...] and txts [A,B,C,B,C...]
        # and arow for lap ts
        if f.endswith('_Stimuli.p') or f.endswith('_BehavioralData.csv'):
            fs.append(f)

    return fs
 
def extractRows(reader):
    rows = []
    for row in reader:
        rows.append(row)
    return rows
       

def loadBeltwayData(basedir,stimFiles, expCode):
    """assumes a standard data file layout"""
    beltwayAlleys = [i-1 for i in [16,17,3,1,5,7,8,10,11]]
    #find file matching the exp you want
    
    for file in stimFiles:
        with open(basedir+file,"r") as csvfile:
            reader = csv.reader(csvfile)
            rows = extractRows(reader)
            #naming convention is counterclockwise ccw and clockwise cw in file
            # but for expcode want 2 letter code so BS is ccw and BR cw using chemistry R/S nomenclature for rotations
            if (rows[0][1] == "CCW" and "BS" in expCode) or (rows[0][1] == "CW" and "BR" in expCode):
                targetFile = file

    with open(basedir+targetFile, "r") as dataFile:
        reader = csv.reader(dataFile)
        stimData = {"stimuli":{alley:[] for alley in beltwayAlleys}, "rewards":{alley:[] for alley in beltwayAlleys}}
        rows = extractRows(reader)
        
        for i in range(len(rows)):
            rows[i] = [x for x in rows[i] if x != '']
            
            
        sessionLen = len(rows[3])
            
        stimData["lapStarts"] = [float(i) for i in rows[3][1:sessionLen]]
        for x,alley in enumerate(beltwayAlleys):
            
            stimData['stimuli'][alley] = rows[4+(x*2)+1][1:sessionLen] # here you take sessionlen-1 because last entry is the lap cued up after last lap when lap over button was pressed ending session
            stimData['rewards'][alley] = [int(i) for i in rows[4+(x*2)][1:sessionLen]]
            

    return stimData
        

def adjustPosCamera(datafile, pos, ts):
     """
     Camera can be oriented differently w.r.t track across days/rats
     if things get rotated. This reads in a .txt file that should be 
     in each day's data dir saying whether x or y should be flipped
     format is e.g. x:640\ny:None
     """
     with open(datafile+"cameraOrientationInfo.txt","r") as cam:
         info = cam.readlines()
     info = [i.rstrip() for i in info]
     info = {info[0][0]:info[0][2:], info[1][0]:info[1][2:]}
     if info['x'] != 'None':
         posx = [int(info['x'])-pos[i][0] for i in ts]
     else:
         posx = [pos[i][0] for i in ts]
    
     if info['y'] != 'None':
         posy = [int(info['y'])-pos[i][1] for i in ts]
     else:
         posy = [pos[i][1] for i in ts]

     return posx, posy
   
def getDaysBehavioralData(datafile, expCode):
    pos = util.read_pos(datafile)
    hdr, rec = util.readNEV(datafile)
    
    ts = np.asarray(sorted(list(pos.keys())))
    posx, posy = adjustPosCamera(datafile, pos, ts)
    position = np.column_stack((ts,posx,posy))

    stimFile = getStimFileName(datafile)
    if "Stimuli" in stimFile:
        stimData = pickle.load(open(datafile+stimFile,"rb")) # single file bc ratterdam full is usually 1 sess
    elif "BehavioralData" in stimFile[0]:
        stimData = loadBeltwayData(datafile, stimFile, expCode) # stimFile is a list for beltway,list of sessions in a day
            
#    if "RF" in expCode:
#        swapTimes = getNLSwapTimes_RF(rec)
#    elif "DF" in expCode:
#        swapTimes = getNLSwapTimes_DF(rec)
#    elif "BR" in expCode or "BS" in expCode:
#        swapTimes = getNLSwapTimes_Beltway(rec, expCode)
        
    if "RF" in expCode or "DF" in expCode:
        alleySwapTS = convertPytoNLts(rec, swapTimes, stimData)
        alleyTracking, alleyVisits,  txtVisits = parseEvents(position, stimData, alleySwapTS, expCode)
        
    elif "BR" in expCode or "BS" in expCode:
        alleySwapTS = getBeltwayLapTimes(datafile)
        
        p_sess,ts_sess = truncateWithinSession(position, ts, alleySwapTS)
        if Def.velocity_filter_thresh > 0:
            p_sess = Filt.velocity_filtering(p_sess)
        alleyTracking, alleyVisits,  txtVisits = parseEvents_Beltway(p_sess, stimData, alleySwapTS, expCode)
        
   
    return alleyTracking, alleyVisits,  txtVisits, p_sess, ts_sess

