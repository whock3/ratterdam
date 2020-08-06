import os, json, pickle, numpy as np, re, math, socket, scipy as sp
from numpy.linalg import norm as npNorm
import scipy.ndimage
from scipy.interpolate import interp1d
from bisect import bisect_left
import scipy.ndimage
from matplotlib.colors import LinearSegmentedColormap
import ratterdam_CoreDataStructures as Core
import ratterdam_ParseBehavior as Parse
import ratterdam_Defaults as Def



def readNEV(datafile):
    "Reads a Neuralynx Events File (.nev) and returns the header and events separately"
    NEVfile = datafile + 'Events.nev'
    f = open(NEVfile,'rb')
    hdr = f.read(16*1024)
    dt =  np.dtype([('nstx',np.int16),('npkt_id',np.int16),('npkt_data_size',np.int16),('qwTimeStamp',np.uint64),('nevent_id',np.int16),('nttl',np.int16),('ncrc',np.int16),('ndummy1',np.int16),
               ('ndummy2',np.int16),('dnExtra',np.int32,8),('EventString',np.string_,128)])
    rec = np.fromfile(f,dtype=dt)
    f.close()
    return hdr,rec

def read_events(datafile):
    '''read in Events.nev, parse to TTL ts, and link those to behavioral events
    
    THis is only for the T-maze data collected for R735
    '''

    #get NL ts for each event
    hdr, rec = readNEV(datafile)
    ttl_code = re.compile(b'0x0000') #pulses come in pairs where one is 0x0000 and other has a nonzero in last pos. I'm grabbing that by convention.
    # for 171104 the nonzero pulse is 0x0010 which breaks my pattern so grabbing the zero one. should be a neglible timing diff. fix if needed. 
    pulsetimes = []
    for event in rec:
        hit = re.search(ttl_code, event[-1])
        if hit is not None:
            pulsetimes.append(event)
    pulsetimes = [i[3] for i in pulsetimes] # redefine to get NL ts

    #pair NL ts with behavioral events
    for f in os.listdir(datafile):
        if f.endswith('_parsedData.json'):
            events = json.load(open(datafile+f,'r'))

    triallen = len(events) - 3  # because metadata,rat,time entries arent trial data
    trial_pulses = [pulsetimes[x:x + 3] for x in range(0, len(pulsetimes), 3)] #chunk into 3 events per trial - begin the ITI, gate open, and decision (agnostic as to correct/incorr)
    for t in range(triallen):
        if t == triallen-1:
            events[str(t)]['Begin_NLts'] = trial_pulses[t][0] #add NL ts as separate entry. keep Ard ts

        else:
            events[str(t)]['Begin_NLts'] = trial_pulses[t][0]  # add NL ts as separate entry. keep Ard ts
            events[str(t)]['GateOpen_NLts'] = trial_pulses[t][1]
            events[str(t)]['Outcome_NLts'] = trial_pulses[t][2]

    return events


def load_data(datafile):
    """
    This is for T-maze, not ratterdam
    """
    os.chdir(datafile)
    events = read_events(datafile)
    feats = Feats.FeatureBankCreator(events)
    feats.createFeatIdx()
    feats.createOutcomeIdx()
    sessionData = pickle.load(open(datafile + 'sessionDataV2.p', 'rb'))
    user = os.getlogin()
    date_interm = os.getcwd().split('\\')[-1].split('_')[0].split('-')  # just splits into [y,m,d]
    d = date_interm[0][2:] + date_interm[1] + date_interm[2]
    return sessionData, feats, events 

def hitTrials(feature, label, features,  labels):
    if label == 'L' or label == 'R':
        idx = np.hstack((np.where(labels == label+'-C'), np.where(labels == label + '-I')))
        return np.intersect1d(features[feature], idx)
    else:
        return np.intersect1d(features[feature], np.where(labels == label))
    
def extract(tt, rmType, sessionData):
    rms = np.asarray([sessionData[tt][i][rmType] for i in range(len(sessionData[tt].keys()))])
    if 'ratemap' in rmType.lower():
        rms = rms.reshape(rms.shape[0],rms.shape[2]) #rms are (len(events),len(track)). 2 dim 
    else:
        if rms.size != 0:
            rms = rms.reshape(rms.shape[0])
    return rms

def read_clust(clustfile):
    '''open a cl-mazeX.X file and read spike times into a list'''
    with open(clustfile) as clust:
        raw_clust_data = clust.read().splitlines()
    spikes = []
    for line in raw_clust_data[13:]:
        spikes.append(float(line.split(',')[-1]))
    return spikes



def read_pos(datafile):
    '''open pos.p.ascii and parse into {ts:[x,y,dir]} '''
    posfile = datafile + 'Pos.p.ascii'
    with open(posfile) as posP:
        raw_position_data = posP.read().splitlines()
        on_data = posP.read().splitlines()
    position_data = {}
    for line in raw_position_data[25:]:
      data = [float(a) for a in line.split(',')]
      position_data[data[0]] = data[1:]
    return position_data


def load_position(datafile):
    """ Returns np array ts,x,y. Not to be confused with
    util.read_pos which just loads the raw dict"""
    pos = read_pos(datafile)
    
    ts = np.asarray(sorted(list(pos.keys())))
    posx = [pos[i][0] for i in ts]
    posy = [480 - pos[i][1] for i in ts]
    position = np.column_stack((ts,posx,posy))
    
    return position


def makeLinTrack():
    
    fstem = interp1d([630,189],[215,215])
    fleft = interp1d([220,208],[415,215])
    fright = interp1d([195,208],[15,215])

    xstem = np.linspace(215,630,200) # DID DIST of stem and L, R basically double it for stem
    xleft = np.linspace(208,220,100)
    xright = np.linspace(195,208,100)

    ystem = fstem(xstem)
    yleft = fleft(xleft)
    yright = fright(xright)
    track = list(zip(reversed(xstem),reversed(ystem))) + \
            list(zip((xleft),(yleft))) + \
            list(zip(reversed(xright),reversed(yright)))

    return track


def project_points(track, points):
    projected = []
    for p in points:
        if p[0] != 0 and p[1] != 0:
            projp = 0
            dist = 10000 # essentially infinite here
            for i,refp in enumerate(track):
                d = math.hypot(refp[0]-p[0],refp[1]-p[1])
                if d < dist:
                    dist = d
                    projp = i
            projected.append(projp)
    if projected == []:
        projected = np.empty((0))
    projected = np.asarray(projected)
    return projected

def smooth(array,smthType='gaussian',smthWin=10):
    '''Smooth using a filter specified by arg. Note the 
    choice of window depends on what filter.'''
    if smthType == 'gaussian':
        return scipy.ndimage.gaussian_filter1d(array,smthWin)
    elif smthType == 'windowed':
         np.asarray([np.mean(array[:,0+(5*i):5+(5*i)]) for i in range(int(array.shape[1]/5))])

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def computeSummaryTrace(trials,traceType='avg'):
    '''
    Take array of trials grouped by txt response in some way
    (e.g. paths to B on R) and compute a 'summary trace' where
    that is def as a single track-len trace +/- error.
    types: "avg" =  smoothed trials, then avg of that +/i sem of that
    '''
    traceTypes = ['avg']
    if traceType not in traceTypes:
        return 'Invalid Summary Trace Type'
    
    if traceType == 'avg':
        trials = smooth(trials)
        avg = np.mean(trials,axis=0)
        sem = scipy.stats.sem(trials)
        return avg, sem
        

def testStat_OvAll(test, other1, other2):
    '''test statistic that takes each summary
    trace (e.g. avg) in turn and computes diff
    of it vs others. return all three in A,B,H order'''
    testStat = test - np.mean(np.vstack((other1, other2)),axis=0)
    return testStat


def testStat_OvO(test, foil, blank):
    '''Compute one vs one test stat. Subtract foil from test.
    So a neg value means mod of foil, potentially. Naming
    test stat after test txt not really true therefore
    Blank param is for third txt. Not used here but for other test stats 
    so keep structure'''
    return test - foil

def computeTestStatistic(testFx, traceTest, traceOther1, traceOther2):
    '''given summary trace (e.g. avg) for each txt (grouped in some manner)
    compute a test statistic trace on each and return them in A,B,H order'''
    testStat = testFx(traceTest, traceOther1, traceOther2)
    return testStat

def gotoKPyDir():
    user = os.getlogin()
    os.chdir(f"c:/Users/{user}/Google Drive/Python_Code/KLab/mts_analysis/")


def getPosFromTs(yourts,position):
    '''position is the full list of ts, yourts are the times you'd like a pos for '''
    adjTs = [takeClosest(position[:,0],i) for i in yourts]
    target = np.empty((0,2))
    for i in adjTs:
        target = np.vstack((target,position[np.where(position[:,0]==i)][:,1:][0]))
    return target

def getTsinInverval(array,start,stop):
    '''left inclusive'''
    idx = np.where(np.logical_and(array >= start, array < stop))
    return array[idx[0]]


def getStimFileName(datafile):
    for f in os.listdir(datafile):
        if f.endswith('_Stimuli.p'):
            return f

def adjustClust(clust, ts, pos):
    ''' adjust for the 640 dim of camera. remove data
    <50 in either x,y as this is spurious. R765 stil dont
    know origin of off-camera streaking.
    and return as a np array'''
    ts = np.asarray(sorted(list(pos.keys())))
    ts = np.asarray(sorted(list(pos.keys())))
    posx = [640 - pos[i][0] for i in ts]
    posy = [pos[i][1] for i in ts]
    for t in ts:
        pos[t][0] = 640 - pos[t][0]
    spikexy = getPosFromTs(clust)
    spikes = np.column_stack((clust,spikexy))
    spikes = spikes[np.where(spikes[:,2] > 50)]
    return spikes

def getVisitPos(alley,visitnum,array,alleyVisits, position):
    '''array is either spike ts or pos ts'''
    tsVisit = getTsinInverval(array, alleyVisits[alley][visitnum][0], alleyVisits[alley][visitnum][1])
    posVisit = getPosFromTs(tsVisit, position)
    return np.column_stack((tsVisit,posVisit))


def list_to_arrays(l):
    arr = np.empty((0,l[0].shape[1]))
    for subarray in l:
        arr = np.vstack((arr,subarray))
    return arr

#def getVisits(alley,txt, alleyByTxt):
#    spikes, pos, idx = alleyByTxt[alley][txt]['spk'],alleyByTxt[alley][txt]['pos'],alleyByTxt[alley][txt]['idx']
#    sv, pv = np.split(spikes,idx[:,0]), np.split(pos,idx[:,1]) # list of subarrays
#    return sv[:-1], pv[:-1]


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

def get_long_alley_dim():
    '''all alleys should have same dim, as theyre same size. so use alley  1 by convention'''
    r,c = alleyBins[0]['rows'].shape[0], alleyBins[0]['cols'].shape[0]
    return max(r,c)

def getAxType(array):
    '''should "axis" = 1 or 2? Which dim is longer'''
    if array.shape[0] > array.shape[1]:
        return 1
    elif array.shape[0] < array.shape[1]:
        return 0

def weird_smooth(U,sigma):
    V=U.copy()
    V[U!=U]=0
    VV=sp.ndimage.gaussian_filter(V,sigma=sigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=sp.ndimage.gaussian_filter(W,sigma=sigma)

    Z=VV/WW
    return Z

def getClustList(datafile):
    """
    Given an absolute path to a data directory
    return a list of clusters in format
    TT{1-16}\\cl-maze1.n
    """
    clustList = []
    for subdir, dirs, fs in os.walk(datafile):
        for f in fs:
            if 'cl-maze' in f and 'OLD' not in f and 'Undefined' not in f:
                clustname = subdir[subdir.index("TT"):] + "\\" + f
                clustList.append(clustname)
    return clustList

def readUnitsToLoad():
    """
    In Google drive > KLab > Ratterdam > Code
    there is ratterdam_UnitsToLoad.txt
    This lists units to be used for whatever present analysis
    Minimal Format is (each line, a unit): rat,day,clust
    E.g: R765,RFD7,TT14\cl-maze1.1 (note how clust is defined, \\ added when readin)
    These are loaded into a list of dicts where each dict
    has keys: 'rat', 'day', 'cluster'
    
    If any entry has more than rat,day,clust that info will be read into
    additional fields of the dict. Note that any additional data in the text
    file should be of key:value format so the readin has that info.
    
    To comment out a unit without deleting text, add leading # to line
    """
    if socket.gethostname() == 'Tolman':
        codeDirBase = 'C:\\Users\\whockei1\\Google Drive\\KnierimLab\\Ratterdam\\Code\\'
    
    elif socket.gethostname() == 'DESKTOP-BECTOJ9':
        codeDirBase = 'C:\\Users\\whock\\Google Drive\\KnierimLab\\Ratterdam\\Code\\'
    
    file = open(codeDirBase+"ratterdam_UnitsToLoad.txt","r")
    unitlist = file.readlines()
    clusts = []
    if len(unitlist) > 0:
        for unit in unitlist:
            if unit[0] != '#': #check if commented out
                unit = unit.split(",")
                rat = unit[0]
                day = unit[1]
                clust = unit[2].rstrip() #removes all right trailing whitespace. neat.
                d = {'rat':rat, 'day':day, 'cluster':clust}
    
                
                #check for additional info, maybe analysis wants a specific alley or something
                if len(unit) > 3: #already split, first 3 are stored above
                    for data in unit[3:]:
                        key, value = data.split(":")
                        d[key] = value.rstrip()
                
                clusts.append(d)
            
    
    return clusts
        
        
def makeCustomColormap(nb=100,name='mymap',c=[]):

    if c ==[]:
        c = [(0,0,0.5),(1,1,0),(1,0,0)]
    mycm = LinearSegmentedColormap.from_list(name,c,N=nb)
    return mycm

def calcSmartMax(array2d, cutoff=0.98, scale=2.5,bins=100):
    """
    Given array where each row is a sample, eg a lin rate map
    find a good max visualiz. value for eg. imshow across all samples
    by getting percentile cutoff and boosting it by scale factor
    
    Bins is tricky, depends on how many rate maps
    go into the analysis. 100 bins is good. 0.98 cutoff.
    """
    frs = []
    for row in array2d:
        frs.extend(row)
    frs = np.asarray(frs)
    frs = frs[np.isfinite(frs)]
    h,b = np.histogram(frs, bins=bins)
    frcum = np.cumsum(h)
    propExp = np.asarray([i/h.sum() for i in frcum])
    try:
        thresh = np.where(propExp < cutoff)[0][-1]
    except:
        thresh = np.where(b == np.median(b))
    return b[thresh]*scale

def checkCloserPoint(p1, p2, pt):
    """
    Given two points p1, p2
    where each is [x,y]
    see which pt is closer to
    (also of form [x,y])
    
    Return string "first" or "second"
    meaning its closer to first point arg
    or second point arg. If equal return "error"
    """
    d1 = npNorm(p1 - pt)
    d2 = npNorm(p2 - pt)
    if d1 < d2:
        return "first"
    elif d2 < d1:
        return "second"
    else:
        return None

def checkVisitSide(visitOccs, bounds, vType):
    """
    Check which side rat entered/exited alley
    visitOccs is [ts,x,y] arr for 1 visit
    bounds is [ul, ll, ur, lr] for alley in question
    vType is whether you want the "entry" side or "exit" side
    Return a label "SW" or "NE"
    """
    if vType == 'entry':
        point = visitOccs[0,1:]
    elif vType == 'exit':
        point = visitOccs[-1,1:]
    ll, ur = bounds[1], bounds[2]
    closerPt = checkCloserPoint(ll, ur, point)
    if closerPt is not None:
        if closerPt == "first":
            side = "SW"
        elif closerPt == "second":
            side = "NE"
    else:
        side = None
    return side

def extractCorners(givenAlleyBounds):
    """Alley bounds gives [[x1, x2], [y1, y2]]. Convert that
    to UL, LL, UR, LL (xn, ym) points n,m <- [1,2]
    ul - x1, y2
    ll - x1, y1
    ur = x2, y2
    lr = x2, y1
    
    Returns ul, ll, ur, lr
    """
    b = givenAlleyBounds # for ease of typing
    ul, ll, ur, lr = [b[0][0], b[1][1]], [b[0][0], b[1][0]], [b[0][1], b[1][1]], [b[0][1], b[1][0]]
    return ul, ll, ur, lr

def midpoint(p1, p2):
    return (p1[0]+p2[0])/2, (p1[1]+p2[1])/2

def stepsmooth(vec, smoothval=0.5):
    nvec = vec[np.where(~np.isnan(vec))]
    svec = np.zeros(nvec.shape[0])
    fvec = np.zeros(vec.shape[0])
    svec[-1] = (nvec[-1]+smoothval*nvec[-2])/2
    svec[0] = (nvec[0]+smoothval*nvec[1])/2
    for i in range(1,nvec.shape[0]-1):
        svec[i] = ((smoothval*nvec[i-1])+nvec[i]+(smoothval*nvec[i+1]))/3
    fvec[np.where(np.isnan(vec))] = np.nan
    fvec[np.where(~np.isnan(vec))] = svec
    return fvec


def makeRM(spikes, position, bins=[30,50]):
    """ this is 2d whole track. uses g"""
    r,c = bins
    hs = np.histogram2d(spikes[:,2],spikes[:,1],bins=[r,c])[0]
    ho = np.histogram2d(position[:,2],position[:,1],bins=[r,c])[0]
    n = (hs*np.reciprocal(ho))*30
    n[np.where(ho==0)] = np.nan
    n = weird_smooth(n, Def.smoothing_2d_sigma)
    n[np.where(ho==0)] = np.nan
    return n

def checkInAlley(array, alley):
    """
    Takes a np.array(ts,x,y) of spikes or occs
    And filters to exlude those outside the alleybounds of alley (pass 1-indexed alley num,
    in full-ratterdam numbering i.e. out of 17).
    Returns np.array of same num cols, possibly fewer rows
    """
    ab = Def.alleyBounds[alley-1]
    xmin, xmax, ymin, ymax = ab[0][0], ab[0][1], ab[1][0], ab[1][1]
    new_array = array[(array[:,1] > xmin)&
              (array[:,1] <= xmax)&
              (array[:,2] > ymin)&
              (array[:,2] <= ymax)]
    return new_array
        
        