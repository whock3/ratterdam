from sklearn import decomposition
from sklearn.neighbors import NearestCentroid
from sklearn import preprocessing
from scipy.spatial.distance import pdist
import numpy as np
import ratterdam_Defaults as Def


def interCentroidDistance(ncObj):
    """
    Input a NearestCentroid object (sklearn)
    and return pairwise Euclidian distances 
    between them
    """
    return pdist(ncObj.centroids_)


def avgDistToCentroid(ncObj,sa,sb,sc):
    """
    Input: ncObj - NearestCentroid object (sklearn)
           sa,sb,sc - vectors which produced centroids
    """
    avga = np.mean(np.linalg.norm(ncObj.centroids_[0]-sa,axis=1))
    avgb = np.mean(np.linalg.norm(ncObj.centroids_[1]-sb,axis=1))
    avgc = np.mean(np.linalg.norm(ncObj.centroids_[2]-sc,axis=1))
    return avga, avgb, avgc


def interSampleDistance(sa,sb,sc):
    """
    Input sa,sb,sc - samples from a,b,c trials
    Return avg Euclidian distance within each set of samples
    """
    avga = np.mean(pdist(sa))
    avgb = np.mean(pdist(sb))
    avgc = np.mean(pdist(sc))
    return avga, avgb, avgc
	
def setupData(unit, alley, features, dimred):
    """
    Input   - unit: Unit class object
            - alley: alley in whole track terms
            - list of features (e.g. COM, Max, etc) to represent a single RM
                or pass 'rm' to use original rm 
            - dimred: if False, no dimensionality reduction
                if a number then thats # pca comps

    Returns: data matrix X (n,f) n samples f features
             label vector Y (n,)
    """
    atrials,btrials,ctrials = np.empty((0,len(features))), np.empty((0,len(features))), np.empty((0,len(features)))
    for i,visit in enumerate(unit.alleys[alley]):
        stim = visit['metadata']['stimulus']
        feats = createVisitSummaryFeatures(unit,alley,i,features)
        if stim == 'A':
            atrials = np.vstack((atrials, feats))
        elif stim == 'B':
            btrials = np.vstack((btrials, feats))
        elif stim == 'C':
            ctrials = np.vstack((ctrials, feats))
            
    X = np.vstack((atrials, btrials, ctrials))
    X[np.where(~np.isfinite(X))] = 0
    X = preprocessing.StandardScaler().fit_transform(X)
    Y = np.asarray([0]*atrials.shape[0] + [1]*btrials.shape[0] + [2]*ctrials.shape[0])
    if dimred:
        if dimred > len(features):
            print("Error - PCA ncomp > original vector size")
        else:
            pca = decomposition.PCA(n_components=dimred)
            pca.fit(X)
            X = pca.transform(X)
    return X, Y


def splitSamplesbyLabel(X,y):
    """
    Given input matrix X with samples of different
    labels, stored in Y, split them into arrays by label
    A=0, B=1, C=2 for texture labels by convention
    """
    a = X[y==0]
    b = X[y==1]
    c = X[y==2]
    return a,b,c


def findCentroids(X,Y):
    """
    Given X (n,f) and Y (n,)
    create NearestCentroid object 
    and return it
    """
    nc = NearestCentroid(metric='euclidean')
    nc.fit(X,Y)
    return nc
	
	
	
"""	
# Intersample distance	
npca=False
X,Y = temp.setupData(unit, alley, npca)
cent = temp.findCentroids(X,Y)
at,bt,ct = temp.splitSamplesbyLabel(X,Y)
ad,bd,cd  = temp.interSampleDistance(at,bt,ct)
mind = min(ad,bd,cd)

nshuff=1000
ss = np.empty((0))
for s in range(nshuff):
    Ys = np.random.permutation(Y)
    scent = temp.findCentroids(X,Ys)
    sat,sbt,sct = temp.splitSamplesbyLabel(X,Ys)
    ssd = np.min(temp.interSampleDistance(sat,sbt,sct))
    ss = np.append(ss,ssd)


plt.hist(ss)
plt.vlines(mind,0,200,'r')
plt.vlines(np.percentile(ss,5),0,200,'k')
plt.title(f"{unit.name} MDC pca comp = {npca}")


# Distance to centroids
npca=9
X,Y = temp.setupData(unit, alley, npca)
cent = temp.findCentroids(X,Y)
at,bt,ct = temp.splitSamplesbyLabel(X,Y)
ad,bd,cd  = temp.avgDistToCentroid(cent,at,bt,ct)
mind = min(ad,bd,cd)

nshuff=1000
ss = np.empty((0))
for s in range(nshuff):
    Ys = np.random.permutation(Y)
    scent = temp.findCentroids(X,Ys)
    sat,sbt,sct = temp.splitSamplesbyLabel(X,Ys)
    ssdc = np.min(temp.avgDistToCentroid(scent,sat,sbt,sct))
    ss = np.append(ss,ssdc)


plt.hist(ss)
plt.vlines(mind,0,200,'r')
plt.vlines(np.percentile(ss,5),0,200,'k')
plt.title(f"{unit.name} MDC pca comp = {npca}")

# intercentroid distance
npca=3
X,Y = temp.setupData(unit, alley, npca)
cent = temp.findCentroids(X,Y)
icd  = temp.interCentroidDistance(cent)

nshuff=1000
ss = np.empty((0))
for s in range(nshuff):
    Ys = np.random.permutation(Y)
    scent = temp.findCentroids(X,Ys)
    sicd = np.max(temp.interCentroidDistance(scent))
    ss = np.append(ss,sicd)


plt.hist(ss)
plt.vlines(max(icd),0,200,'r')
plt.vlines(np.percentile(ss,95),0,200,'k')
plt.title(f"{unit.name} ICD pca comp = {npca}")

"""