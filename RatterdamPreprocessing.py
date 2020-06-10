import numpy as np, matplotlib.pyplot as plt, functools, pickle, json, datetime
# alley :[[xmin, xmax], [ymin, ymax]]
alleyBounds = {1:[[395, 495],[352, 390]],
               2:[[395, 495],[218, 255]],
               3:[[495, 535],[255, 352]],
               4:[[360, 395],[255, 352]],
               5:[[260, 360],[352, 390]],
               6:[[216, 260],[255, 352]],
               7:[[120, 216],[352, 390]],
               8:[[87, 120],[255, 352]],
               9:[[120, 216],[218, 255]],
               10:[[87, 120],[115, 218]],
               11:[[120, 216],[75, 115]],
               12:[[216, 260],[115, 218]],
               13:[[260, 360],[218, 255]],
               14:[[260, 360],[75, 115]],
               15:[[360, 395],[115, 218]],
               16:[[395, 495],[75, 115]],
               17:[[495, 535],[115, 218]]}

def getNLSwapTimes(rec):
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


def convertPytoNLts(rec, nltimes, pydata):
    '''nl times cleaned first''' 
    pySwapTimes = [i[0] for i in pydata['order']]
    adjPy = [i - pySwapTimes[0]for i in pySwapTimes]
    projs = [min(enumerate(nltimes), key=lambda x: abs(x[1]-i)) for i in [j+nltimes[0] for j in adjPy]]
    swapTS = [nltimes[i[0]] for i in projs]
    beginTS = [i[3] for i in rec if b'begin session' in  i[10]][0]
    endTS = [i[3] for i in rec if b'end session' in i[10]][0]
    alleySwapTS = {i:[] for i in range(1,18)}
    for i in range(len(swapTS)):
        alleySwapTS[pydata['order'][i][1]].append(swapTS[i])
    for i in range(1,18):
        alleySwapTS[i] = [beginTS] + alleySwapTS[i] + [endTS]
    return alleySwapTS
    
def getSingleAlleyTxtResponse(alley):
    spikesbyTxt = {'A':np.empty((0,2)), 'B':np.empty((0,2)), 'C':np.empty((0,2))}
    for i in range(len(alleySwapTS[alley])-1):
        txt = pydata['stimuli'][alley][i][0]
        z = functools.reduce(np.logical_and,(spikes > alleySwapTS[alley][i], spikes < alleySwapTS[alley][i+1], spikex > alleyBounds[alley][0][0], spikex < alleyBounds[alley][0][1], spikey > alleyBounds[alley][1][0], spikey < alleyBounds[alley][1][1]))
        q = np.stack((spikex[z],spikey[z]),axis=-1)
        spikesbyTxt[txt] = np.vstack((spikesbyTxt[txt],q))
    return spikesbyTxt    

pos = util.read_pos(datafile)
ts = list(pos.keys())
posx = [pos[i][0] for i in ts]
posy = [480 - pos[i][1] for i in ts]
spikes = np.asarray(spikes)
spikex = np.asarray([pos[takeClosest(ts,i)][0] for i in spikes])
spikey = np.asarray([480 - pos[takeClosest(ts,i)][1] for i in spikes])


def getSingleAlleyTxtPos(alley):
    posbyTxt = {'A':np.empty((0,2)), 'B':np.empty((0,2)), 'C':np.empty((0,2))}
    for i in range(len(alleySwapTS[alley])-1):
        txt = pydata['stimuli'][alley][i][0]
        z = functools.reduce(np.logical_and,(ts > alleySwapTS[alley][i], ts < alleySwapTS[alley][i+1], posx > alleyBounds[alley][0][0], posx < alleyBounds[alley][0][1], posy > alleyBounds[alley][1][0], posy < alleyBounds[alley][1][1]))
        q = np.stack((posx[z],posy[z]),axis=-1)
        posbyTxt[txt] = np.vstack((posbyTxt[txt],q))
    return posbyTxt

def alleyRM(i,txt,orien):
    """occ normed """
    if orien == 'v':
        bins = [50,75]
        ax = 0
    elif orien == 'h':
        bins = [75, 50]
        ax = 1
    spk = np.histogram2d(allalleys[i][txt][:,0], allalleys[i][txt][:,1], bins = bins)[0]
    pos = np.histogram2d(allalleysPos[i][txt][:,0], allalleysPos[i][txt][:,1], bins = bins)[0]
    rm = spk/gauss(pos, 1)
    return rm

allalleys = {i: getSingleAlleyTxtResponse(i) for i in range(1,18)}
allalleysPos = {i: getSingleAlleyTxtPos(i) for i in range(1,18)}
alleysByTxt = {'A':np.empty((0,2)), 'B':np.empty((0,2)), 'C': np.empty((0,2))}
alleysByTxtPos = {'A':np.empty((0,2)), 'B':np.empty((0,2)), 'C': np.empty((0,2))}
for txt in ['A', 'B', 'C']:
    for i in range(1,18):
        alleysByTxt[txt] = np.vstack((alleysByTxt[txt], allalleys[i][txt]))
        alleysbyTxtPos[txt] = np.vstack((alleysByTxtPos[txt], allalleysPos[i][txt]))


for i, t in enumerate(['A','B','C']):
    plt.figure(i+1)
    spkhist = np.nan_to_num(np.histogram2d(alleysByTxt[t][:,0], alleysByTxt[t][:,1],bins=[np.linspace(80,550,150), np.linspace(70,400,150)])[0])
    poshist = np.nan_to_num(np.histogram2d(alleysByTxtPos[t][:,0], alleysByTxtPos[t][:,1], bins=[np.linspace(80,550,150), np.linspace(70,400,150)])[0])
    plt.imshow(spkhist/gauss(poshist,2))
    plt.colorbar()
    plt.title(t,fontsize=30)

a = 1
alley = {'A':alleyRM(a,'A','h'), 'B':alleyRM(a,'B','h'), 'C':alleyRM(a,'C','h')}
for i,t in enumerate(['A','B','C']):
    plt.figure(i+1)
    plt.imshow(alley[t])
    plt.colorbar()
    plt.title(t)
