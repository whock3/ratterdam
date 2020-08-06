# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:15:18 2019

@author: whockei1

GUI For Ratterdam Beltway
"""

# imports
from matplotlib import pyplot as plt
import os, numpy as np, time
import matplotlib.gridspec as gridspec

def readInData():
    data = open(basepath+fname)
    line = data.read()
    data.close()
    return line

def turnOffStuff():
    for metric in ['strInfo', 'lapTimes', 'stimInfo']:
        ax = beltway_structure['metrics'][metric]['ax']
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        
        
def show_plot(figure_id=None): 
    """found on SO. https://stackoverflow.com/questions/8202228/make-matplotlib-plotting-window-pop-up-as-the-active-one"""
    if figure_id is None:
        fig = plt.gcf()
    else:
        # do this even if figure_id == 0
        fig = plt.figure(num=figure_id)

    plt.show()
    plt.pause(1e-9)
    fig.canvas.manager.window.activateWindow()
    fig.canvas.manager.window.raise_()


def getTxtRatios():
    a, b, c = rtxt_hists['A'], rtxt_hists['B'], rtxt_hists['C']
    total = a+b+c
    return a/total, b/total, c/total


basepath = 'C:\\Users\\admin\\Desktop\\ratterdam_code\\'
fname = 'beltway_gui_data.txt'
modtime = os.path.getmtime(basepath+fname)
colorlookup = {'B':'lightgrey', 'C':'darkolivegreen', 'A':'saddlebrown'}


plt.rc('axes', linewidth=4)
aspectRatio = 10/7 #hardcoded
scl = 10 # smaller size
gsframe = gridspec.GridSpec(7,10) # if intersections are 1x1 units and 1x2 units. Not totally to scale
regions = [1,2,3,4,5,6,7,8,9,'startbox']
beltway_structure = {i:{} for i in regions}
fig = plt.figure(figsize=(int(scl*aspectRatio),scl))
# tedious but straightforward
beltway_structure[1]['ax'] = plt.subplot(gsframe[6,1:3])
beltway_structure[2]['ax'] = plt.subplot(gsframe[4:6, 0])
beltway_structure[3]['ax'] = plt.subplot(gsframe[1:3, 0])
beltway_structure[4]['ax'] = plt.subplot(gsframe[0, 1:3])
beltway_structure[5]['ax'] = plt.subplot(gsframe[0, 4:6])
beltway_structure[6]['ax'] = plt.subplot(gsframe[0, 7:9])
beltway_structure[7]['ax'] = plt.subplot(gsframe[1:3, 9])
beltway_structure[8]['ax'] = plt.subplot(gsframe[4:6, 9])
beltway_structure[9]['ax'] = plt.subplot(gsframe[6, 7:9])
beltway_structure['startbox']['ax'] = plt.subplot(gsframe[6,4:6])

for i in regions:
    beltway_structure[i]['ax'].set_xticklabels([])
    beltway_structure[i]['ax'].set_yticklabels([])

beltway_structure['metrics'] = {}
for metric in ['strInfo', 'lapTimes', 'stimInfo']:
    beltway_structure['metrics'][metric] = {}
beltway_structure['metrics']['strInfo']['ax']= plt.subplot(gsframe[2,2])
beltway_structure['metrics']['lapTimes']['ax'] = plt.subplot(gsframe[2,4:9])
beltway_structure['metrics']['stimInfo']['ax'] = plt.subplot(gsframe[4,2])
turnOffStuff()
plt.suptitle("RATTERDAM GUI", fontsize=22)


txt_hist = []
r_hist = []
lapTimes = []
rtxt_hists = {'A':0, 'B':0, 'C':0}
lapNum = 1
exit = False
initTime = time.time()
lastModTime = initTime
while not exit:
    modTime = os.path.getmtime(basepath+fname)
    if modTime > lastModTime:
        lapTimes.append(modTime - lastModTime)
        data = readInData()
        txts = [i[0] for i in data.split(",")[:-1]]
        rs = [i[1] for i in data.split(",")[:-1]]
        txt_hist.append(txts)
        r_hist.append(rs)

        for i, txt in enumerate(txts):
            beltway_structure[i+1]['ax'].set_facecolor(colorlookup[txt])
            if rs[i] == '1':
                rtxt_hists[txt] += 1
                beltway_structure[i+1]['ax'].spines['bottom'].set_color('red')
                beltway_structure[i+1]['ax'].spines['top'].set_color('red')
                beltway_structure[i+1]['ax'].spines['right'].set_color('red')
                beltway_structure[i+1]['ax'].spines['left'].set_color('red')

        beltway_structure['metrics']['lapTimes']['ax'].plot(lapTimes,marker='o',linewidth=2,color='k')

        text = f"Lap {lapNum}"
        beltway_structure['metrics']['strInfo']['ax'].clear()
        turnOffStuff()
        beltway_structure['metrics']['strInfo']['ax'].text(0,0.5,text,size=16)

        a,b,c = getTxtRatios()
        beltway_structure['metrics']['stimInfo']['ax'].bar([1,2,3],[a,b,c],color=['saddlebrown','lightgrey','darkolivegreen'])

        #fig.canvas.draw()
        show_plot()
        plt.pause(0.1)

        lastModTime = modTime
        lapNum += 1

        if data[-1] == 'X':
            exit = True
