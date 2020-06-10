# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:30:28 2019

@author: whockei1
"""

for ev in bounds:
    plt.vlines(ev[0],0,700)
    plt.annotate(ev[1],xy=(ev[0],800),xytext=(ev[0],800),rotation=90)
plt.ylim([0,800])
plt.scatter(position[:,0],position[:,1],c='g')
plt.scatter(position[:,0],position[:,2],c='r')

for i in range(len(bounds)-1):
    seg = position[np.where(np.logical_and(position[:,0] > bounds[i][0], position[:,0] < bounds[i+1][0]))]
    x,y = seg[:,1], seg[:,2]
    x = [i for i in x if i > 10]
    y = [i for i in y if i < 470]
    print(bounds[i][1])
    print(f"mean x: {round(np.mean(x),2)}")
    print(f"mean y: {round(np.mean(y),2)}")
    print("----------------------------")


    
ab = {}
for i in range(len(bounds)-1):
    seg = position[np.where(np.logical_and(position[:,0] > bounds[i][0], position[:,0] < bounds[i+1][0]))]
    x,y = seg[:,1], seg[:,2]
    x = np.mean([i for i in x if i > 10])
    y = np.mean([i for i in y if i < 470 and i > 50])
    x = float('%.2f'%(x))
    y = float('%.2f'%(y))
    meta = str(bounds[i][1]).split("_")
    if "alley" in meta[0]:
        alley = int(meta[1])-1
        if alley not in ab:
            ab[alley] = [[0,0], [0,0]]
        case = meta[2]
        if "UL" in case:
            ab[alley][0][0] = x
            ab[alley][1][1] = y
        elif "LL" in case:
            ab[alley][0][0] = x
            ab[alley][1][1] = y
        elif "UR" in case:
            ab[alley][0][1] = x
            ab[alley][1][1] = y
        elif "LR" in case:
            ab[alley][0][1] = x
            ab[alley][1][0] = y