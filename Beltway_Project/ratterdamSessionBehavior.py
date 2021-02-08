import numpy as np, matplotlib.pyplot as plt, datetime, random, os
import json, collections, itertools

data = json.load(open("R765_180503_15-47.json","r"))
rewards = [i for i in data['events'] if i == 'reward']
times = [x for i,x in enumerate(data['ts']) if data['events'][i] == 'reward']
sides = [i for i in data['sides']]

AlleyStateGraph = {0:{'A':[3, 4], 'B':[2]},
                                1:{'A':[3, 12, 14], 'B':[2, 16]},
                                2:{'A':[1, 16], 'B':[0]},
                                3:{'A':[1, 12, 14], 'B':[0, 4]},
                                4:{'A':[5, 6], 'B':[0, 3]},
                                5:{'A':[8, 11, 12], 'B':[4, 6]},
                                6:{'A':[7], 'B':[4, 5]},
                                7:{'A':[8, 9], 'B':[6]},
                                8:{'A':[7, 9], 'B':[5, 11, 12]},
                                9:{'A':[10], 'B':[7, 8]},
                                10:{'A':[9], 'B':[11, 13]},
                                11:{'A':[13, 10], 'B':[5, 8, 12]},
                                12:{'A':[5, 8, 11], 'B':[3, 14, 1]},
                                13:{'A':[10, 11], 'B':[14, 15]},
                                14:{'A':[13, 15], 'B':[1, 3, 12]},
                                15:{'A':[13, 14], 'B':[16]},
                                16:{'A':[15], 'B':[1, 2]}
        }

rawvisits = []
rawvisits.append((0, data['alleys'][0]))
for idx in range(1,len(data['alleys'])):
    if data['alleys'][idx] != data['alleys'][idx-1]:
        rawvisits.append((idx, data['alleys'][idx]))

passes = []
for visit in rawvisits:
    if isPassorNot(visit):
        passes.append(visit)
    

def isPassorNot(visit):
    match  = False
    initBeam = sides[visit[0]]
    for i,r in enumerate(data['alleys'][visit[0]:]):
        if r == other(initBeam):
            potMatch = i
            break
    if not any([i in data['alleys'][visit[0]:potMatch] for i in visit[1]+AlleyStateGraph[visit[1]][other(initBeam)]]):
        travseral = True
    return traversal
            
