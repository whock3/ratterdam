# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 23:33:13 2021

@author: whockei1

Scrap file
"""
#%% Functions to load
def checkCheating(unit, turn, turn_np1):
    """
    deprecated. solving problem a different way
    
    Check to see if any of the unit's fields overlap with either of the non-local
    regions involved in this turn (i.e. the pre and post alleys were decoding dir
     from) and only if not do we let this unit contribute data to this turn
    """
    threshold = 0.10
    anyExcessiveOverlap = False
    for perim in unit.perimeters:
        for a in [turn['Alley-'], turn_np1['Alley+']]:
            alley = ratborders.alleyInterBounds[str(a)]
            xmin,xmax = alley[0][0],alley[0][1]
            ymin,ymax = alley[1][0],alley[1][1]
            alleyPerim = np.array([[xmin, ymin], [xmax, ymin],
                                           [xmax, ymax], [xmin, ymax]])
            alleySize = (xmax-xmin) * (ymax-ymin)
            overlap1 = repPC.overlap(perim, alleyPerim)
            if overlap1 > alleySize*threshold:
                anyExcessiveOverlap = True
            
    return anyExcessiveOverlap



#%% Turn based rate maps - setup
# Group data by allo dir pre (4 plots) and allo dir post (4 plots) +/- 1.5 turns
#(the half turn is bc the turns end at the intersection so punch out a bit more
# in time to get the full +/1 1 turn)
#Reminder about code order: N,E,S,W,  F,R,B,L


# turnRMS = {'Pre':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
#                   }, 
           
#            'Post':{'1':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '2':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '3':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},
#                   '4':{'Spikes':np.empty((0,3)), 'Pos':np.empty((0,3))},  
#                   }
#            }

# for field in unit.fields:
#     for visit in field:
#         turnIdx = np.argmin(np.abs(turns['Ts exit'].astype(np.double)-visit[0]))
#         try:
#             # first/last turn has no previous/next turn, so just ignore it
#             if turnIdx > 2 and turnIdx < turns.shape[0]-2:
#                 spikesPre = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
#                 spikesPost = unit.spikes[(unit.spikes[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.spikes[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

#                 turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Spikes'],spikesPre))
#                 turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Spikes'],spikesPost))
    
#                 occPre = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx-2]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx]['Ts exit']))]
#                 occPost = unit.position[(unit.position[:,0]>float(turns.iloc[turnIdx]['Ts exit']))&(unit.position[:,0]<=float(turns.iloc[turnIdx+2]['Ts exit']))]

#                 turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'] = np.vstack((turnRMS['Pre'][turns.iloc[turnIdx]['Allo-']]['Pos'],occPre))
#                 turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'] = np.vstack((turnRMS['Post'][turns.iloc[turnIdx]['Allo+']]['Pos'],occPost))
#         except:
#             print("Likely invalid code found")

# for d in ['Pre', 'Post']:
#     fig, ax = plt.subplots(2,2)
#     for i,code in enumerate([('1','N'),('2','E'),('3','S'),('4','W')]):
#         s,o = turnRMS[d][code[0]]['Spikes'], turnRMS[d][code[0]]['Pos']
#         rm = util.makeRM(s,o)
#         fig.axes[i].imshow(rm,aspect='auto',interpolation='None',cmap=cmap,vmax=7,origin='lower')
#         fig.axes[i].set_title(code[1])
#     plt.suptitle(f"{d}-turn Bearing")
    