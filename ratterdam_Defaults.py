"""
This file contains default parameters for main ratterdam scripts.
Values can be changed to alter analysis pipeline behavior
"""
import sys
import newAlleyBounds as nab

sys.path.insert(0, 'E:\\UserData\\Documents\\GitHub\\ratterdam\\')

alleyBounds_byEye = {0:[[150, 250],[350, 400]],
               1:[[150, 250],[215, 250]],
               2:[[115, 150],[250, 350]],
               3:[[250, 285],[250, 350]],
               4:[[285, 385],[350, 400]],
               5:[[385, 415],[250, 350]],
               6:[[415, 515],[350, 400]],
               7:[[515, 560],[250, 350]],
               8:[[415, 515],[215, 250]],
               9:[[515, 560],[115, 215]],
               10:[[415, 515],[75, 115]],
               11:[[385, 415],[115, 215]],
               12:[[285, 385],[215, 250]],
               13:[[285, 385],[75, 115]],
               14:[[250, 285],[115, 215]],
               15:[[150, 250],[75, 115]],
               16:[[115, 150],[115, 215]]}

#format is [x1, x2], [y1, y2]
alleyBounds_DEP = {0:[[165, 210],[345, 385]],
                   1:[[150, 215],[210, 250]],
                   2:[[105, 140],[272, 340]],
                   3:[[230, 270],[275, 340]],
                   4:[[295, 359],[352, 390]],
                   5:[[370, 410],[260, 327]],
                   6:[[430, 492],[347, 390]],
                   7:[[500, 540],[260, 325]],
                   8:[[420, 480],[210, 254]],
                   9:[[500, 550],[140, 205]],
                   10:[[415, 480],[75, 115]],
                   11:[[370, 410],[135, 205]],
                   12:[[280, 350],[210, 250]],
                   13:[[280, 350],[75, 115]],
                   14:[[230, 275],[132, 205]],
                   15:[[165, 230],[75, 115]],
                   16:[[105, 145],[115, 185]]}

alleyBounds_fromVideo = {0:[[150, 240],[310, 355]],
               1:[[150, 240],[190, 230]],
               2:[[115, 160],[230, 315]],
               3:[[240, 285],[230, 315]],
               4:[[285, 375],[310, 355]],
               5:[[375, 425],[230, 315]],
               6:[[425, 510],[310, 355]],
               7:[[505, 555],[230, 315]],
               8:[[425, 510],[190, 230]],
               9:[[505, 555],[110, 190]],
               10:[[425, 510],[65, 110]],
               11:[[375, 425],[110, 190]],
               12:[[285, 375],[190, 230]],
               13:[[285, 375],[65, 110]],
               14:[[240, 285],[110, 190]],
               15:[[150, 240],[65, 110]],
               16:[[115, 160],[110, 190]]}

alleyBounds_fromPlot2 = {0:[[150, 240],[340, 400]],
             1:[[150, 240],[205, 265]],
             2:[[100, 150],[265, 340]],
             3:[[240, 290],[265, 340]],
             4:[[290, 375],[340, 400]],
             5:[[375, 425],[265, 340]],
             6:[[425, 505],[340, 400]],
             7:[[505, 555],[265, 340]],
             8:[[425, 505],[205, 265]],
             9:[[505, 555],[130, 205]],
             10:[[425, 505],[70, 130]],
             11:[[375, 425],[130, 205]],
             12:[[290, 375],[205, 265]],
             13:[[290, 375],[70, 130]],
             14:[[240, 290],[130, 205]],
             15:[[150, 240],[70, 130]],
             16:[[100, 150],[130, 205]]}

alleyBounds = alleyBounds_fromPlot2

alleyBounds_Macaulay = {0: [[164.36529411764707, 262.81470588235294], [367.5, 391.91]],
 1: [[182.77, 239.37], [215.3452941176471, 283.3547058823529]],
 2: [[129.96, 151.99], [255.75529411764705, 354.98470588235296]],
 3: [[231.25529411764705, 298.61470588235295], [278.05, 332.4]],
 4: [[291.43, 349.74], [341.25529411764705, 412.3447058823529]],
 5: [[360.86529411764707, 428.49470588235295], [277.6, 333.56]],
 6: [[438.1, 489.48], [343.2352941176471, 409.0247058823529]],
 7: [[481.99529411764706, 547.4247058823529], [295.3, 348.86]],
 8: [[438.97, 492.41], [215.79529411764707, 284.4147058823529]],
 9: [[481.4752941176471, 544.2747058823529], [170.81, 226.03]],
 10: [[437.75, 490.05], [92.78529411764706, 159.77470588235292]],
 11: [[363.01529411764704, 427.7647058823529], [146.28, 203.03]],
 12: [[313.93, 370.47], [213.54529411764707, 284.42470588235295]],
 13: [[313.43, 369.34], [86.63529411764706, 155.18470588235294]],
 14: [[235.48529411764707, 301.31470588235294], [143.64, 202.75]],
 15: [[186.3, 241.63], [85.55529411764705, 152.06470588235294]],
 16: [[111.01529411764706, 178.05470588235295], [146.5, 202.9]]}

alleyBounds_retroflective_onIRs = {0: [[178, 232], [357, 400]],
 1: [[177, 233], [225, 272]],
 2: [[123, 169], [283, 346]],
 3: [[245, 299], [285, 349]],
 4: [[306, 365], [357, 402]],
 5: [[376, 426], [288, 349]],
 6: [[436, 486], [360, 402]],
 7: [[501, 540], [290, 348]],
 8: [[439, 492], [235, 280]],
 9: [[503, 545], [166, 225]],
 10: [[437, 494], [109, 152]],
 11: [[381, 431], [158, 220]],
 12: [[310, 372], [225, 279]],
 13: [[310, 372], [99, 146]],
 14: [[250, 300], [149, 215]],
 15: [[183, 241], [95, 143]],
 16: [[125, 174], [151, 213]]}


alleyBounds_retroflective_onRim = {0: [[170, 235], [342, 415]],
 1: [[170, 235], [210, 292]],
 2: [[107, 174], [285, 350]],
 3: [[230, 310], [285, 350]],
 4: [[310, 365], [342, 415]],
 5: [[365, 445], [285, 350]],
 6: [[445, 503], [342, 415]],
 7: [[495, 556], [285, 350]],
 8: [[445, 503], [210, 292]],
 9: [[500, 560], [160, 215]],
 10: [[445, 503], [97, 158]],
 11: [[365, 445], [160, 215]],
 12: [[310, 365], [210, 292]],
 13: [[310, 365], [97, 158]],
 14: [[230, 310], [160, 215]],
 15: [[170, 235], [97, 158]],
 16: [[107, 174], [160, 215]]}



alleyBounds_retroflective_onRim_manuallyshifted1={0: [[160, 235], [332, 415]],
 1: [[160, 235], [200, 292]],
 2: [[97, 174], [275, 350]],
 3: [[220, 310], [275, 350]],
 4: [[300, 365], [332, 415]],
 5: [[355, 445], [275, 350]],
 6: [[435, 503], [332, 415]],
 7: [[485, 556], [275, 350]],
 8: [[435, 503], [200, 292]],
 9: [[490, 560], [150, 215]],
 10: [[435, 503], [87, 158]],
 11: [[355, 445], [150, 215]],
 12: [[300, 365], [200, 292]],
 13: [[300, 365], [87, 158]],
 14: [[220, 310], [150, 215]],
 15: [[160, 235], [87, 158]],
 16: [[97, 174], [150, 215]]}


alleyBounds_retroflective_onRim_manuallyshifted2={0: [[145, 250], [332, 415]],
 1: [[160, 235], [200, 292]],
 2: [[97, 174], [260, 365]],
 3: [[220, 310], [275, 350]],
 4: [[285, 380], [332, 415]],
 5: [[355, 445], [275, 350]],
 6: [[420, 518], [332, 415]],
 7: [[485, 556], [260, 365]],
 8: [[435, 503], [200, 292]],
 9: [[490, 560], [135, 230]],
 10: [[420, 518], [87, 158]],
 11: [[355, 445], [150, 215]],
 12: [[300, 365], [200, 292]],
 13: [[300, 365], [87, 158]],
 14: [[220, 310], [150, 215]],
 15: [[145, 250], [87, 158]],
 16: [[97, 174], [135, 230]]}

#format is [n1...n8] where each is loc of IR/alley point that defines intersection. ni = [xi, yi]
intersectionBounds = {'A':[[105, 385], [165, 385], [165, 385], [140, 340], [105, 340]],
                      'B':[[210, 385], [295, 390], [295, 352], [270, 340], [230, 340], [210, 345]],
                      'C':[[359, 390], [430, 390], [430, 347], [410, 327], [370, 327], [359, 352]],
                      'D':[[492, 390], [540, 390], [540, 325], [500, 325], [492, 347]],
                      'E':[[105, 272], [140, 272], [150, 250], [150, 210], [145, 185], [105, 185]],
                      'F':[[230, 275], [270, 275], [280, 250], [280, 210], [275, 205], [230, 205], [215, 210], [215, 250]],
                      'G':[[370, 260], [410, 260], [420, 254], [420, 210], [410, 205], [370, 205], [350, 210], [350, 250]],
                      'H':[[500, 260], [540, 260], [550, 205], [500, 205], [480, 210], [480, 254]],
                      'I':[[105, 115], [145, 115], [165, 115], [165, 75], [105, 75]],
                      'J':[[230, 132], [275, 132], [280, 115], [280, 75], [230, 75], [230, 115]],
                      'K':[[370, 135], [410, 135], [415, 115], [415, 75], [350, 75], [350, 115]],
                      'L':[[500, 140], [550, 140], [550, 75], [480, 75], [480, 115]]
                      }


#order of points below: left intersection midpoint, UL, LL, UR, LR, right intersection midpoint
alleyBounds_PolygonTriEnds = {0:[(),(),(),(),(),()],}


dataPathNub = "D:\\Ratterdam\\"
figurePathNub = "C:\\Users\\whock\\Google Drive\\KnierimLab\\Ratterdam\\Figures\\"

velocity_filter_thresh = 1.5#units are cm/s
ptsCm_krieger= 4.85 # number of camera points equiv to a cm

ptsCm_macaulay = 4.72
ptsCm = ptsCm_macaulay


#old bins = 15,30
cmPerBin = 1.5 # this is how many cm each bin should be long. This is chosen here by user. 
rewardZoneLength = 36 # size in cm of reward zone. This is hardcoded based on the alleybounds used
                      # and as of 2-12-21 using IAI (intersection-alley-intersection) def thats ~36cm long
longbin = round(rewardZoneLength/cmPerBin)
singleAlleyBins = [longbin+1,8]   # used [13,x] for single alley 
wholeAlleyBins = [round(int(480/ptsCm)),int(round(640/ptsCm))]

# Smoothing values to use
smoothing_2d_sigma = 2
smoothing_1d_sigma = 0.5# if gaussian smoothing this is sigma. If step smoothing where you average a bin by its nextdoor neighbors it is the scaling
                        #factor for those neigbors. i.e. (smoothing_1d_sigma*bin[i-1] + bin[i] + smoothing_1d_sigma*bin[i+1]) / 3

# Define visit parameters in milliseconds
# This pertains to visit alg used thru 12/18 wherein you group ts into chunks
# and say a chunk is a visit if it meets the following thresholds:
default_visit_gap= 30*1e6 #Gap between successive visits, i.e. exclude mult quick scans into alley
default_visit_duration = 0.1*1e6 #Minimum duration of visit


# imshow scale saturdation percentile cutoffs
wholetrack_imshow_pct_cutoff = 95
singlealley_imshow_pct_cutoff = 97

# toggle as to whether to include rewarded trials when constructing the Core.Unit() data structure
includeRewards = 2# flag to decide whether rewarded trials will be included in analysis. 0 - no rewards, 1 - only rewards, 2 - all trials


beltwayAlleys = [16,17,3,1,5,7,8,10,11]
beltwayAlleyLookup = {16:1, 17:2, 3:3, 1:4, 5:5, 7:6, 8:7, 10:8, 11:9}


fieldOverlapThresh_min = 0.3 # threshold for repetition analysis, deciding what 
                    # proportion of size of region the field overlaps with before
                    # its said to overlap that region. e.g. if thresh is 0.1
                    # then the field overlaps a region if it impinges on an area
                    # equal to 10% or more of the size of that region in surface area 
                    
                    
fieldOverlapThresh_max = 1.0   # 21-12-08, adding a max overlap argument for
                                # analyses where I just want to look at portions
                                #of fields encroaching slightly on an alley
                                # (eg 10-30% overlap). For existing analysis
                                #looking at fields overlapping significantly w 
                                # alley the max will simply be 1 


    
# Codes are, in order, North,East,South,West
# and Forward,Right,Back,Left (i.e. cw round circle centered on rat
#in both cases) for allo, ego respectively 

allocodedict = {'1':'N','2':'E','3':'S','4':'W','0':'X'}
egocodedict = {'1':'S','2':'R','3':'B','4':'L','0':'X'}


# Defaults for plotting, particularly for posters. Equates label sizes etc

ticksize = 52
ylabelsize = 64
xlabelsize = 64


# 2022-03-11 toggle for creating Unit() class objects for each neuron
# Choose whether all fields the watershed-based algorithm comes up with (toggle is True)
# or whether we filter based on my manual curation (fieldInclusionList file in each recording day directory
# and the manuallyRedrawnField folders) (toggle is False)

includeAllDetectedFields = False

