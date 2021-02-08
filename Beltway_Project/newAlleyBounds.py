# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:12:41 2020
@author: Ruo-Yah Lai
"""
import numpy as np
from collections import namedtuple


R781V = [118, 158, 249, 289, 382, 422, 502, 543]
R781H = [95, 138, 234, 275, 368, 409]
R808V = [100, 147, 236, 285, 375, 418, 504, 546]
R808H = [82, 135, 215, 269, 356, 404]
R859V = np.array([22, 32, 49, 61, 78, 88, 104, 113])*4.72
R859H = np.array([18, 28, 46, 56, 74, 84])*4.72


#R781
v1, v2, v3, v4, v5, v6, v7, v8 = R781V
h1, h2, h3, h4, h5, h6 = R781H

cityBlocks = np.array([[[v2,h2], [v3,h2], [v3,h3], [v2,h3], [v2,h2]],
                       [[v2,h4], [v3,h4], [v3,h5], [v2,h5], [v2,h4]],
                       [[v4,h2], [v5,h2], [v5,h3], [v4,h3], [v4,h2]],
                       [[v4,h4], [v5,h4], [v5,h5], [v4,h5], [v4,h4]],
                       [[v6,h2], [v7,h2], [v7,h3], [v6,h3], [v6,h2]],
                       [[v6,h4], [v7,h4], [v7,h5], [v6,h5], [v6,h4]]])

alleyLines = np.array([[[v1,h1], [v1,h6]],
                       [[v2,h1], [v2,h6]],
                       [[v3,h1], [v3,h6]],
                       [[v4,h1], [v4,h6]],
                       [[v5,h1], [v5,h6]],
                       [[v6,h1], [v6,h6]],
                       [[v7,h1], [v7,h6]],
                       [[v8,h1], [v8,h6]],
                       [[v1,h1], [v8,h1]],
                       [[v1,h2], [v8,h2]],
                       [[v1,h3], [v8,h3]],
                       [[v1,h4], [v8,h4]],
                       [[v1,h5], [v8,h5]],
                       [[v1,h6], [v8,h6]],])

alleyInterBounds = {'0':[[v2, v3],[h5, h6]],
              '1':[[v2, v3],[h3, h4]],
              '2':[[v1, v2],[h4, h5]],
              '3':[[v3, v4],[h4, h5]],
              '4':[[v4, v5],[h5, h6]],
              '5':[[v5, v6],[h4, h5]],
              '6':[[v6, v7],[h5, h6]],
              '7':[[v7, v8],[h4, h5]],
              '8':[[v6, v7],[h3, h4]],
              '9':[[v7, v8],[h2, h3]],
              '10':[[v6, v7],[h1, h2]],
              '11':[[v5, v6],[h2, h3]],
              '12':[[v4, v5],[h3, h4]],
              '13':[[v4, v5],[h1, h2]],
              '14':[[v3, v4],[h2, h3]],
              '15':[[v2, v3],[h1, h2]],
              '16':[[v1, v2],[h2, h3]],
              'A':[[v1, v2],[h5,h6]],
              'B':[[v3, v4],[h5,h6]],
              'C':[[v5, v6],[h5,h6]],
              'D':[[v7, v8],[h5,h6]],
              'E':[[v1, v2],[h3,h4]],
              'F':[[v3, v4],[h3,h4]],
              'G':[[v5, v6],[h3,h4]],
              'H':[[v7, v8],[h3,h4]],
              'I':[[v1, v2],[h1,h2]],
              'J':[[v3, v4],[h1,h2]],
              'K':[[v5, v6],[h1,h2]],
              'L':[[v7, v8],[h1,h2]]}

#intersection + alley + intersection
IAI = {"A0B":[[v1, v4],[h5, h6]],
       "E1F":[[v1, v4],[h3, h4]],
       "A2E":[[v1, v2],[h3, h6]],
       "B3F":[[v3, v4],[h3, h6]],
       "B4C":[[v3, v6],[h5, h6]],
       "C5G":[[v5, v6],[h3, h6]],
       "C6D":[[v5, v8],[h5, h6]],
       "D7H":[[v7, v8],[h3, h6]],
       "G8H":[[v5, v8],[h3, h4]],
       "H9L":[[v7, v8],[h1, h4]],
       "K10L":[[v5, v8],[h1, h2]],
       "G11K":[[v5, v6],[h1, h4]],
       "F12G":[[v3, v6],[h3, h4]],
       "J13K":[[v3, v6],[h1, h2]],
       "F14J":[[v3, v4],[h1, h4]],
       "I15J":[[v1, v4],[h1, h2]],
       "E16I":[[v1, v2],[h1, h4]]}

alleyBounds = namedtuple("alleyBounds", "cityBlocks alleyLines alleyInterBounds IAI")
R781 = alleyBounds(cityBlocks, alleyLines, alleyInterBounds, IAI)

#R808
v1, v2, v3, v4, v5, v6, v7, v8 = R808V
h1, h2, h3, h4, h5, h6 = R808H

cityBlocks = np.array([[[v2,h2], [v3,h2], [v3,h3], [v2,h3], [v2,h2]],
                       [[v2,h4], [v3,h4], [v3,h5], [v2,h5], [v2,h4]],
                       [[v4,h2], [v5,h2], [v5,h3], [v4,h3], [v4,h2]],
                       [[v4,h4], [v5,h4], [v5,h5], [v4,h5], [v4,h4]],
                       [[v6,h2], [v7,h2], [v7,h3], [v6,h3], [v6,h2]],
                       [[v6,h4], [v7,h4], [v7,h5], [v6,h5], [v6,h4]]])

alleyLines = np.array([[[v1,h1], [v1,h6]],
                       [[v2,h1], [v2,h6]],
                       [[v3,h1], [v3,h6]],
                       [[v4,h1], [v4,h6]],
                       [[v5,h1], [v5,h6]],
                       [[v6,h1], [v6,h6]],
                       [[v7,h1], [v7,h6]],
                       [[v8,h1], [v8,h6]],
                       [[v1,h1], [v8,h1]],
                       [[v1,h2], [v8,h2]],
                       [[v1,h3], [v8,h3]],
                       [[v1,h4], [v8,h4]],
                       [[v1,h5], [v8,h5]],
                       [[v1,h6], [v8,h6]],])

alleyInterBounds = {'0':[[v2, v3],[h5, h6]],
              '1':[[v2, v3],[h3, h4]],
              '2':[[v1, v2],[h4, h5]],
              '3':[[v3, v4],[h4, h5]],
              '4':[[v4, v5],[h5, h6]],
              '5':[[v5, v6],[h4, h5]],
              '6':[[v6, v7],[h5, h6]],
              '7':[[v7, v8],[h4, h5]],
              '8':[[v6, v7],[h3, h4]],
              '9':[[v7, v8],[h2, h3]],
              '10':[[v6, v7],[h1, h2]],
              '11':[[v5, v6],[h2, h3]],
              '12':[[v4, v5],[h3, h4]],
              '13':[[v4, v5],[h1, h2]],
              '14':[[v3, v4],[h2, h3]],
              '15':[[v2, v3],[h1, h2]],
              '16':[[v1, v2],[h2, h3]],
              'A':[[v1, v2],[h5,h6]],
              'B':[[v3, v4],[h5,h6]],
              'C':[[v5, v6],[h5,h6]],
              'D':[[v7, v8],[h5,h6]],
              'E':[[v1, v2],[h3,h4]],
              'F':[[v3, v4],[h3,h4]],
              'G':[[v5, v6],[h3,h4]],
              'H':[[v7, v8],[h3,h4]],
              'I':[[v1, v2],[h1,h2]],
              'J':[[v3, v4],[h1,h2]],
              'K':[[v5, v6],[h1,h2]],
              'L':[[v7, v8],[h1,h2]]}

#intersection + alley + intersection
IAI = {"A0B":[[v1, v4],[h5, h6]],
       "E1F":[[v1, v4],[h3, h4]],
       "A2E":[[v1, v2],[h3, h6]],
       "B3F":[[v3, v4],[h3, h6]],
       "B4C":[[v3, v6],[h5, h6]],
       "C5G":[[v5, v6],[h3, h6]],
       "C6D":[[v5, v8],[h5, h6]],
       "D7H":[[v7, v8],[h3, h6]],
       "G8H":[[v5, v8],[h3, h4]],
       "H9L":[[v7, v8],[h1, h4]],
       "K10L":[[v5, v8],[h1, h2]],
       "G11K":[[v5, v6],[h1, h4]],
       "F12G":[[v3, v6],[h3, h4]],
       "J13K":[[v3, v6],[h1, h2]],
       "F14J":[[v3, v4],[h1, h4]],
       "I15J":[[v1, v4],[h1, h2]],
       "E16I":[[v1, v2],[h1, h4]]}

R808 = alleyBounds(cityBlocks, alleyLines, alleyInterBounds, IAI)

#R859
v1, v2, v3, v4, v5, v6, v7, v8 = R859V
h1, h2, h3, h4, h5, h6 = R859H

cityBlocks = np.array([[[v2,h2], [v3,h2], [v3,h3], [v2,h3], [v2,h2]],
                       [[v2,h4], [v3,h4], [v3,h5], [v2,h5], [v2,h4]],
                       [[v4,h2], [v5,h2], [v5,h3], [v4,h3], [v4,h2]],
                       [[v4,h4], [v5,h4], [v5,h5], [v4,h5], [v4,h4]],
                       [[v6,h2], [v7,h2], [v7,h3], [v6,h3], [v6,h2]],
                       [[v6,h4], [v7,h4], [v7,h5], [v6,h5], [v6,h4]]])

alleyLines = np.array([[[v1,h1], [v1,h6]],
                       [[v2,h1], [v2,h6]],
                       [[v3,h1], [v3,h6]],
                       [[v4,h1], [v4,h6]],
                       [[v5,h1], [v5,h6]],
                       [[v6,h1], [v6,h6]],
                       [[v7,h1], [v7,h6]],
                       [[v8,h1], [v8,h6]],
                       [[v1,h1], [v8,h1]],
                       [[v1,h2], [v8,h2]],
                       [[v1,h3], [v8,h3]],
                       [[v1,h4], [v8,h4]],
                       [[v1,h5], [v8,h5]],
                       [[v1,h6], [v8,h6]],])

alleyInterBounds = {'0':[[v2, v3],[h5, h6]],
              '1':[[v2, v3],[h3, h4]],
              '2':[[v1, v2],[h4, h5]],
              '3':[[v3, v4],[h4, h5]],
              '4':[[v4, v5],[h5, h6]],
              '5':[[v5, v6],[h4, h5]],
              '6':[[v6, v7],[h5, h6]],
              '7':[[v7, v8],[h4, h5]],
              '8':[[v6, v7],[h3, h4]],
              '9':[[v7, v8],[h2, h3]],
              '10':[[v6, v7],[h1, h2]],
              '11':[[v5, v6],[h2, h3]],
              '12':[[v4, v5],[h3, h4]],
              '13':[[v4, v5],[h1, h2]],
              '14':[[v3, v4],[h2, h3]],
              '15':[[v2, v3],[h1, h2]],
              '16':[[v1, v2],[h2, h3]],
              'A':[[v1, v2],[h5,h6]],
              'B':[[v3, v4],[h5,h6]],
              'C':[[v5, v6],[h5,h6]],
              'D':[[v7, v8],[h5,h6]],
              'E':[[v1, v2],[h3,h4]],
              'F':[[v3, v4],[h3,h4]],
              'G':[[v5, v6],[h3,h4]],
              'H':[[v7, v8],[h3,h4]],
              'I':[[v1, v2],[h1,h2]],
              'J':[[v3, v4],[h1,h2]],
              'K':[[v5, v6],[h1,h2]],
              'L':[[v7, v8],[h1,h2]]}

#intersection + alley + intersection
IAI = {"A0B":[[v1, v4],[h5, h6]],
       "E1F":[[v1, v4],[h3, h4]],
       "A2E":[[v1, v2],[h3, h6]],
       "B3F":[[v3, v4],[h3, h6]],
       "B4C":[[v3, v6],[h5, h6]],
       "C5G":[[v5, v6],[h3, h6]],
       "C6D":[[v5, v8],[h5, h6]],
       "D7H":[[v7, v8],[h3, h6]],
       "G8H":[[v5, v8],[h3, h4]],
       "H9L":[[v7, v8],[h1, h4]],
       "K10L":[[v5, v8],[h1, h2]],
       "G11K":[[v5, v6],[h1, h4]],
       "F12G":[[v3, v6],[h3, h4]],
       "J13K":[[v3, v6],[h1, h2]],
       "F14J":[[v3, v4],[h1, h4]],
       "I15J":[[v1, v4],[h1, h2]],
       "E16I":[[v1, v2],[h1, h4]]}

R859 = alleyBounds(cityBlocks, alleyLines, alleyInterBounds, IAI)


alleyInterType = {'A':['intersection'],
                  'B':['intersection'],
                  'C':['intersection'],
                  'D':['intersection'],
                  'E':['intersection'],
                  'F':['intersection'],
                  'G':['intersection'],
                  'H':['intersection'],
                  'I':['intersection'],
                  'J':['intersection'],
                  'K':['intersection'],
                  'L':['intersection'],
                  '0':['alley','horizontal'],
                  '1':['alley','horizontal'],
                  '4':['alley','horizontal'],
                  '6':['alley','horizontal'],
                  '8':['alley','horizontal'],
                  '10':['alley','horizontal'],
                  '12':['alley','horizontal'],
                  '13':['alley','horizontal'],
                  '15':['alley','horizontal'],
                  '2':['alley','vertical'],
                  '3':['alley','vertical'],
                  '5':['alley','vertical'],
                  '7':['alley','vertical'],
                  '9':['alley','vertical'],
                  '11':['alley','vertical'],
                  '14':['alley','vertical'],
                  '16':['alley','vertical']}