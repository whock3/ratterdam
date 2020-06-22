# -*- coding: utf-8 -*-
"""

"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def rect(bottom_left, top_right):
    return patches.Rectangle(tuple(bottom_left), (top_right[0]-bottom_left[0]), (top_right[1]-bottom_left[1]), fill=False)

fig,ax=plt.subplots(1)
ax.add_patch(rect([145, 250], [332, 415]))
ax.add_patch(rect([160, 235], [200, 292]))
ax.add_patch(rect([97, 174], [260, 365]))
ax.add_patch(rect([220, 310], [275, 350]))
ax.add_patch(rect([285, 380], [332, 415]))
ax.add_patch(rect([355, 445], [275, 350]))
ax.add_patch(rect([420, 518], [332, 415]))
ax.add_patch(rect([485, 556], [260, 365]))
ax.add_patch(rect([435, 503], [200, 292]))
ax.add_patch(rect([490, 560], [135, 230]))
ax.add_patch(rect([420, 518], [87, 158]))
ax.add_patch(rect([355, 445], [150, 215]))
ax.add_patch(rect([300, 365], [200, 292]))
ax.add_patch(rect([300, 365], [87, 158]))
ax.add_patch(rect([220, 310], [150, 215]))
ax.add_patch(rect([145, 250], [87, 158]))
ax.add_patch(rect([97, 174], [135, 230]))
ax.set_xlim(0, 50)
ax.set_ylim(0, 60)
