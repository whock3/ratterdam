# -*- coding: utf-8 -*-
"""

"""

from RDP_vectorized import RDP
import matplotlib.pyplot as plt
import numpy as np
import csv
with open("pos.p.ascii", 'r') as csv_file:
    data_iter = csv.reader(csv_file)
    data = [data for data in data_iter]
data_np1 = np.array(data[50574:-17], dtype = "float64")
a = np.array([1, 4.72, 4.72])
data_np2 = data_np1[:,0:3]/a[None,:]
data_np3 = data_np2[np.all(data_np2 > np.array([0, 0, 0]), axis=1)]

count = 1

for i in np.linspace(0.1, 2.3, 10):
    var_sum = 0
    for j in np.arange(3000, 98800, 3000):
        ResultList = data_np3[j-3000]
        var_sum = var_sum + np.var(RDP(data_np3[j-3000:j], i))
        
    if count == 1:
        var = [i, var_sum/32]
    else:
        var = np.vstack([var, [i, var_sum/32]])
    count += 1

plt.plot(data_np3[23000:23300, 0], data_np3[23000:23300, 1])
plt.show()
RList = RDP(data_np3[23000:23300], 1)
plt.plot(RList[:,0], RList[:,1])
plt.show()