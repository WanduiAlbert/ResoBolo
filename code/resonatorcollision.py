#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from itertools import permutations, combinations

measured_freqs = np.array([299.8, 308.7, 315.7 , 316.6, 339.8, 342.8, 397.4,
    400.6, 404.3, 411.37, 411.45, 412.99, 416.3, 421.8, 424.6, 441.9, 445.7,
    460.4, 477.189, 477.198, 478.9, 485.2, 493.3, 495.4])
#measured_freqs = np.array([299.8, 308.7, 315.7 , 316.6, 339.8, 342.8, 397.4,
#    400.6, 404.3, 411.4, 412.9, 416.3, 421.8, 423.2, 424.6, 441.9, 445.7,455.4,
#    460.4, 477.19, 478.9, 485.2, 493.3, 495.4])

measured_Qcs = np.array([26155, 38302, 39604, 273236, 308894, 117245, 85398,
    227683, 106499, 55007, 62470, 153821, 306648, 79960, 315739, 16307, 12628,
    64543, 1716716, 29189, 49033, 68830, 86019, 23951])
#measured_Qcs = np.array([26155, 38302, 39604, 273236, 308894, 117245, 85398,
#    227683, 106499, 55007, 62470, 153821, 306648, 79960, 315739, 16307, 12628,
#    64543, 29189, 49033, 68830, 86019, 23951])

true_freqs = np.array([306, 318, 321, 324, 327, 330, 333, 336, 339, 351, 354,
    369, 372, 381, 384, 387, 390, 393, 396, 399, 402, 405, 408, 417, 420, 435,
    438, 450, 453, 456, 459, 462, 465, 468, 471, 483])


#plt.hist(np.diff(true_freqs), histtype='step', linewidth=2, edgecolor='b',
#        density=True, label="Predicted")
#plt.hist(np.diff(measured_freqs), histtype='step', linewidth=2, edgecolor='r',
#        density=True, label="Measured")
#plt.grid()
#plt.xlabel("Frequency Shift")
#plt.ylabel("Normalized Counts")
#plt.legend()
#plt.show()

index = np.arange(true_freqs.size)
sampler = combinations(index, measured_freqs.size)

N = 4000000
diff = np.zeros(N)
label = np.arange(N)
for i, samp in enumerate(sampler):
    print (i)
    if i >= N: break
    samp_freqs = true_freqs[list(samp)]
    diff[i] = np.std((samp_freqs - measured_freqs))

plt.plot(label, diff)
plt.show()
