import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

w = 2
g = 1.5 * np.arange(1,7)
s = 2.0 * np.arange(1, 7)

Npair = 20
L = w + Npair * 2 * (s + w)
W = 76.5

Nsq = L * (W - g)/w/(s + w)

C = np.array([0.1656, 0.1341, 0.1168, 0.1088, 0.1025, 0.09805])

Cs = C/Nsq * 1e4 #to convert to pF


fig, ax = plt.subplots(figsize=(12,12))
ax.plot(s, Cs, 'bs-')
ax.set_xlabel(r'Spacing s [$\mu m$]')
ax.set_ylabel(r'Surface Capacitance ($ \times 10^{-4} $) [pF/sq]')
ax.grid(which='both')
ax.axis('tight')


fig2, ax2 = plt.subplots(figsize=(12,12))
ax2.plot(g, Cs, 'bs-')
ax2.set_xlabel(r'End Gap g [$\mu m$]')
ax2.set_ylabel(r'Surface Capacitance ($ \times 10^{-4} $) [pF/sq]')
ax2.grid(which='both')
ax2.axis('tight')

plt.show()