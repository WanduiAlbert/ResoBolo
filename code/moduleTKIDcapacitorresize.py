import numpy as np
import matplotlib.pyplot as plt
from math import pi

N0 = 502
h0 = 1006
h1 = 550
N1 = 650

pF = 1e-12
nH = 1e-9
MHz = 1e6

N = np.r_[300:1000:500j]
h = np.r_[200:1200:500j]

X, Y = np.meshgrid(N, h)

C = 25.10*(X/N0 + Y/h0 - 1) + 0.62
Lpar = 1.54*(X/N0)**2 + 0.50*(X/N0) + 2.00*(Y/h0) - 0.22
L0 = 10
Ltotal = L0 + Lpar

fr = 1./(2*pi)/np.sqrt(Ltotal*nH*C*pF)/MHz
fr[fr > 700] = np.nan
#print (fr[::10, ::10])

fig, ax = plt.subplots(figsize=(10,10))
img = ax.pcolormesh(X, Y, fr, cmap='RdBu', alpha=0.5)
contours = ax.contour(X, Y, fr, 20, colors='black')
plt.clabel(contours, inline=True, fontsize=18)
ax.axhline(h0, color='w', ls='--')
ax.axhline(h1, color='w', ls='--')
ax.axvline(N0, color='w', ls='--')
ax.set_xlabel('N')
ax.set_ylabel('h [um]')
fig.colorbar(img, label='Frequency [MHz]')
plt.savefig('frequency_vs_N_vs_h.png')
plt.show()


C1 = 25.10*(N1/N0 + h1/h0 - 1) + 0.62
Lpar1 = 1.54*(N1/N0)**2 + 0.50*(N1/N0) + 2.00*(h1/h0) - 0.22

print (C1, Lpar1)
