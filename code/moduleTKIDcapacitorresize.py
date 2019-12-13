import numpy as np
import matplotlib.pyplot as plt
from math import pi

N0 = 502
h0 = 1006

pF = 1e-12
nH = 1e-9
MHz = 1e6

N = np.r_[500:1000:100j]
h = np.r_[100:1000:100j]

X, Y = np.meshgrid(N, h)

C = 25.10*(X/N0 + Y/h0 - 1) + 0.62
Lpar = 1.54*(X/N0)**2 + 0.50*(X/N0) + 2.00*(Y/h0) - 0.22
L0 = 10
Ltotal = L0 + Lpar

fr = 1./(2*pi)/np.sqrt(Ltotal*nH*C*pF)/MHz
print (fr[::10, ::10])

fig, ax = plt.subplots(figsize=(10,10))
#img = ax.pcolormesh(X, Y, fr, alpha=0.5)
img2 = ax.contour(X, Y, fr)
ax.set_xlabel('N')
ax.set_ylabel('h [um]')
#fig.colorbar(img)
plt.savefig('frequency_vs_N_vs_h.png')

