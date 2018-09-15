#! /usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi

MHz = 1e6
Qr = 20000
fr = 300 * MHz
f = np.r_[295:305:10000j]*MHz
x0 = (f-fr)/fr
y0 = Qr*x0

# Now solve the cubic equation to determine the resonator detuning
a =-12 
k2 = np.sqrt((y0**3/27 + y0/12 + a/8)**2 - (y0**2/9 - 1/12)**3, dtype=np.complex128)
k1 = np.power(a/8 + y0/12 + k2 + y0**3/27, 1./3)
eps = (-1 + 3**0.5 * 1j)/2

y1 = y0/3 + (y0**2/9-1/12)/k1 + k1
y2 = y0/3 + (y0**2/9-1/12)/eps/k1 + eps*k1
y3 = y0/3 + (y0**2/9-1/12)/eps**2/k1 + eps**2*k1

#thresh = 1e-4
#low_to_high = False
#if low_to_high:
#    y = y2.real
#    mask = (np.abs(y2.imag) >= thresh)
#    y[mask] = y1.real[mask]
#else:
#    y = y1.real
#    mask = (np.abs(y1.imag) >= thresh)
#    y[mask] = y2.real[mask]

#plt.figure()
#plt.plot(y0, y)
#plt.grid()
#plt.axis('tight')
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(y0, y1.real, label='y1 Real')
plt.plot(y0, y2.real, label='y2 Real')
plt.plot(y0, y3.real, label='y3 Real')
plt.plot(y0, y1.imag,ls="dashed", label='y1 Imag')
plt.plot(y0, y2.imag,ls="dashed", label='y2 Imag')
plt.plot(y0, y3.imag,ls="dashed", label='y3 Imag')
plt.xlabel('y0 [Linewidths]')
plt.ylabel('y [Linewidths]')
plt.xlim([-15,15])
plt.ylim([-15,15])
plt.legend()
plt.grid()
#plt.axis('tight')
#plt.savefig('cubic_roots.png')
plt.show()
