#! /usr/bin/env python3

import numpy as np
from scipy.special import k0, i0
from scipy.constants import k, c, h
from math import pi

import matplotlib.pyplot as plt

MHz = 1e6

fr = 310*MHz
Tc = 1.4
delta = 1.76*k*Tc
alphak = 0.51

T = np.r_[90:400:1000j]*1e-3

q = h*fr/(2*k*T)
xMB_full = -alphak*np.sqrt(pi*k*T/2/delta)*np.exp(-delta/k/T)
xMB_full *= (1 + np.sqrt(2*delta/(pi*k*T)))*np.exp(-q)*i0(q)
xMB_full *= 1e6

xMB_approx = -alphak*np.sqrt(pi*k*T/2/delta)*(1 +
		1/pi*np.sqrt(2*delta/h/fr))*np.exp(-delta/k/T)
xMB_approx *= 1e6


dQiMB_full = 4*alphak/pi*np.exp(-delta/k/T)*np.sinh(q)*k0(q)
dQiMB_approx = 2*alphak*np.sqrt(h*fr/pi/k/T)*np.exp(-delta/k/T)
dQiloss = 1e-5
dQiTLS = 5e-5*np.tanh(q)

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(T*1e3, xMB_full, 'k', label='full')
ax.plot(T*1e3, xMB_approx, 'k--', label='approx')
ax.set_xlabel('Temperature [mK]')
ax.legend(loc='upper right')
ax.set_ylabel('Frequency Shift [ppm]')
ax.grid(which='both')



fig, ax = plt.subplots(figsize=(10,10))
ax.semilogy(T*1e3, dQiMB_full + dQiTLS, 'k', label='full')
ax.semilogy(T*1e3, dQiMB_approx + dQiTLS, 'k--', label='approx')
ax.set_xlabel('Temperature [mK]')
ax.set_ylabel('1/Qi_MB')
ax.legend(loc='upper right')
ax.grid(which='both')
plt.show()
