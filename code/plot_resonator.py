#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi

MHz = 1e6
kHz = 1e3

Qr = 10000
Qc = 20000
phic = 0.3
Qe = Qc*np.cos(phic)*np.exp(1j*phic)
fr = 336*MHz

df = np.r_[-0.1:0.1:1000j]*MHz

x = df/fr

dQe = 1./Qe
dQr = 1./Qr
dQc = 1./Qc

rho1 = dQc/2/dQr
rho2 = np.abs(dQe)/2/dQr

S21 = 1. - dQe/(dQr + 2.j*x)
S21dB = 20*np.log10(np.abs(S21))

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(df/kHz, np.abs(S21), 'b')
#ax.plot(df/kHz, np.real(S21), 'r')
#ax.plot(df/kHz, np.imag(S21), 'g')
ax.grid()
ax.set_xlabel('$f - f_r$ [kHz]')
ax.set_ylabel('$|S_{21}|$')
#text = '$Q_r$ = %d\n$Q_c$ = %d\n$\phi_c$ = %1.1f'%(Qr,Qc,phic)
#ax.annotate(text, (-100,-2), (-100,-3))
plt.savefig('S21_mag.png')


fig, ax = plt.subplots(figsize=(10,10))
x = np.r_[1-2*rho2:1.05:100j]
y = -np.tan(phic)*x + (1+0*rho1)*np.tan(phic)
ax.plot(S21.real, S21.imag, 'b', linewidth=2)
ax.scatter(1-rho2*np.cos(phic), rho2*np.sin(phic), color='k', s=40)
ax.plot(x,y,linestyle='--', color='k')
ax.hlines(0, 1-2*rho2,1.05, colors='k',linestyles='--')
ax.axis('square')
ax.grid()
ax.set_xlabel('I')
ax.set_ylabel('Q')
#text = '$Q_r$ = %d\n$Q_c$ = %d\n$\phi_c$ = %1.1f'%(Qr,Qc,phic)
#ax.annotate(text, (0.6,0.1), (0.6,0))
plt.savefig('IQ.png')
