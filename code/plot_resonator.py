#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi

MHz = 1e6
kHz = 1e3

Qr = 20000
Qc = 40000
phic = 0.5
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

S21sym = 1. - dQc/(dQr + 2.j*x)
S21symdB = 20*np.log10(np.abs(S21sym))

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(df/kHz, np.abs(S21sym), 'r', label='\n$\phi_c$ = 0')
ax.plot(df/kHz, np.abs(S21), 'b', label='\n$\phi_c$ = %1.1f'%phic)
#ax.plot(df/kHz, np.real(S21), 'r')
#ax.plot(df/kHz, np.imag(S21), 'g')
ax.grid()
ax.set_xlabel('$f - f_r$ [kHz]')
ax.set_ylabel('$|S_{21}|$')
ax.legend(loc='lower left')
text = '$Q_r$ = %d\n$Q_c$ = %d'%(Qr,Qc)
ax.annotate(text, (-100, 0.75), (-100, 0.73), fontsize=25)
plt.savefig('S21_mag.png')
#plt.show()

fig, ax = plt.subplots(figsize=(10,10))
x = np.r_[1-2*rho1:1.0:100j]
y = -np.tan(phic)*x + np.tan(phic)
ax.plot(S21sym.real, S21sym.imag, 'r', linewidth=2, label='\n$\phi_c$ = 0')
ax.plot(S21.real, S21.imag, 'b', linewidth=2, label='\n$\phi_c$ = %1.1f'%phic)
ax.plot(x,y,linestyle='-', color='k')
ax.hlines(0, 1-2*rho1, 1.0, colors='k',linestyles='-')
ax.scatter(1-rho2*np.cos(phic), rho2*np.sin(phic), color='k', s=40)
ax.scatter(1-rho1, 0, color='k', s=40)
ax.scatter(1-2*rho1, 0, color='r', s=60)
ax.scatter(1-2*rho1, 2*rho1*np.tan(phic), color='b', s=60)
text = '$\phi_c$'
ax.annotate(text, (0.9, 0.0), (0.9, 0.01), fontsize=25)
text = '$f\ =\ f_r$'
ax.annotate(text, (1-2*rho1, 0.0),
        (1-2*rho1 + 0.08, 0.08), fontsize=25, color='r',
        arrowprops={'width':2, 'color':'r'})
ax.annotate(text, (1-2*rho1, 2*rho1*np.tan(phic)),
        (1-2*rho1 + 0.1, 2*rho1*np.tan(phic) + 0.06), fontsize=25, color='b',
        arrowprops={'width':2, 'color':'b'})
text = r'$Q_r/Q_c$'
ax.annotate(text, (1-rho1, 0.0), (1-rho1, -0.03), fontsize=25)
text = r'$Q_r/\left|\hat{Q}_c\right|$'
ax.annotate(text, (1-rho2*np.cos(phic), rho2*np.sin(phic)),
        (1-rho2*np.cos(phic), rho2*np.sin(phic) + 0.02), fontsize=25)
ax.legend(loc='lower left')
ax.axis('square')
ax.grid()
ax.set_xlabel('I')
ax.set_ylabel('Q')
#text = '$Q_r$ = %d\n$Q_c$ = %d\n$\phi_c$ = %1.1f'%(Qr,Qc,phic)
#ax.annotate(text, (0.6,0.1), (0.6,0))
plt.savefig('IQ.png')
plt.show()
