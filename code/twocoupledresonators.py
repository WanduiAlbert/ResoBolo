#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import scipy.linalg as linalg
import reso_fit

N = 10

nH = 1e-9
pF = 1e-12
MHz = 1e6

L = 3*nH
L12 = 18*nH
C = 40*pF
Ca = 0.2*pF
Cb = 0.21*pF
R = 1e6
R1 = 1e6
R2 = 1e6

Z0 = 50

f = np.r_[400:600:920000j]*MHz
omega = 2*pi*f

numS21  = 2j*(2*L+L12)*R1*R2
numS21 += -2*L*(L+L12)*(R1+R2)*omega
numS21 += - 2j*L*(L*L12 + (2*C + Ca + Cb)*(L+L12)*R1*R2)*omega**2
numS21 += 2*L**2*L12*((C + Ca)*R1 + (C + Cb)*R2)*omega**3
numS21 += 2j*(Ca + C)*(C + Cb)*L**2*L12*R1*R2*omega**4

denomS21  = -2j*(2*L + L12)*R1*R2
denomS21 += (2*L*(L+L12)*(R1+R2) + (Ca + Cb)*(2*L+L12)*R1*R2*Z0)*omega
denomS21 += 1j*L*(2*L*L12 + 2*(2*C+Ca+Cb)*(L+L12)*R1*R2 +
        (Ca+Cb)*(L+L12)*(R1+R2)*Z0)*omega**2
denomS21 += -L*(2*L*L12*((C+Ca)*R1 + (C+Cb)*R2) + ((Ca+Cb)*L*L12 + 2*(Ca*Cb*L12
    + C*(Ca+Cb)*(L+L12))*R1*R2)*Z0)*omega**3
denomS21 += -1j*L**2*L12*(2*(C+Ca)*(C+Cb)*R1*R2 + (Ca*Cb + C*(Ca+Cb))*(R1+R2)*Z0)*omega**4
denomS21 += C*(2*Ca*Cb + C*(Ca+Cb))*L**2*L12*R1*R2*Z0*omega**5

S21 = numS21/denomS21
re = S21.real
im = S21.imag

nt = f.size//2
feven = f[np.argmin(np.abs(S21)[:nt])]
fodd = f[nt + np.argmin(np.abs(S21)[nt:])]
feven = 1./np.sqrt(C*L)/2/pi
fodd = np.sqrt((1/L + 2./L12)/C)/2/pi

print (feven/MHz, fodd/MHz)

mask_evn = np.abs((f - feven)/MHz) < 2
f0_evn,Qi_evn,Qr_evn, Qe_evn_re, Qe_evn_im, a_evn,ymodel_evn_re,ymodel_evn_im = reso_fit.do_fit(f[mask_evn], re[mask_evn],im[mask_evn],plot=False,get_cov=False,verbose=False)
Qe_evn = Qe_evn_re + 1j*Qe_evn_im
dQe_evn = 1./Qe_evn
Qc_evn = 1./np.real(dQe_evn)
phic_evn = np.arctan2(Qe_evn_im, Qe_evn_re)
ffine_evn = np.r_[f[mask_evn][0]:f[mask_evn][-1]:10000j]
ymodel_evn = np.sqrt(ymodel_evn_re**2 + ymodel_evn_im**2)

mask_odd = np.abs((f - fodd)/MHz) < 2
f0_odd,Qi_odd,Qr_odd, Qe_odd_re, Qe_odd_im, a_odd,ymodel_odd_re,ymodel_odd_im = reso_fit.do_fit(f[mask_odd], re[mask_odd],im[mask_odd],plot=False,get_cov=False,verbose=False)
Qe_odd = Qe_odd_re + 1j*Qe_odd_im
dQe_odd = 1./Qe_odd
Qc_odd = 1./np.real(dQe_odd)
phic_odd = np.arctan2(Qe_odd_im, Qe_odd_re)
ffine_odd = np.r_[f[mask_odd][0]:f[mask_odd][-1]:10000j]
ymodel_odd = np.sqrt(ymodel_odd_re**2 + ymodel_odd_im**2)

print (f0_evn, Qc_evn, phic_evn)
print (f0_odd, Qc_odd, phic_odd)

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(f/MHz, np.abs(S21))
ax.plot(ffine_evn/MHz, ymodel_evn, 'r')
ax.plot(ffine_odd/MHz, ymodel_odd, 'r')
ax.vlines(f0_evn, 0, 1, color='k', ls='dashed')
ax.vlines(f0_odd, 0, 1, color='k', ls='dashed')
ax.grid()
ax.set_xlabel('Frequency [MHz]')
ax.set_ylabel('|S21|')
plt.show()



