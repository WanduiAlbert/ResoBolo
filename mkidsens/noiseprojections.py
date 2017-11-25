#! /usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from tkid_params import *

# First thing is I want to plot out the responsivity of the bolometer as a
# function of frequency

print ("The bolometer time constant", tau_b)

nu = np.r_[0:1000:1000j] * u.Hz

r = ((chi_c * chi_qp * beta/4) * tau_qp/tau_th * kappa/P_b * 1/(1 + 1j *
    2*np.pi*nu*tau_b)).to(1/u.pW)

S_opt = ((2 * (1 + n_opt) * h * nu_opt * P_opt * np.ones_like(r.value,
    dtype=np.float64))**0.5).to(u.aW/u.Hz**0.5)

S_ph = ((4 * k_B * T_b**2 * G_b * np.ones_like(r.value,
  dtype=np.float64))**0.5).to(u.aW/u.Hz**0.5)

S_amp = (((k_B * T_amp/P_g)/np.abs(r)**2)**0.5).to(u.aW/u.Hz**0.5)

S_shot = 2 * n_qp * V_sc * (1/tau_max + 1/tau_qp)

S_gr = ((((chi_c*chi_qp*beta/4)*(tau_qp/(n_qp * V_sc)))**2 *
  S_shot/np.abs(r)**2)**0.5).to(u.aW/u.Hz**0.5)

fig, ax = plt.subplots(figsize=(10,10))
ax.loglog(nu, S_opt, 'k--', label='Optical')
ax.loglog(nu, S_ph, 'b', label='Phonon')
ax.loglog(nu, S_amp, 'r--', label='Amplifier')
ax.loglog(nu, S_gr, 'g', label='Gen-Recomb')

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'NEP [aW/rtHz]')
ax.grid(which='both')
ax.legend(loc='best')

plt.savefig(r'NEP.png')
plt.show()

exit()

fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'Responsivity [Hz/fW]')
ax.legend(loc = 'best')
ax.grid(which = 'both')
ax.axis('tight')

fig2, ax2 = plt.subplots(figsize=(10,10))
ax3 = ax2.twinx()
ax3.semilogx(nu, np.angle(r), 'r', label=r'phase')
ax2.semilogx(nu, np.abs(r), 'k', label=r'abs')
ax2.set_xlabel(r'Frequency [Hz]')
ax2.set_ylabel(r'Responsivity [Hz/fW]')
ax2.grid(which = 'both')
ax2.axis('tight')
ax3.set_ylabel(r'Phase')
plt.show()
