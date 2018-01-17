#! /usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from tkid_params import *

# First thing is I want to plot out the responsivity of the bolometer as a
# function of frequency

print ("The bolometer time constant", tau_b)


nu = np.r_[0:1000:5000j] * u.Hz

# I'll ignore the bolometer roll off factor for now
H = np.ones_like(nu.value)
#H =  1/(1 + 1j * 2 * np.pi * nu * tau_b)
#H = np.abs(H)

# r = dx/dPopt
r = ((chi_qp * beta/2/Q_i) * tau_qp/tau_th * kappa/P_b).to(1/u.pW) * H

# S_opt = (2 * (1 + n_opt) * h * nu_opt * P_opt * H).to(u.aW**2/u.Hz)

S_ph = (4 * chi_ph * k_B * T_b**2 * G_b * H ).to(u.aW**2/u.Hz)

NEP_ph = S_ph ** 0.5

print ("The generated power is {0:2.2f}".format(P_g))

S_amp = (k_B * T_amp/P_g).to(1/u.Hz)

NEP_amp = (S_amp**0.5/r).to(u.aW/u.Hz**0.5) * H


S_shot = 2 * n_qp * V_sc * (1/tau_max + 1/tau_qp) * H

NEP_gr = (S_shot**0.5 * (tau_th/(n_qp * V_sc) * P_b/kappa)).to(u.aW/u.Hz**0.5)

NEP_total = (NEP_gr**2 + NEP_amp**2 + NEP_ph**2)**0.5
# S_gr = ((((chi_c*chi_qp*beta/4)*(tau_qp/(n_qp * V_sc)))**2 *
#   S_shot/np.abs(r)**2)**0.5).to(u.aW/u.Hz**0.5)

r_f = (r * f_r).to(u.kHz/u.pW)

print(r_f[0])

fig, ax = plt.subplots(figsize=(10,10))
# ax.loglog(nu, S_opt, 'k--', label='Optical')
ax.loglog(nu, NEP_ph, 'b', label='Phonon')
ax.loglog(nu, NEP_amp, 'r--', label='Amplifier')
ax.loglog(nu, NEP_gr, 'g', label='Gen-Recomb')
ax.loglog(nu, NEP_total, 'k', label='Total')

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'NEP [aW/rtHz]')
# ax.set_ylim([10, 1000])
ax.set_xlim([0.1, 1000])
ax.grid(which='major', axis='both')
ax.legend(loc='best')

plt.savefig(r'NEP.png')
plt.show()

# exit()


# r *= (f_r * Q_c/2/Q_r**2)
# r.to(u.kHz/u.pW)
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'Frequency [Hz]')
# ax.set_ylabel(r'Responsivity [Hz/fW]')
# ax.legend(loc = 'best')
# ax.grid(which = 'both')
# ax.axis('tight')

# fig2, ax2 = plt.subplots(figsize=(10,10))
# # ax3 = ax2.twinx()
# # ax3.semilogx(nu, np.angle(r), 'r', label=r'phase')
# ax2.semilogx(nu, np.abs(r_f), 'k', label=r'abs')
# ax2.set_xlabel(r'Frequency [Hz]')
# ax2.set_ylabel(r'Responsivity [kHz/pW]')
# ax2.grid(which = 'both')
# ax2.axis('tight')
# # ax3.set_ylabel(r'Phase')
# plt.show()
