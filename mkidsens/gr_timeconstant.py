#! /usr/bin/env python3

"""
My goal is to estimate the relative sizes of the various power transport
mechanisms in our devices to see if they are large enough to be considered.
The three key terms are
1. Electron-phonon coupling
2. Kapitza boundary coupling
3. Bolometric coupling

"""
import numpy as np
from math import pi
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy.constants import h, k, c
from scipy.special import kn, iv
from mpl_toolkits.axes_grid1 import make_axes_locatable

np.set_printoptions(precision=4)

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)
hbar = h/2/pi
MHz = 1e6
kHz = 1e3
pW = 1e-12
aW = 1e-18
um  = 1e-6
microsecond = 1e-6
mJ = 1e-3
mol = 1
Kelvin = 1
cm = 1e-2
g = 1
Ohm = 1
# Material properties of the Aluminum superconductor
gamma = 1.35 * mJ/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g/mol
rho = 1.15 * Ohm * cm

Tc = 1.329
Delta = 1.764*k*Tc
alphak = 0.4

PdBm = -90
Preadout = 1e-3*10**(PdBm/10)
eta_read = 0.00
Gamma_gen = eta_read*Preadout/Delta

print (Preadout/pW)
Tbath = 80e-3
K = 140*pW #/K^3
Nsq = 16200
w = 1*um
t = 0.05*um
A = t*w
l = Nsq*w
Vsc = l * A


# Want to fix T and then play around with n_qp_star and tau_max
#T = np.r_[90:800:100j]*1e-3
T = 320e-3
N_0 = ((3 * (gamma * (density/A_r)))/(2*pi**2 * k**2))
n_th = 2 * N_0 * np.sqrt(2*pi* k * T* Delta)*np.exp(-Delta/(k*T))
#tau_max = 955.1 * microsecond
tau_max = np.r_[600:1200:200j] * microsecond
#n_qp_star = 378 * 1/um**3
n_qp_star = np.r_[300:800:200j] * 1/um**3
Tmax, Nqp_star = np.meshgrid(tau_max, n_qp_star)
n_qp = np.sqrt((n_th + Nqp_star)**2 + 2*Gamma_gen*Nqp_star*Tmax/Vsc) - Nqp_star
tau = Tmax/(1 + n_qp/Nqp_star)
f3dB = 1/2/pi/tau

GT = 13.6*pW
kappa = 0.5 + Delta/(k*T)
NEPgr = 2*GT/kappa * np.sqrt(tau/(n_qp*Vsc))


limits = [tau_max[0]/microsecond, tau_max[-1]/microsecond, n_qp_star[0]*um**3,
        n_qp_star[-1]*um**3]
plt.figure(figsize=(12,10))
plt.ylabel('n_qp_star [um^-3]')
plt.xlabel('tau_max [us]')
plt.title("NEP and f3dB at T=%1.1fmK"%(T*1e3))
CS = plt.contour(Tmax/microsecond, Nqp_star*um**3, NEPgr/aW, 30, colors='black')
plt.clabel(CS, inline=True, fontsize=10)
im = plt.imshow(f3dB/kHz, origin='lower', cmap='viridis', alpha=0.8, extent=limits)
plt.scatter([955.1], [359], s=25, color='k', marker='X')
ax=plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
#cbar.set_label('NEP [aW/rtHz]')#, rotation=270)
cbar.set_label('f3dB [kHz]')#, rotation=270)


#fig2, ax2 = plt.subplots(figsize=(10,10))
#CS2 = ax.contourf(Tmax/microsecond, Nqp_star*um**3, NEPgr/aW, 30)
#ax2.clabel(CS2, inline=True, fontsize=10)
#cbar2 = fig2.colorbar(CS2)
#ax2.set_ylabel('n_qp_star [um^-3]')
#ax2.set_xlabel('tau_max [us]')
#cbar2.set_label('NEP_gr [aW/rtHz]')#, rotation=270)
plt.show()

tau_max = np.r_[600:1200:200j] * microsecond
n_qp_star = np.r_[300:800:200j] * 1/um**3
R = 1./(tau_max * n_qp_star)

tau = 1/R/n_th
f3dB = 1/2/pi/tau
NEPgr = 2*GT/(kappa*n_th) / np.sqrt(R*Vsc)
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(R/um**3, NEPgr/aW )
ax2 = ax.twinx()
ax2.plot(R/um**3, f3dB/kHz,'r')
ax.set_xlabel('R x10^-6 [um^3 us^-1]')
ax.grid(which='both')
ax.set_ylabel('NEP [aW/rtHz]')
ax2.set_ylabel('f3dB [kHz]')
ax.set_title("GR noise and rolloff frequency predictions with 4pW loading")
plt.savefig("GR_and_f3dB_4pW_prediction.png")
#plt.plot(T*1e3, P_cooling/pW, 'k--')
plt.show()

