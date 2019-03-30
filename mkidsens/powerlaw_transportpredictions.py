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

np.set_printoptions(precision=4)

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)
hbar = h/2/pi
MHz = 1e6
kHz = 1e3
pW = 1e-12
um  = 1e-6
microsecond = 1e6
mJ = 1e-3
mol = 1
Kelvin = 1
cm = 1e-2
g = 1
Ohm = 1
# Material properties of the Aluminum superconductor
#tamax = 500 * microsecond
#n_qp_star = 100 * 1/um**3
gamma = 1.35 * mJ/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g/mol
rho = 1.15 * Ohm * cm

Tc = 1.329
Delta = 1.764*k*Tc
alphak = 0.4

Tbath = 80e-3
K = 140*pW #/K^3
Nsq = 16200
Asq = 1*um**2
t = 0.05*um
Sig_kap = 850
Sig_eph = 0.2e9
beta = 3
T = np.r_[90:800:100j]*1e-3
N_0 = ((3 * (gamma * (density/A_r)))/(2*pi**2 * k**2))
n_th = 2 * N_0 * np.sqrt(2*pi* k * T* Delta)*np.exp(-Delta/(k*T))

fr0 = 337.4*MHz
fg = fr0 + 10*kHz
eta = h * fg / (2*k*T)
S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k*T))*np.sinh(eta)*K_0(eta)
S_2 = 1 + np.sqrt(2*Delta/(pi*k*T))*np.exp(-eta)*I_0(eta)
x = -(alphak * S_2 * n_th)/(4 * N_0 * Delta)
Qi = (2 * N_0 * Delta)/(alphak * S_1 * n_th)

fr = fr0*(1 + x)



Qc = 16000
phic = 0.0
Qr = 1./(1./Qc + 1./Qi)
Qe = Qc*np.cos(phic)*np.exp(1j*phic)

Z0 = 50

#f0 = 337.4*MHz
w0 = 2*pi*fr
dws = 2*pi*(np.arange(0, 800, 50) - 400)*kHz
P_cooling = K*(T**beta - Tbath**beta)

dw = 0
#w = np.r_[-10:10:100j]*kHz

#Zl = 0.5*Z0*(Qe/Qr)*(1 - 2j*Qr*dw/w0 - Qr/Qe)
#Rl = np.real(Zl)
#I_load = Vsource * (Zl/(Z0/Zl + 2))

Vsource = 6.3e-4#14e-6
Preadout = Vsource**2/4/Z0
print (Preadout/pW)

#Pdiss = np.abs(I_load)**2 * Rl
plt.figure(figsize=(10,10))
for i, dw in enumerate(dws):
    #print (dw/2/pi/kHz)
    #if i > 0: continue
    dx = dw/w0
    Pdiss = Preadout*(2*Qr**2/Qi/Qc)/(1 + 4*Qr**2*dx**2)

    plt.semilogy(T*1e3, Pdiss/pW, label="%1.1f"%(dw/kHz))
#plt.plot(T*1e3, P_cooling/pW, 'k--')
plt.show()


