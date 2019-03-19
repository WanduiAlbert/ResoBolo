#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.constants import h, k
from scipy.special import kn, iv
import pdb
K0 = lambda x: kn(0, x)
I0 = lambda x: iv(0, x)

MHz = 1e6
GHz = 1e9
kHz = 1e3
pW = 1e-12
eV = 1.6e-19 # Joules
um = 1e-6
N0 = 1.72e10/um**3/eV
alpha_k = 0.46
mK = 1e-3

K = 151.2*pW
n = 2
Tc = 1.329
def get_islandtemperature(Tbath, P):
	return (P/K + Tbath**(n+1))**(1./(n+1))

def get_xQMB(T, nu):

	delta = 3.5*k*Tc/2
	q = (h*nu/2/k/T)
	nqp = 2*N0*np.sqrt(2*pi*delta*k*T)*np.exp(-delta/(k*T))
	S1 = 2/pi * np.sqrt(2*delta/(pi*k*T))*np.sinh(q)*K0(q)
	S2 = 1 + np.sqrt(2*delta/(pi*k*T))*np.exp(-q)*I0(q)

	x = -alpha_k * S2 * nqp / (4*N0*delta)
	Qi = (2*N0*delta)/(alpha_k*S1*nqp)
	return x, 1./Qi

def get_responsivity(T, nu):
	G = K*(n+1)*T**n
	delta = 3.5*k*Tc/2
	kappa = 0.5 + delta/(k*T)
	q = (h*nu/2/k/T)
	nqp = 2*N0*np.sqrt(2*pi*delta*k*T)*np.exp(-delta/(k*T))
	S2 = 1 + np.sqrt(2*delta/(pi*k*T))*np.exp(-q)*I0(q)
	x = alpha_k * S2 * nqp / (4*N0*delta)
	S = (nu * x * kappa)/(G * T)
	return S

def get_loading(T, nu0, bw):
	return k*T*nu0*bw


load_temps = np.arange(300.)#[2, 18, 77, 300]
nu0 = 150*GHz
bw = 0.15

fr0 = 337.4*MHz
spacing = 5*MHz
reso_freqs = fr0# + np.arange(3)*spacing
Tbath = 0.08
Qc = 16000
dQc = 1./(Qc)# + np.random.randn(3)*1200)
delta_nu = np.r_[-2:2:5000j]*MHz

plt.figure(1)
for Tr in load_temps[::15]:
	power = get_loading(Tr, nu0, bw)
	print ("The loading power is %1.3fpW"%(power/pW))
	Tis = get_islandtemperature(Tbath, power)
	print ("The island temperature is %1.3fmK"%(Tis/mK))
	xMB, dQi = get_xQMB(Tis, fr0)
	fr = fr0*(1+xMB)
	dQr = dQi + dQc
	x = delta_nu/fr
	S21 = (1 - dQc/(dQr + 2j*x))
	freq = (delta_nu+fr)
	magS21 = 20*np.log10(np.abs(S21))
	plt.plot(delta_nu/MHz, magS21, label="%1.1fK Loading"%(Tr))


plt.xlabel("Frequency [MHz]")
plt.ylabel("|S21|")
plt.title("Resonator Loading vs Load Temperature")
plt.grid()
plt.legend(loc="best")
plt.show()

reso_bw = np.zeros_like(load_temps)
island_temps = np.zeros_like(load_temps)
resps = np.zeros_like(load_temps)
crosstalk_bw = np.zeros_like(load_temps)
heater_powers = np.zeros_like(load_temps)
ct_thresh = 0.01
for i, Tr in enumerate(load_temps):
	#pdb.set_trace()
	power = get_loading(Tr, nu0, bw)
	heater_powers[i] = power
	Tis = get_islandtemperature(Tbath, power)
	xMB, dQi = get_xQMB(Tis, fr0)
	resps[i] = get_responsivity(Tis, fr0)
	fr = fr0*(1+xMB)
	dQr = dQi + dQc
	reso_bw[i] = (0.5*fr*dQr)
	crosstalk_bw[i] = (0.5*fr*dQr)*np.sqrt((1-ct_thresh)/ct_thresh)
	island_temps[i] = Tis
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, reso_bw/kHz)
ax.set_xlabel("Load Temperature [K]")
ax.set_ylabel("Resonator Bandwidth [kHz]")
ax2 = ax.twiny()
ax2.plot(island_temps/mK, reso_bw/kHz, linestyle='None')
ax2.set_xlabel("Island Temperature [mK]")
#ax2.set_xticks(island_temps[::30]/mK)
ax.grid()
plt.savefig("resonatorbandwidth_vs_loading.png")
plt.show()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, resps*pW/kHz)
ax.set_xlabel("Load Temperature [K]")
ax.set_ylabel("Responsivity [kHz/pW]")
ax2 = ax.twiny()
ax2.plot(island_temps/mK, resps*pW/kHz, linestyle='None')
ax2.set_xlabel("Island Temperature [mK]")
#ax2.set_xticks(island_temps[::30]/mK)
ax.grid()
plt.savefig("resonator_responsivity_vs_loading.png")
plt.show()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, crosstalk_bw/MHz)
ax.set_xlabel("Load Temperature [K]")
ax.set_ylabel("Resonator Spacing [MHz]")
ax.set_title("Resonator Spacing for crosstalk < 1% level")
ax2 = ax.twiny()
ax2.plot(island_temps/mK, crosstalk_bw/MHz, linestyle='None')
ax2.set_xlabel("Island Temperature [mK]")
#ax2.set_xticks(island_temps[::30]/mK)
ax.grid()
plt.savefig("crosstalkbandwidth_vs_loading.png")
plt.show()

x = heater_powers/K
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, island_temps/mK)
ax.plot(load_temps[:20], (Tbath + 1./(n+1)*x[:20]*Tbath**-2)/mK, 'r--')
ax.plot(load_temps, (x**(1/(n+1)) + 1./(n+1)*x**(-n/(n+1))*Tbath**(n+1))/mK, 'k--')
ax.set_xlabel("Load Temperature [K]")
ax.set_ylabel("Island Temperature [mK]")
ax.set_title("Island vs Load Temperature with 0.25K bath temperature")
ax.grid()
plt.savefig("islandtemp_vs_loading.png")
plt.show()
