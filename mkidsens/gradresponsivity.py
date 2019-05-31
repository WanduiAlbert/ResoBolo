#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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

K = 141.2*pW
n = 2
Tc = 1.329
def get_islandtemperature(Tbath, P):
	return (P/K + Tbath**(n+1))**(1./(n+1))

def get_islandpower(T, Tbath):
	return K*(T**(n+1) - Tbath**(n+1))

def get_rjtemp(P, nu0, bw):
    return P/(k*nu0*bw)

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


fr0 = 337.4*MHz
spacing = 5*MHz
reso_freqs = fr0# + np.arange(3)*spacing
Tbath = 0.25
Qc = 16000
dQc = 1./(Qc)# + np.random.randn(3)*1200)
delta_nu = np.r_[-2:2:5000j]*MHz

loading = np.r_[0:15:30j]*pW
dP = (loading[1]-loading[0])

Tis = get_islandtemperature(Tbath, loading)
S = get_responsivity(Tis, fr0)
p = np.polyfit(loading/pW, S/kHz*pW, 2)
print (p)
Sfit = np.polyval(p, loading/pW)
gradp = [2*p[0], p[1]]

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(loading/pW, S/kHz*pW)
ax.plot(loading/pW, Sfit)
ax.grid()
ax.set_xlabel('P [pW]')
ax.set_ylabel('Responsivity [kHz/pW]')
plt.savefig('responsivity_vs_P.png')
#plt.show()

gradS = np.gradient(S)/dP
gradSfit = gradp[0]*(loading/pW) + gradp[1]
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(loading/pW, gradS/kHz*pW**2)
ax.plot(loading/pW, gradSfit)
ax.grid()
ax.set_xlabel('P [pW]')
ax.set_ylabel('dS/dP [$\mathrm{kHz}/\mathrm{pW}^2$]')
plt.savefig('gradresponsivity_vs_P.png')
plt.show()

#xMB, dQi = get_xQMB(Tis, fr0)
#fr = fr0*(1+xMB)
#dQr = dQi + dQc
#x = delta_nu/fr
#S21 = (1 - dQc/(dQr + 2j*x))
#freq = (delta_nu+fr)
#magS21 = 20*np.log10(np.abs(S21))
#plt.plot(freq/MHz, magS21, label="%1.1fK Loading"%(Tr))



