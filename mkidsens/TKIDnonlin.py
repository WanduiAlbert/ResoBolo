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
ppm = 1e-6

K = 122*pW
n = 1.862
Tc = 1.284
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
Tbath = 0.25
Qc = 16000

P = 5.5 + np.r_[-0.5:0.5:100j]
P *= pW
PdBm = -90
Pg = 1e-3*10**(PdBm/10)
Ptotal = P + 0.5*Pg

T = get_islandtemperature(Tbath, Ptotal)
x, dQi = get_xQMB(T, fr0)
fr = fr0*(1 + x)

P0 = np.mean(P)
#df = fr - fr0
p = np.polyfit((P-P0)/pW, x/ppm, 2)
print (p)
xfit = np.polyval(p, (P-P0)/pW)*ppm
offset = p[2]*ppm
#x -= offset

linear_x = ((P-P0)/pW)*p[1]*ppm + offset
quadratic_x = ((P-P0)/pW)**2*p[0]*ppm
responsivity = get_responsivity(T, fr0)

fig, ax = plt.subplots(figsize=(10,10))
#ax.plot(P/pW, x/ppm, label='Detector Response')
#ax.plot(P/pW, (x - xfit)/ppm)
ax.plot(P/pW, linear_x/ppm, label='Linear Response')
#ax.plot(P/pW, (x - linear_x)/ppm)
ax.grid()
ax.set_xlabel('Power [pW]')
ax.set_ylabel('x [ppm]')
#ax.legend(loc='upper right')
plt.savefig('linear_response.png')
#plt.show()


fig, ax = plt.subplots(figsize=(10,10))
#ax.plot(P/pW, x/ppm)
ax.plot(P/pW, quadratic_x/ppm)
#ax.plot(P/pW, (x - xfit)/ppm)
ax.grid()
ax.set_xlabel('Power [pW]')
ax.set_ylabel('Deviation in x [ppm]')
plt.savefig('deviation_from_linear_response.png')
plt.show()
