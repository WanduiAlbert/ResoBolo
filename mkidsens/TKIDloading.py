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


load_temps = np.arange(301.)#[2, 18, 77, 300]
nu0 = 150*GHz
bw = 0.25

fr0 = 337.4*MHz
spacing = 5*MHz
reso_freqs = fr0# + np.arange(3)*spacing
Tbath = 0.25
Qc = 16000
dQc = 1./(Qc)# + np.random.randn(3)*1200)
delta_nu = np.r_[-2:2:5000j]*MHz

plt.figure(1)
for Tr in load_temps[::30]:
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
    plt.plot(freq/MHz, magS21, label="%1.1fK Loading"%(Tr))


plt.xlabel("Frequency [MHz]")
plt.ylabel("|S21| [dB]")
plt.title("Resonator Loading vs Load Temperature")
plt.grid()
plt.legend(loc="best")
plt.savefig("337.4resonator_loaded.png")
plt.close()

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
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Bandwidth [kHz]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
islandticks = np.array([3,5,6,7,8,9,10])*100
tickpos = get_rjtemp(get_islandpower(islandticks*mK, Tbath), nu0, bw)
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks)
ax.grid()
plt.savefig("resonatorbandwidth_vs_loading.png")
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, resps*pW/kHz)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Responsivity [kHz/pW]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
islandticks = np.array([3,5,6,7,8,9,10])*100
tickpos = get_rjtemp(get_islandpower(islandticks*mK, Tbath), nu0, bw)
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks)
ax.grid()
plt.savefig("resonator_responsivity_vs_loading.png")
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, crosstalk_bw/MHz)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Spacing\nfor < 1% crosstalk [MHz]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
islandticks = np.array([3,5,6,7,8,9,10])*100
tickpos = get_rjtemp(get_islandpower(islandticks*mK, Tbath), nu0, bw)
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks)
ax.grid()
plt.savefig("crosstalkbandwidth_vs_loading.png")
plt.close()

x = heater_powers/K + 1e-10
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, island_temps/mK)
#ax.plot(load_temps[:20], (Tbath + 1./(n+1)*x[:20]*Tbath**-2)/mK, 'r--')
#ax.plot(load_temps, (x**(1/(n+1)) + 1./(n+1)*x**(-n/(n+1))*Tbath**(n+1))/mK, 'k--')
ax.set_xlabel("Rayleigh Jeans Temperature [K]")
ax.set_ylabel("Resonator Bandwidth [kHz]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
islandticks = np.array([3,5,6,7,8,9,10])*100
tickpos = get_rjtemp(get_islandpower(islandticks*mK, Tbath), nu0, bw)
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks)
ax.grid()
plt.savefig("islandtemp_vs_loading.png")
plt.close()


# Now lets look at the attenuation length needed for the neutral density filter
C = 20*np.log10(np.exp(1))*2*pi/300
att = np.logspace(0, 5, 100)
attdB = 10*np.log10(att)
n_imag = 0.038 # For CR-110 eccosorb
thickness = attdB/C/n_imag*(GHz/nu0)

plt.figure(figsize=(10,10))
plt.plot(attdB, thickness)
plt.grid()
plt.xlim(0, 50)
plt.ylim(0, 50)
plt.xlabel("Attenuation [dB]")
plt.ylabel("Thickness [mm]")
plt.title("CR-110 eccosorb at 150 GHz")
plt.savefig("att_vs_eccosorbthickness.png")
