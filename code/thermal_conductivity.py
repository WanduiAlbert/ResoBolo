#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from scipy import optimize
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.special import expi
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter, LogFormatter, LogLocator, StrMethodFormatter
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
Tc = 1.284
mK = 1e-3

cm = 1e-2
um = 1e-6
nm = 1e-9
Kelvin = 1
pW = 1e-12
MHz = 1e6

# Bolometer dimensions
Tc_leg = 1.2
L0 = 2.45e-8
rho = 1.255e-8 # residual resistivity of Al
rho = 5.432e-8 # residual resistivity of Nb
Nlegs = 6
l = 300 *um
A_metal = 120 * nm * (6*um*(Nlegs - 2) + 10*um*2)
A_SiN = 300 * nm * (6*um*(Nlegs - 2) + 10*um*2)
G0 = 55*pW#/K
T0 = 0.38
n = 3
#K_SiN = (l/A_SiN)*G0*T0**(1-n)
K_SiN = 141*pW#/K
# I'll actually make my own choice of K0 and b here
b = 1.8 * Tc_leg
K_metal = b**2*(L0/rho)*6/pi**2
#K_metal = b**2*L0*6/pi**2
#K *= 1e-7
print (K_metal*A_metal/l)
print (K_SiN)

def thermal_model(T, A, b):
    return A*(b/T)*np.exp(-b/T)

def pure_lattice(T, K, n):
    return K*T**n

def thermal_model_lattice(T, A, b, K, n):
    k1 = (A_metal/l)*(b/T)*np.exp(-b/T)
    k2 = K_SiN * T**n
    return k1 + k2
    #return 1./(1./k1 + 1./k2)

def get_SC_conductivity(T):
    return -(A_metal/l)*K_metal*expi(-b/T)

def get_SiN_conductivity(T):
    return K_SiN*T**n

def get_SiN_G(T):
    return K_SiN*n*T**(n-1)

def get_SC_G(T):
    return (A_metal/l)*K_metal*np.exp(-b/T)/T

#fn = 'Thermal_Conductivity_of_Tin_0.1_to_1K.csv'
#fn = 'Thermal_Conductivity_of_Zinc_0.1_to_1K.csv'
#fn = 'Thermal_Conductivity_of_Al_0.1_to_1K.csv'
#fn = 'Thermal_Conductivity_of_Niobium_0.1_to_1K.csv'
#data = np.loadtxt(fn, delimiter=',', skiprows=1).T
#T = np.r_[0.1:1:1000j]
#data[data==0.] = np.nan
#labels = ['S.J. Laredo (1995) pure','G.M. Graham (1958) pure',
#        'S.J. Laredo (1995) 0.35% In', 'N.V. Zavaritskii (1958a) 0.002% impurity']
#labels = ['N.V. Zavaritskii (1961a) perpendicular',
#        'N.V. Zavaritskii (1961a) parallel',
#        'N.V. Zavaritskii (1958b) 0.0001% impurities',
#        'N.V. Zavaritskii (1958b) normal with 6mT longitudinal B']
#labels = ['N.V. Zavaritskii (1958b)',
#        'L. Passell (1958)',
#        'C.B. Satterthwaite (1962) RRR = 430']
#labels = ['Connolly & Mendelssohn (1962) RRR=60.5',
#        'Carlson & Satterthwaite (1970)',
#        'Connolly & Mendelssohn (1962) RRR=120']
#K0s = []
#bs = []
#plt.figure(1, figsize=(10,10))
#plt.figure(2, figsize=(10,10))
#for i in range(0, 6, 2):
#    mask = np.logical_not(np.isnan(data[i]))
#    plt.figure(1)
#    plt.scatter(data[i][mask], data[i + 1][mask], label=labels[i//2])
#    # power law fit to the conductivity data
#    #p = np.polyfit(np.log(data[i][mask]), np.log(data[i+1][mask]), 1)
#    #fit = np.exp(np.polyval(p, np.log(T)))
#    # Exponential fit to the conductivity data
#    p = np.polyfit(1./data[i][mask], np.log(data[i+1][mask]), 1)
#    #fit = np.exp(np.polyval(p, 1./T))
#    #print (p)
#    #p0 = [K , b]#, K_SiN]
#    p0 = [0.1 , 3]
#    popt, pcov = optimize.curve_fit(pure_lattice, data[i][mask],
#            data[i+1][mask], p0=p0, method='lm')
#    print (popt)
#    #fit = thermal_model_lattice(T, *popt)
#    #fit = thermal_model(T, *popt)
#    fit = pure_lattice(T, *popt)
#    K0s.append(popt[0]*popt[1]/cm/Kelvin)
#    bs.append(popt[1])
#    #K0s.append(np.exp(p[1])/cm/Kelvin)
#    #bs.append(-p[0])
#    plt.plot(T, fit, color=colors[i//2], ls='--')
#    plt.figure(2)
#    plt.scatter(data[i][mask], data[i + 1][mask], label=labels[i//2])
#    plt.plot(T, fit, color=colors[i//2], ls='--')
##plt.scatter(data[2], data[3], label=)
##plt.scatter(data[4], data[5], label=)
##plt.scatter(data[6], data[7], label=)
#plt.figure(1)
#plt.plot(T, K*np.exp(-b/T)/T, 'r-', label='Model')
#plt.grid()
#plt.yscale('log')
#plt.ylim(bottom=1e-8)
#plt.legend(loc='best')
#plt.xlabel('Temperature [K]')
#plt.ylabel('Thermal Conductivity of Nb [W/m/K]')
#plt.savefig('niobium_conductivity.png')
##plt.ylabel('Thermal Conductivity of Al [W/m/K]')
##plt.savefig('aluminum_conductivity.png')
##plt.ylabel('Thermal Conductivity of Zn [W/m/K]')
##plt.savefig('zinc_conductivity.png')
##plt.ylabel('Thermal Conductivity of Sn [W/m/K]')
##plt.savefig('tin_conductivity.png')
#plt.close()
#plt.figure(2)
#plt.plot(1./T, K*np.exp(-b/T)/T, 'r-', label='model')
#plt.grid()
#plt.ylim(bottom=1e-8)
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(loc='best')
#plt.xlabel('T [K]')
##plt.ylabel('Thermal Conductivity of Al [W/m/K]')
##plt.savefig('aluminum_conductivity_vs_1overT.png')
#plt.ylabel('Thermal Conductivity of Nb [W/m/K]')
#plt.savefig('niobium_conductivity_vs_1overT.png')
#plt.close()
#print (K0s)
#print (np.array(bs)/Tc_leg)
#exit()
#b = np.mean(bs)
#K = K0s[-1]* A/l#*1e-1


print (b/Tc, K_SiN/pW, K_metal/pW)



def get_conductance(Tbath, T):
    #return K/T*np.exp(-b/T)
    T_fine = np.r_[Tbath:1:10000j]
    P_fine = (get_SC_conductivity(T_fine) - get_SC_conductivity(Tbath))
    #P_fine = 0
    P_fine += (get_SiN_conductivity(T_fine) - get_SiN_conductivity(Tbath))
    spl = UnivariateSpline(T_fine, P_fine, k=3, s=0)
    spl_deriv = spl.derivative()
    return spl_deriv(T)

def get_islandtemperature(Tbath, P):
    T_fine = np.r_[Tbath:1:10000j]
    P_fine = (get_SC_conductivity(T_fine) - get_SC_conductivity(Tbath))
    #P_fine = 0
    P_fine += (get_SiN_conductivity(T_fine) - get_SiN_conductivity(Tbath))
    #print (K/pW, b)
    #plt.figure(figsize=(10,10))
    #plt.plot(T_fine, P_fine/pW)
    #plt.xlabel('Temperature [K]')
    #plt.ylabel('Power [pW]')
    #plt.grid()
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.savefig('PvsT.png')
    #exit()
    spl = UnivariateSpline(P_fine, T_fine, k=3, s=0)
    #interpolator = interp1d(P_fine, T_fine, kind='cubic')
    #return interpolator(P)
    return spl(P)
    #lhs = P/K + Tbath*np.exp(-b/Tbath) + b*expi(-b/Tbath)
    #lhs /= b
    #return b/get_root(lhs)


def get_islandpower(T, Tbath):
    P = (get_SC_conductivity(T) - get_SC_conductivity(Tbath))
    #P = 0
    P += (get_SiN_conductivity(T) - get_SiN_conductivity(Tbath))
    return P
    #return K*(T**(n+1) - Tbath**(n+1))

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
    #G = K*(n+1)*T**n
    G = get_conductance(Tbath, T)
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
Tbath = 0.35
Qc = 16000
dQc = 1./(Qc)# + np.random.randn(3)*1200)
delta_nu = np.r_[-2:2:5000j]*MHz

plt.figure(1)
for Tr in load_temps[::30]:
    power = get_loading(Tr, nu0, bw)
    print ("The loading power is %1.3fpW"%(power/pW))
    Tis = get_islandtemperature(Tbath, power)
    print ("The island temperature is %1.3fmK"%(Tis/mK))
    #exit()
    xMB, dQi = get_xQMB(Tis, fr0)
    fr = fr0*(1+xMB)
    dQr = dQi + dQc
    x = delta_nu/fr
    S21 = (1 - dQc/(dQr + 2j*x))
    freq = (delta_nu+fr)
    magS21 = 20*np.log10(np.abs(S21))
    plt.plot(freq/MHz, magS21, label="%1.1fK Loading"%(Tr))

#exit()

plt.xlabel("Frequency [MHz]")
plt.ylabel("|S21| [dB]")
plt.title("Resonator Loading vs Load Temperature")
plt.grid()
plt.legend(loc="best")
plt.savefig("337.4MHz_resonator_loaded.png")
plt.close()

reso_bw = np.zeros_like(load_temps)
island_temps = np.zeros_like(load_temps)
resps = np.zeros_like(load_temps)
crosstalk_bw = np.zeros_like(load_temps)
heater_powers = np.zeros_like(load_temps)
conductances = np.zeros_like(load_temps)
ct_thresh = 0.01
for i, Tr in enumerate(load_temps):
    #pdb.set_trace()
    power = get_loading(Tr, nu0, bw)
    heater_powers[i] = power
    Tis = get_islandtemperature(Tbath, power)
    G = get_conductance(Tbath, Tis)
    xMB, dQi = get_xQMB(Tis, fr0)
    resps[i] = get_responsivity(Tis, fr0)
    fr = fr0*(1+xMB)
    dQr = dQi + dQc
    reso_bw[i] = (0.5*fr*dQr)
    crosstalk_bw[i] = (0.5*fr*dQr)*np.sqrt((1-ct_thresh)/ct_thresh)
    island_temps[i] = Tis
    conductances[i] = G
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, reso_bw/kHz)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Bandwidth [kHz]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
#islandticks = np.r_[3:8:6j]*100
islandticks = np.array([3,6.0,6.5,6.75])*100
tickpos = get_rjtemp(get_islandpower(islandticks*mK, Tbath), nu0, bw)
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
print (tickpos)
ax2.set_xticklabels(islandticks.astype(int))
ax2.set_xlim(ax.get_xlim())
ax.grid()
plt.savefig("resonatorbandwidth_vs_loading.png")
plt.close()
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, resps*pW/kHz)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Responsivity [kHz/pW]")
ax2 = ax.twiny()
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax.set_xlim((0, 300))
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks.astype(int))
ax2.set_xlim(ax.get_xlim())
ax.grid()
plt.savefig("resonator_responsivity_vs_loading.png")
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, crosstalk_bw/MHz)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Resonator Spacing\nfor < 1% crosstalk [MHz]")
ax.set_xlim((0, 300))
ax2 = ax.twiny()
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks.astype(int))
ax2.set_xlim(ax.get_xlim())
ax.grid()
plt.savefig("crosstalkbandwidth_vs_loading.png")
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(load_temps, island_temps/mK)
ax.set_xlabel("Trj [K]")
ax.set_ylabel("Island Temperature [mK]")
ax.set_xlim((0, 300))
powerticks = np.r_[0:150:6j]
#islandticks = np.array([3,5,6,7,8,9,10])*100
ptickpos = get_rjtemp(powerticks*pW, nu0, bw)
ax2 = ax.twiny()
ax2.set_xlabel("Heater Power [pW]")
ax2.set_xticks(ptickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(powerticks.astype(int))
ax2.set_xlim(ax.get_xlim())
ax.grid()
plt.savefig("islandtemp_vs_loading.png")
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(island_temps/mK, conductances/pW)
ax.plot(island_temps/mK, get_SiN_G(island_temps)/pW, ls='--',color='g')
ax.plot(island_temps/mK, get_SC_G(island_temps)/pW, ls='--',color='r')
ax.set_xlabel("Island Temperature [mK]")
ax.set_ylabel("G [pW/K]")
#ax.set_xticklabels([250,300,350,400,450,500,550,600])
ax.set_yscale('log')
ax.set_ylim(bottom=0.1)
#ax.set_xscale('log')
#ax.set_xlim((250, 600))
#loc = LogLocator(subs='all', numticks=10)
#ax.xaxis.set_major_locator(loc)
#ax.xaxis.set_minor_locator(loc)
#formatter = ScalarFormatter()
#formatter.set_scientific(False)
#ax.yaxis.set_major_formatter(formatter)
#ax.xaxis.set_major_formatter(formatter)
#ax.xaxis.set_minor_formatter(formatter)
powerticks = np.array([0, 10,20,30,50,150])
itickpos = get_islandtemperature(Tbath, powerticks*pW)/mK
print (powerticks)
print (itickpos)
ax2 = ax.twiny()
ax2.set_xlabel("Heater Power [pW]")
ax2.set_xticks(itickpos)
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(powerticks)
ax2.set_xlim(ax.get_xlim())
ax.grid(which='both')
plt.savefig("islandtemp_vs_G.png")
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(heater_powers/pW, conductances/pW)
#ax.set_xlabel("Island Temperature [mK]")
ax.set_xlabel("Heater Power [pW]")
ax.set_ylabel("G [pW/K]")
ax.set_yscale('log')
#ax.set_xscale('log')
ax.grid()
islandticks = np.array([300, 400, 500, 600, 700])
islandticks = np.array([3,6.0,6.5,6.75])*100
tickpos = get_islandpower(islandticks*mK, Tbath)/pW
ax2 = ax.twiny()
ax2.set_xlabel("Island Temperature [mK]")
ax2.set_xticks(tickpos)
print (tickpos)
#ax2.set_xscale('log')
ax2.xaxis.set_major_formatter(FormatStrFormatter("%1.0f"))
ax2.set_xticklabels(islandticks.astype(int))
ax2.set_xlim(ax.get_xlim())
plt.savefig("heaterpower_vs_G.png")
plt.close()
