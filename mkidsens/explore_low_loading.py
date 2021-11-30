#! /usr/bin/env python3

# My version of the resonator code. Most of the calculations are adapted from the ipython notebook that I have in this folder.
# I initialize the resonator bolometer with some of its parameters so I can easily adjust it.

import numpy as np
from math import pi
from scipy.constants import h,k,c,m_e,N_A
import matplotlib.pyplot as plt
from scipy.special import kn, iv

pW = 1e-12
MHz = 1e6
kHz = 1e3
Hz = 1
um = 1e-6
us = 1e-6
Kelvin = 1
microsecond = 1e-6
mJ = 1e-3
mol = 1
cm = 1e-2
g = 1
uOhm = 1e-6
pF = 1e-12
nH = 1e-9
aW = 1e-18

k_B = k

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

hbar = h/(2*pi)
Q_int = 1e8
T_amp = 3.8 * Kelvin
eta_read = 0.01


# Inductor Material Properties
## For Hafnium
#T_c = 0.395 * Kelvin
#gamma = 2.15 * mJ/mol/Kelvin**2
#density = 13.31 * g/cm**3
#A_r = 178.486 * g/mol
#rho = 97 * uOhm * cm

## For Aluminum
T_c = 0.5 * Kelvin
gamma = 1.35 * mJ/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g/mol
rho = 1.15 * uOhm * cm

Delta = 1.764 * k_B * T_c



# Material properties of the Aluminum superconductor
tau_max = 955.1 * microsecond
#tau_max = 500 * u.microsecond
n_qp_star = 100 * 1/um**3
#tau_max = 688.1 * u.microsecond
#n_qp_star = 663 * 1/u.um**3
#tau_max = 500.1 * u.microsecond
#n_qp_star = 100 * 1/u.um**3
#R = 2.1 *u.um**3 / u.s
R = 1./(tau_max * n_qp_star)
#n_qp_star = 518 * 1/u.um**3
#tau_max = 1./(R * n_qp_star)
#tau_max = 488.1 * u.microsecond
#R *= 0.5

# Physical properties of the superconductor + resonator capacitor
t = 0.05 * um
w_trace = 1 * um #width of inductor trace
s_trace = 1 * um #spacing between inductor meanders
N_sq = 16200  # Number of squares
l = N_sq * w_trace # Total length of the inductor (approximately. More exact is +21um + ...)
A = t * w_trace # cross-sectional area
V_sc = l * A
L_g = 6.1 * nH # From simulations
Rs = (rho/t) # surface resistance in ohms/sq
L_k = (h * Rs/(2 * np.pi**2 * Delta) * N_sq) # Kinetic inductance contribution

L = L_g + L_k # total inductance
alphak = (L_k/L)



P_opt = 1.00*pW #((Vdc/Rb)**2 * Rh).to(u.pW)

n = 2.000 # conductivity index = beta + 1
K_leg = 5 * pW/Kelvin**(n+1)

f_r = 336*MHz
f_g = 336*MHz
omega_r = 2*pi*f_r

Pg_dBm = -110
Pg = 1e-3 * 10**(Pg_dBm/10)

plt.figure(123, figsize=(10,10))
plt.figure(321, figsize=(10,10))
plt.figure(213, figsize=(10,10))

T_0 = 50e-3 * Kelvin # Previously 0.23K Temperature of the thermal bath
Tstart = 65e-3

T = np.r_[Tstart:0.2:1000j]
kappa = 0.5 + Delta/(k_B*T)
#chi_ph = 0.658 # flink factor for phonon noise
chi_ph = (n+1)/(2*n+3)*((T_0/T)**(2*n+3)-1)/((T_0/T)**(n+1)-1)
G = K_leg * (n+1)*T**n

def get_loading(T):
    return K_leg*((T*1e-3)**(n+1) - T_0**(n+1))/pW


def get_islandtemp(P):
    return 1e3*(P*pW/K_leg + T_0**(n+1))**(1./(n+1))

eta = h * f_g / (2*k_B * T)
S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k_B*T))*np.sinh(eta)*K_0(eta)
S_2 = 1 + np.sqrt(2*Delta/(pi*k_B*T)) * np.exp(-eta) * I_0(eta)
beta = S_2/S_1

N_0 = (3 * (gamma * (density/A_r)))/(2*pi**2 * k_B**2)
Gamma_gen = 0
n_th = 2*N_0 * np.sqrt(2*pi* k_B * T* Delta)*np.exp(-Delta/(k_B*T))

# Quality factors
n_qp = np.sqrt((n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star
tau_qp = tau_max/(1 + n_qp/n_qp_star)
Q_qp = (2 * N_0 * Delta)/(alphak * S_1 * n_qp)
Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * T))/np.sinh(eta)/K_0(eta)
#Q_c = 3000
Q_i = 1./(1/Q_qp + 1./Q_int)
Q_c = Q_i
Q_r  = 1./(1./Q_c + 1./Q_i)
chi_c = 4*Q_r**2/(Q_c*Q_i)
chi_g = 1
x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)
P_diss = (chi_g*chi_c/2) * Pg
Pc_tls = 1e-3*10**(-95.3/10)
Ec_tls = (Pc_tls*Q_i/omega_r)
Nc_tls = (Ec_tls/(hbar*omega_r))
E_stored = (P_diss*Q_i/omega_r)
N_ph = (E_stored/(hbar*omega_r)) # number of #microwave photons
#print (N_ph)
#print (Nc_tls)
#exit()
#responsivity
S = f_r * x * kappa / (G * T)

NEP_ph = np.sqrt(4*chi_ph*k_B*T**2*G)
NEP_amp = (G*T/(kappa*x))*(2/Q_i)*np.sqrt(k_B*T_amp/Pg)
#NEP_amp = (f_r/S)*(Q_c/(2*Q_r**2))*np.sqrt(k_B*T_amp/Pg)
NEP_gr = (2*G*T/n_qp/kappa)/np.sqrt(R*V_sc)
# TLS NEP
alpha_tls = 0.5
beta_tls = 2
kappatls0 = 3.2e-16/Hz*Kelvin**beta_tls*Hz**0.5/np.sqrt(Nc_tls)
#print (kappatls0)
#exit()
nu_tls = 1 * Hz
# TLS spectrum at 1 Hz
Stls = kappatls0/np.sqrt(1 + N_ph/Nc_tls)*T**(-beta_tls)* nu_tls**(-alpha_tls)
#print (Stls)
NEP_tls = (Stls**0.5*f_r/S)

NEP_total = np.sqrt(NEP_ph**2 + NEP_amp**2 + NEP_gr**2 + NEP_tls**2)


#fig, ax = plt.subplots(figsize=(10,10))
#ax.plot(T*1e3, NEP_total/aW, 'r',label='Total')
#ax.plot(T*1e3, NEP_ph/aW, color='green', ls='dashed',label='Phonon')
#ax.plot(T*1e3, NEP_gr/aW, color='blue', ls='dotted', label='GR')
#ax.plot(T*1e3, NEP_amp/aW, color='black', ls='-.', label='Amplifier')
#ax.plot(T*1e3, NEP_tls/aW, color='cyan', ls='dashed', label='TLS @ 1 Hz')
#ax2 = ax.secondary_xaxis('top', functions=(get_loading, get_islandtemp))
#ax2.set_xlabel('Island Loading [pW]')
#ax.set_xlim(left=Tstart*1e3)
#ax.grid()
#ax.set_xlabel('Island Temperature [mK]')
#ax.set_ylabel('NEP [aW/$\sqrt{\mathrm{Hz}}$]')
#ax.legend(loc='upper left', title='NEP')
##ax.set_title(r"Tc = %1.1f K"%(T_c))
##ax2 = ax.twinx()
##ax2.plot(T*1e3, S/kHz*pW, color='blue')
##ax2.set_ylabel('Responsivity [KHz/pW]')
##plt.savefig('responsivityNEP_vs_temperature_%1.1fK.pdf'%(T_c))
#plt.savefig('low_loading_responsivityNEP_vs_temperature_%1.1fK.png'%(T_c))

fig, ax = plt.subplots(figsize=(10,10))
ax.semilogy(T*1e3, NEP_total/aW, 'r',label='Total')
ax.semilogy(T*1e3, NEP_ph/aW, color='green', ls='dashed',label='Phonon')
ax.semilogy(T*1e3, NEP_gr/aW, color='blue', ls='dotted', label='GR')
ax.semilogy(T*1e3, NEP_amp/aW, color='black', ls='-.', label='Amplifier')
ax.semilogy(T*1e3, NEP_tls/aW, color='cyan', ls='dashed', label='TLS @ 1 Hz')
ax2 = ax.secondary_xaxis('top', functions=(get_loading, get_islandtemp))
ax2.set_xlabel('Island Loading [pW]')
ax.set_xlim(left=Tstart*1e3)
ax.annotate('$T_{\mathrm{bath}} = 250$ mK', (350, 60),
        **{'color':'black', 'fontsize':20, 'ha':'center',
            'va':'center', 'bbox': dict(boxstyle="round",
                fc='white', ec='black', pad=0.2)})
ax.grid()
ax.set_xlabel('Island Temperature [mK]')
ax.set_ylabel('NEP [aW/$\sqrt{\mathrm{Hz}}$]')
ax.legend(loc='lower right', title='NEP')
#ax2 = ax.twinx()
#ax2.semilogy(T*1e3, S/kHz*pW, color='blue')
#ax2.set_ylabel('Responsivity [KHz/pW]')
#plt.savefig('low_loading_responsivityNEP_vs_temperature_%1.1fK_log.pdf'%(T_c))
plt.savefig('low_loading_responsivityNEP_vs_temperature_%1.1fK_log.png'%(T_c))

plt.figure(123)
ax = plt.gca()
ax.plot(T*1e3, Q_r, label='%1.1fK'%(T_c))
ax2 = ax.secondary_xaxis('top', functions=(get_loading, get_islandtemp))
ax2.set_xlabel('Island Loading [pW]')
#ax.set_xlim(left=Tstart*1e3)
ax.grid()
ax.set_xlabel('Island Temperature [mK]')
ax.set_ylabel('$Q_r$')
ax.legend(title='$T_c$', loc='best')
#plt.savefig('low_loading_Qr_vs_temperature.pdf')
plt.savefig('low_loading_Qr_vs_temperature.png')

plt.figure(321)
ax = plt.gca()
ax.plot(T*1e3, NEP_total/aW, label='%1.1fK'%(T_c))
ax2 = ax.secondary_xaxis('top', functions=(get_loading, get_islandtemp))
ax2.set_xlabel('Island Loading [pW]')
ax.grid()
ax.set_xlabel('Island Temperature [mK]')
ax.set_ylabel('NEP [aW/$\sqrt{\mathrm{Hz}}$]')
ax.legend(loc='lower right', title='$T_c$')
#plt.savefig('low_loading_NEP_vs_temperature.pdf')
plt.savefig('low_loading_NEP_vs_temperature.png')

plt.figure(213)
ax = plt.gca()
ax.plot(T*1e3, S/kHz*pW, label='%1.1fK'%(T_c))
ax2 = ax.secondary_xaxis('top', functions=(get_loading, get_islandtemp))
ax2.set_xlabel('Island Loading [pW]')
ax.grid()
ax.set_xlabel('Island Temperature [mK]')
ax.set_ylabel('Responsivity [kHz/pW]')
ax.legend(loc='upper left', title='$T_c$')
#plt.savefig('low_loading_responsivity_vs_temperature.pdf')
plt.savefig('low_loading_responsivity_vs_temperature.png')

plt.show()

#Tcs = np.r_[0.8:2:1000j]
#Ts = np.r_[0.25:1:1000j]
#t, tc = np.meshgrid(Ts, Tcs)
#g = K_leg * (n+1)*t**n
#Delta = 1.764 * k_B * tc
#
#eta = h * f_g / (2*k_B * t)
#S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k_B*t))*np.sinh(eta)*K_0(eta)
#S_2 = 1 + np.sqrt(2*Delta/(pi*k_B*t)) * np.exp(-eta) * I_0(eta)
#
#n_th = 2*N_0 * np.sqrt(2*pi* k_B * t* Delta)*np.exp(-Delta/(k_B*t))
#
## Quality factors
#n_qp = np.sqrt((n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star
#tau_qp = tau_max/(1 + n_qp/n_qp_star)
#Q_qp = (2 * N_0 * Delta)/(alphak * S_1 * n_qp)
#Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * t))/np.sinh(eta)/K_0(eta)
#Q_c = 22400
#Q_i = 1./(1/Q_qp + 1./Q_int)
#Q_r  = 1./(1./Q_c + 1./Q_i)
#chi_c = 4*Q_r**2/(Q_c*Q_i)
#x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)
#
##responsivity
#S = f_r * x * kappa / (g * t)
#
#NEP_ph = np.sqrt(4*chi_ph*k_B*t**2*g)
#NEP_amp = (2/S)*np.sqrt(k_B*T_amp/Pg)/(chi_c*Q_i)
#NEP_gr = (2*g*T/n_qp/kappa)/np.sqrt(R*V_sc)
#
#NEP_total = np.sqrt(NEP_ph**2 + NEP_amp**2 + NEP_gr**2)
#
#plt.figure(figsize=(10,10))
#im = plt.imshow(S/kHz*pW, origin='lower', interpolation='bilinear')
#plt.xlabel('Island Temperature [mK]')
#plt.ylabel('Tc [K]')
#im.set_extent([Ts[0], Ts[-1], Tcs[0], Tcs[-1]])
#plt.colorbar()
#plt.savefig('responsivity_vsTandTc.png')
#
#plt.figure(figsize=(10,10))
#im = plt.imshow(NEP_gr/aW, origin='lower', interpolation='bilinear')
#plt.xlabel('Island Temperature [mK]')
#plt.ylabel('Tc [K]')
#im.set_extent([Ts[0], Ts[-1], Tcs[0], Tcs[-1]])
#plt.colorbar()
#plt.savefig('NEP_vsTandTc.png')