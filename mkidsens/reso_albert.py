#! /usr/bin/env python3

# My version of the resonator code. Most of the calculations are adapted from the ipython notebook that I have in this folder.
# I initialize the resonator bolometer with some of its parameters so I can easily adjust it.

import numpy as np
import astropy.units as u
from astropy.constants import h,k_B, c, m_e
import matplotlib.pyplot as plt
from scipy.special import kn, iv

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

hbar = h/(2*np.pi)
dBm = u.dB(u.mW)

# Resonator Parameters
gamma_s = 1
Q_int = 187519
# f_g = 205.128 * u.MHz
f_g = 328 * u.MHz
T_amp = 3 * u.Kelvin
eta_read = 0.1

# I'll keep these here for now but I'm more interested in doing calculations for a dark run.
nu_opt = 150 * u.GHz
dnu_opt = 0.25
eta_opt = 0.25
T_rj = 31* u.Kelvin #0.1 * u.Kelvin
N_pol = 1

# We are using the heater pad to mimic the optical load on the bolometer
Rh = 0.106 * u.Ohm
Rb = 854 * u.kOhm
Vdc = 7.0 * u.Volt

P_opt = ((Vdc/Rb)**2 * Rh).to(u.pW)
print ("The optical power is {0:2.2f}".format(P_opt))

gamma_leg = 1.6#2.65
K_leg = 120 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.32 * u.Kelvin
T_0 = 0.23 * u.Kelvin # Temperature of the thermal bath
C_b = 1 * u.picoJoule/u.Kelvin

print ("Resonator parameters set.")

# Material properties of the Aluminum superconductor
tau_max = 500 * u.microsecond
n_qp_star = 100 * 1/u.um**3
gamma = 1.35 * u.mJ/u.mol/u.Kelvin**2
density = 2.7 * u.g/u.cm**3
A_r = 26.98 * u.g/u.mol
rho = 1.15 * u.uOhm * u.cm

#T_b = 0.38 * u.Kelvin # Temperature of the bolometer
# Calculated parameters

Delta = (1.764 * k_B * T_c).to('J')
# P_opt = (eta_opt * N_pol * dnu_opt * nu_opt * k_B * T_rj).to('pW')
# print (P_opt)
# P_read = 3.0 * u.pW
P_read = (-100 * dBm).to(u.pW)

print ("The readout power is {0:2.2f}".format(P_read))
x = P_read/P_opt

# Determine T_b by balancing the input and the output power to the resobolo
# T_b= ((((1 + x)* P_opt)/K_leg + T_0**gamma_leg)**(1./gamma_leg)).to('K')
# print(T_b)
T_b = 60 * u.mK

# Physical properties of the superconductor + resonator capacitor
t = 0.05 * u.um
w_trace = 1 * u.um #width of inductor trace
s_trace = 1 * u.um #spacing between inductor meanders
N_sq = 16200  # Number of squares
l = N_sq * w_trace # Total length of the inductor (approximately. More exact is +21um + ...)
A = t * w_trace # cross-sectional area
V_sc = l * A
L_g = 6.1 * u.nH # From simulations

Rs = (rho/t).to('Ohm') # surface resistance in ohms/sq
print ("The surface resistance", Rs)
L_k = (h * Rs/(2 * np.pi**2 * Delta) * N_sq).to('nH') # Kinetic inductance contribution
print ("The kinetic inductance", L_k)
L = L_g + L_k # total inductance

print ("The total inductance", L)
alpha = (L_k/L).to(1)

f_r = 329 * u.MHz
omega_r = (2*np.pi*f_r).to(1/u.s)
C = (1/(L * omega_r**2)).to(u.pF)

# C_i = 27.98 * u.pF
C_c1 = 0.19 * u.pF
C_c = C_c1 # The same now but will change once I add capacitance to ground

# C = C_c1 + C_i
C_i = C - C_c
Z0 = 50 * u.Ohm # Characteristic impedance of the line
# omega_r = (1./np.sqrt(L*C)).to('1/s')
# f_r = (omega_r/(2*np.pi)).to('MHz')

eta = (h * f_g / (2 * k_B * T_b)).to(1).value # Weird conversion because astropy
S_1 = ((2/np.pi)*np.sqrt(2*Delta/(np.pi*k_B*T_b))*np.sinh(eta * u.rad)*K_0(eta)).to(1)
S_2 = (1 + np.sqrt(2*Delta/(np.pi*k_B*T_b)) * np.exp(-eta) * I_0(eta)).to(1)
print (S_1, S_2)
beta = S_1/S_2

N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')
Gamma_gen = (eta_read * P_read/Delta).to('1/s')
n_th = (2*N_0 * np.sqrt(2*np.pi* k_B * T_b* Delta)*np.exp(-Delta/(k_B*T_b))).to('1/um^3')
P_crit = (np.pi * N_0 * Delta**3 * A**2 * Z0/(hbar * rho)).to(dBm)
print ("Critical readout power: {0:2.2f} ".format(P_crit))
print ("Quasiparticle parameters calculated.\n")

# Quality factors
Gamma_th = ((n_th * V_sc/tau_max)* (1 + 0.5 * n_th/n_qp_star)).to('1/s')
n_qp = (np.sqrt((n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star).to('1/um^3')
tau_th = (tau_max/(1 + n_th/n_qp_star)).to('us')
tau_qp = (tau_max/(1 + n_qp/n_qp_star)).to('us')
Q_qp = ((2 * N_0 * Delta)/(alpha * gamma_s * S_1 * n_qp)).to(1)
Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * T_b))/np.sinh(eta)/K_0(eta)
Q_c = (2 * C_i/(omega_r * C_c**2 * Z0)).to(1)
Q_i = 1./(1/Q_qp + 1./Q_int)
Q_r  = 1./(1./Q_c + 1./Q_i)

print ("Q factor from qp losses", Q_qp)
print ("Resonant Frequency", f_r)
print ("Internal Q factor", Q_i)
print ("Coupling Q factor", Q_c)
print ("Resonator Q factor", Q_r)
print ("Kinetic Inductance Fraction", alpha)
print ("surface resistance", Rs)

dx = ((f_g - f_r)/f_r).to(1)
S_21 = 1 - (Q_r/Q_c)* 1./(1 + 2 * 1j * Q_r * dx)
df_r = (f_r/(2*Q_r)).to('kHz')
df_g = (f_g - f_r).to('kHz')

chi_c = (4 * Q_r**2)/(Q_i * Q_c)
chi_g = 1./(1 + (df_g/df_r)**2)
chi_qp = Q_i/Q_qp

print("")
print ("Resonator Bandwidth", df_r)
print ("Coupling efficiency", chi_c)
print ("Detuning efficiency", chi_g)
print ("Fraction of Q_i from qp losses", chi_qp)

# Now we include the NEP estimates for the resobolo currently
kappa = (1/2 + Delta/k_B/T_b).to(1)
P_leg = P_opt * (1 + x) # Total power into the resobolo thermal link
G_b = (gamma_leg * P_leg/T_b).to('pW/K') # Conductance of the resobolo
P_b = G_b * T_b
tau_b = (C_b/G_b).to('s')

print("")
print("Conductance of Island ", G_b)
print ("QP time constant", tau_qp)
print ("Thermal recombination time constant", tau_th)
print ("Bolometer Time constant", tau_b.to('ms'))

s = ((chi_c* chi_qp/4) * beta * (tau_qp/tau_th) * (kappa/P_b)).to('1/pW') # ignoring the
#roll off factor due to the bolometer responsivity
sx = s/8/Q_c # frequency responsivity

print ("")
print ("resobolo responsivity below 30Hz optical fluctuations", sx)

# Optical NEP
n_opt = 1/(np.exp(h*nu_opt/k_B/T_b) - 1).to(1)
NEP_opt = ((2 * h*nu_opt * P_opt * (1 + n_opt))**0.5).to('aW/Hz^(0.5)')

# Phonon NEP
NEP_ph  = ((4 * k_B * T_b**2 * G_b)**0.5).to('aW/Hz^(0.5)')

# TLS NEP
NEP_TLS = ((P_b * tau_th * 2*Q_i/(kappa * tau_qp * chi_qp * beta)) *(2 *
    1e-19 /u.Hz)**0.5).to('aW/Hz^(0.5)')

# Amplifier NEP
P_g = (2/chi_c/chi_g)* P_read
NEP_amp = ((k_B * T_amp/P_g)**0.5 * P_b * tau_th * 4/(kappa * tau_qp * chi_c *
  chi_qp * beta)).to('aW/Hz^(0.5)')

# Shot NEP
NEP_shot = ((2 * n_qp * V_sc * (1/tau_qp + 1/tau_max))**0.5 * P_b *
    tau_th/kappa/n_qp/V_sc).to('aW/Hz^(0.5)')

NEP_total = NEP_opt + NEP_ph + NEP_TLS + NEP_amp + NEP_shot

print ("")
print ("Optical NEP", NEP_opt)
print ("Phonon NEP", NEP_ph)
print ("TLS NEP", NEP_TLS)
print ("Amplifier NEP", NEP_amp)
print ("Shot NEP", NEP_shot)
print ("Total NEP", NEP_total)
