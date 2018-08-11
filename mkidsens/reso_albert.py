#! /usr/bin/env python3

# My version of the resonator code. Most of the calculations are adapted from the ipython notebook that I have in this folder.
# I initialize the resonator bolometer with some of its parameters so I can easily adjust it.

import numpy as np
import astropy.units as u
from astropy.constants import h,k_B, c, m_e, N_A
import matplotlib.pyplot as plt
from scipy.special import kn, iv

np.set_printoptions(precision=4)

datatype = u.quantity.Quantity

def niceprint(*args):
  __builtins__.print(*("{0:0.4f}".format(a) if isinstance(a, datatype) or\
      isinstance(a, float) else a for a in args))

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

hbar = h/(2*np.pi)
dBm = u.dB(u.mW)

# Resonator Parameters
gamma_s = 1
Q_int = 3e5
# f_g = 205.128 * u.MHz
f_g = 328.8 * u.MHz
T_amp = 3 * u.Kelvin
eta_read = 0.1
chi_ph = 0.658 # flink factor for phonon noise

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

gamma_leg = 2.65 # conductivity index = beta + 1
K_leg = 120 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.32 * u.Kelvin
T_0 = 0.06 * u.Kelvin # Previously 0.23K Temperature of the thermal bath


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

#Let's try and make an array of P_read
# P_read = (np.linspace(-140,-60,100)*dBm).to(u.pW)
P_read = (-90  * dBm).to(u.pW)

x = P_read/P_opt

# Determine T_b by balancing the input and the output power to the resobolo
T_b= ((((1 + x)* P_opt)/K_leg + T_0**gamma_leg)**(1./gamma_leg)).to('K')
niceprint("The temperature of the island", T_b)
# T_b = 60 * u.mK

# From the Specific heat and thermal conductivity of low-stress amorphous
# Si–N membranes paper by B.L. Zink*, F. Hellman. Si-N membrane-based microcalorimetry:
# Heat capacity and thermal conductivity of thin films
# A_sn = 21 * u.g/u.mol
# rho_sn = 2.9 * u.g/u.cm**3
# T_D = 985 * u.K # Debye temperature of amorphous Si-N
# V_island = 480 * u.um * 150 * u.um * 0.25 * u.um #Assuming 500nm island thickness
# N = (rho_sn * V_island/A_sn) * N_A
# C_b = ((12*np.pi**4/5) * N * k_B * (T_b/T_D)**3).to(u.pJ/u.Kelvin)
# print (C_b)
C_b = 0.09 * u.pJ/u.Kelvin
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
niceprint ("The surface resistance", Rs)
L_k = (h * Rs/(2 * np.pi**2 * Delta) * N_sq).to('nH') # Kinetic inductance contribution
niceprint("Kinetic Inductance per Square", (h * Rs/(2 * np.pi**2 * Delta)).to(u.pH))
niceprint ("The kinetic inductance", L_k)
L = L_g + L_k # total inductance

niceprint ("The total inductance", L)
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
niceprint ("The parameter eta has a value ", eta)
S_1 = ((2/np.pi)*np.sqrt(2*Delta/(np.pi*k_B*T_b))*np.sinh(eta * u.rad)*K_0(eta)).to(1)
S_2 = (1 + np.sqrt(2*Delta/(np.pi*k_B*T_b)) * np.exp(-eta) * I_0(eta)).to(1)
# niceprint (S_1, S_2)
beta = S_2/S_1

N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')
Gamma_gen = (eta_read * P_read/Delta).to('1/s')
n_th = (2*N_0 * np.sqrt(2*np.pi* k_B * T_b* Delta)*np.exp(-Delta/(k_B*T_b))).to('1/um^3')
E_crit = 2 * N_0 * Delta**2 * V_sc
niceprint ("Quasiparticle parameters calculated.\n")

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
# Overwrite to simulate more realistic devices
Q_i = 48494
Q_c = 155298
Q_r = 1./(1./Q_c + 1./Q_i)
P_crit = (0.8 * (2*np.pi*f_r)* E_crit * Q_c/2/Q_r**3).to(u.pW)

niceprint("")
niceprint ("n_qp", n_qp)
niceprint ("n_th ", n_th)


niceprint("")
niceprint ("Q factor from qp losses, Q_qp ", Q_qp)
niceprint ("Resonant Frequency, f_r ", f_r)
niceprint ("Internal Q factor, Q_i ", Q_i)
niceprint ("Coupling Q factor, Q_c ", Q_c)
niceprint ("Resonator Q factor, Q_r ", Q_r)
niceprint ("Kinetic Inductance Fraction, alpha ", alpha)
niceprint ("Beta ", beta)
niceprint ("surface resistance, Rs", Rs)

dx = ((f_g - f_r)/f_r).to(1)
S_21 = 1 - (Q_r/Q_c)* 1./(1 + 2 * 1j * Q_r * dx)
df_r = (f_r/(2*Q_r)).to('kHz')
df_g = (f_g - f_r).to('kHz')

chi_c = (4 * Q_r**2)/(Q_i * Q_c)
chi_g = 1./(1 + (df_g/df_r)**2)
chi_qp = Q_i/Q_qp
P_g = (2/chi_c) * P_read

niceprint("")
niceprint ("Resonator Bandwidth", df_r)
niceprint ("Coupling efficiency", chi_c)
niceprint ("Detuning efficiency", chi_g)
niceprint ("Fraction of Q_i from qp losses", chi_qp)

# Now we include the NEP estimates for the resobolo currently
kappa = (1/2 + Delta/k_B/T_b).to(1)
P_leg = P_opt * (1 + x) # Total power into the resobolo thermal link
gamma_g = (K_leg * gamma_leg * T_b**gamma_leg / P_leg).to(1)
G_b = (gamma_g * P_leg/T_b).to('pW/K') # Conductance of the resobolo
P_b = G_b * T_b
tau_b = (C_b/G_b).to('us')
f_b = (1/tau_b).to(u.kHz) # roll off frequency for bolometer
niceprint("kappa", kappa)

niceprint ("")
niceprint ("dln n_qp/ d ln T_b", kappa)
niceprint ("Conductance of Island ", G_b)
niceprint ("Heat Capacity of Island ", C_b)
niceprint ("Quasiparticle lifetime, tau_qp ", tau_qp)
niceprint ("Equilibrium qp lifetime, tau_th ", tau_th)
niceprint ("Bolometer Time constant", tau_b)
niceprint ("Bolometer Roll off frequency", f_b)

niceprint ("")
niceprint ("The optical power is ", P_opt)
niceprint ("The readout power is ", P_read)
niceprint ("The generated power is ", P_g)
niceprint ("Critical readout power ", P_crit)
niceprint ("The island power, P_b ", P_b)

# Calculating the responsivities on resonance
# r = dx/dPopt
r = (0.5 * (chi_qp * beta/Q_i) * tau_qp/tau_th * kappa/P_b).to(1/u.pW)
r_f = (f_r * r).to(u.kHz/u.pW)

niceprint ("")
niceprint ("resobolo responsivity ignoring bolometer rolloff", r_f)

# Phonon NEP
S_ph = (4 * chi_ph * k_B * T_b**2 * G_b ).to(u.aW**2/u.Hz)

NEP_ph = S_ph ** 0.5

# Amplifier NEP
S_amp = (k_B * T_amp/P_g).to(1/u.Hz)

# NEP_amp = (2 * S_amp**0.5/r/chi_c/Q_i).to(u.aW/u.Hz**0.5)
NEP_amp = (2 * S_amp**0.5 * Q_c/Q_r**2/r).to(u.aW/u.Hz**0.5)


# Shot NEP
S_shot = 2 * n_qp * V_sc * (1/tau_max + 1/tau_qp)

NEP_gr = (S_shot**0.5 * (tau_th/(n_qp * V_sc) * P_b/kappa)).to(u.aW/u.Hz**0.5)

# Total NEP
NEP_total = (NEP_gr**2 + NEP_amp**2 + NEP_ph**2)**0.5

niceprint ("")
niceprint ("Phonon NEP", NEP_ph)
niceprint ("Amplifier NEP", NEP_amp)
niceprint ("Shot NEP", NEP_gr)
niceprint ("Total NEP", NEP_total)

nu = np.logspace(-1, 6, 5000) * u.Hz

# I'll ignore the bolometer roll off factor for now
#H = np.ones_like(nu.value)
ones = np.ones_like(nu.value)
H =  1/(1 + 1j * 2 * np.pi * nu * tau_b)
H = np.abs(H)

NEP_ph = NEP_ph * ones
NEP_amp = NEP_amp * H
NEP_gr = NEP_gr * H
NEP_total = (NEP_gr**2 + NEP_amp**2 + NEP_ph**2)**0.5

fig, ax = plt.subplots(figsize=(10,10))
ax.loglog(nu, NEP_ph, 'b', label='Phonon')
ax.loglog(nu, NEP_amp, 'r--', label='Amplifier')
ax.loglog(nu, NEP_gr, 'g', label='Gen-Recomb')
ax.loglog(nu, NEP_total, 'k', label='Total')

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'NEP [aW/rtHz]')
ax.set_ylim([1, 1000])
ax.set_xlim([0.1, 1e6])
ax.grid(which='major', axis='both')
ax.legend(loc='best')

plt.savefig(r'NEP.png')
