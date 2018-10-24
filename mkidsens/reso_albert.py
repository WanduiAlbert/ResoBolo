#! /usr/bin/env python3

# My version of the resonator code. Most of the calculations are adapted from
# the ipython notebook that I have in this folder.
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
f_g = 469.9 * u.MHz
T_amp = 5.22 * u.Kelvin
eta_read = 0.1
chi_ph = 0.658 # flink factor for phonon noise

# I'll keep these here for now but I'm more interested in doing calculations for a dark run.
nu_opt = 150 * u.GHz
dnu_opt = 0.25
eta_opt = 0.25
T_rj = 31* u.Kelvin #0.1 * u.Kelvin
N_pol = 1

# We are using the heater pad to mimic the optical load on the bolometer
#Rh = 0.106 * u.Ohm
#Rb = 854 * u.kOhm
#Vdc = 7.0 * u.Volt

P_opt = 5.00*u.pW #((Vdc/Rb)**2 * Rh).to(u.pW)

gamma_leg = 2.975 # conductivity index = beta + 1
K_leg = 120.660 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.385 * u.Kelvin
T_0 = 0.08 * u.Kelvin # Previously 0.23K Temperature of the thermal bath

# Material properties of the Aluminum superconductors
tau_max = 500 * u.microsecond
n_qp_star = 100 * 1/u.um**3
gamma = 1.35 * u.mJ/u.mol/u.Kelvin**2
density = 2.7 * u.g/u.cm**3
A_r = 26.98 * u.g/u.mol
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')

rho = 1.15 * u.uOhm * u.cm
L_g = 6.1 * u.nH # From simulations
Z0 = 50 * u.Ohm # Characteristic impedance of the line
C_c1 = 0.2878 * u.pF
C_c = C_c1/2.

class TKIDBolometer:

	def init(self, fr, f_g, Tc,  Cc, T_0,
			P_opt, P_read, K_leg, gamma_leg, L_g=6.1, **kwargs):
		self.C_b = 0.14 * u.pJ/u.Kelvin # From TKID paper
		self.fr = fr * u.MHz
		self.f_g = f_g * u.MHz
		self.Tc = Tc
		self.L_g = L_g * u.nH
		self.Cc1 = Cc
		self.C_c = self.C_c/2
		self.T_0 = T_0
		self.P_opt = P_opt * u.pW
		self.P_read = (P_read * dBm).to(u.pW)
		self.gamma_leg = gamma_leg
		self.K_leg = K_leg * u.picoWatt/u.Kelvin**self.gamma_leg
		try:
			self.t = kwargs['t'] *u.um
			self.w_trace = kwargs['w_trace'] *u.um
			self.s_trace = kwargs['s_trace'] *u.um
			self.N_sq = kwargs['N_sq']
		except KeyError:
			print ("Using the default inductor dimensions")
			self.t = 0.05 *u.um
			self.w_trace = 1 * u.um
			self.s_trace = 1 * u.um
			self.N_sq = 16200
		print ("Resonator parameters set.")
		self.Delta = (1.764 * k_B * T_c).to('J')

	def operating_point(self):
		self.x = self.P_read/self.P_opt
		self.T_b= ((((1 + self.x)* self.P_opt)/self.K_leg +
			self.T_0**self.gamma_leg)**(1./self.gamma_leg)).to('K')
		niceprint("The temperature of the island", self.T_b)
		self.l = self.N_sq * self.w_trace
		self.A = self.t * self.w_trace # cross-sectional area
		self.V_sc = self.l * self.A

		self.Rs = (rho/self.t).to('Ohm') # surface resistance in ohms/sq
		niceprint ("The surface resistance", self.Rs)
		self.L_k_psq = (self.h * self.Rs/(2 * np.pi**2 * self.Delta)).to(u.pH)
		self.L_k =  (self.L_k_psq * self.N_sq).to('nH') # Kinetic inductance contribution
		niceprint("Kinetic Inductance per Square", self.L_k_psq)
		niceprint ("The kinetic inductance", self.L_k)
		self.L = self.L_g + self.L_k # total inductance

		niceprint ("The total inductance", self.L)
		self.alpha = (self.L_k/self.L).to(1)

		self.omega_r = (2*np.pi*self.f_r).to(1/u.s)
		self.C = (1/(self.L * self.omega_r**2)).to(u.pF)

		self.C_i = self.C - self.C_c

		self.eta = (h * self.f_g / (2 * k_B * self.T_b)).to(1).value # Weird conversion because astropy
		niceprint ("The parameter eta has a value ", self.eta)
		self.S_1 = ((2/np.pi)*np.sqrt(2*self.Delta/(np.pi*k_B*self.T_b))*
				np.sinh(self.eta * u.rad)*K_0(self.eta)).to(1)
		self.S_2 = (1 + np.sqrt(2*self.Delta/(np.pi*k_B*self.T_b))*
				np.exp(-self.eta) * I_0(self.eta)).to(1)
		# niceprint (S_1, S_2)
		self.beta = self.S_2/self.S_1

		self.Gamma_gen = (eta_read * self.P_read/self.Delta).to('1/s')
		self.n_th = (2*N_0 * np.sqrt(2*np.pi* k_B * self.T_b* self.Delta)*np.exp(-self.Delta/(k_B*self.T_b))).to('1/um^3')
		self.E_crit = 2 * N_0 * self.Delta**2 * self.V_sc
		niceprint ("Quasiparticle parameters calculated.\n")

		# Quality factors
		self.Gamma_th = ((self.n_th * self.V_sc/tau_max)* (1 + 0.5 * self.n_th/n_qp_star)).to('1/s')
		self.n_qp = (np.sqrt((self.n_th + n_qp_star)**2 + (2*self.Gamma_gen*n_qp_star*tau_max/self.V_sc)) - n_qp_star).to('1/um^3')
		self.tau_th = (tau_max/(1 + self.n_th/n_qp_star)).to('us')
		self.tau_qp = (tau_max/(1 + self.n_qp/n_qp_star)).to('us')
		self.Q_qp = ((2 * N_0 * self.Delta)/(self.alpha * gamma_s * self.S_1 * self.n_qp)).to(1)
		self.Q_sigma = (np.pi/4)*np.exp(self.Delta/(k_B * self.T_b))/np.sinh(self.eta)/K_0(self.eta)
		self.Q_c = (2 * self.C_i/(self.omega_r * self.C_c**2 * Z0)).to(1)
		self.Q_i = 1./(1/self.Q_qp + 1./Q_int)
		self.Q_r  = 1./(1./self.Q_c + 1./self.Q_i)
		# Overwrite to simulate more realistic devices
		#Q_i = 48494
		#Q_c = 155298
		#Q_r = 1./(1./Q_c + 1./Q_i)
		self.P_crit = (0.8 * (2*np.pi*self.f_r)* self.E_crit * self.Q_c/2/self.Q_r**3).to(u.pW)

		niceprint("")
		niceprint ("n_qp", self.n_qp)
		niceprint ("n_th ", self.n_th)


		niceprint("")
		niceprint ("Q factor from qp losses, Q_qp ", self.Q_qp)
		niceprint ("Resonant Frequency, f_r ", self.f_r)
		niceprint ("Internal Q factor, Q_i ", self.Q_i)
		niceprint ("Coupling Q factor, Q_c ", self.Q_c)
		niceprint ("Resonator Q factor, Q_r ", self.Q_r)
		niceprint ("Kinetic Inductance Fraction, alpha ", self.alpha)
		niceprint ("Beta ", self.beta)
		niceprint ("surface resistance, Rs", self.Rs)

		self.dx = ((self.f_g - self.f_r)/self.f_r).to(1)
		self.S_21 = 1 - (self.Q_r/self.Q_c)* 1./(1 + 2 * 1j * self.Q_r * self.dx)
		self.df_r = (self.f_r/(2*self.Q_r)).to('kHz')
		self.df_g = (self.f_g - self.f_r).to('kHz')

		self.chi_c = (4 * self.Q_r**2)/(self.Q_i * self.Q_c)
		self.chi_g = 1./(1 + (self.df_g/self.df_r)**2)
		self.chi_qp = self.Q_i/self.Q_qp
		P_g = (2/chi_c) * P_read

		niceprint("")
		niceprint ("Resonator Bandwidth", self.df_r)
		niceprint ("Coupling efficiency", self.chi_c)
		niceprint ("Detuning efficiency", self.chi_g)
		niceprint ("Fraction of Q_i from qp losses", self.chi_qp)

		# Now we include the NEP estimates for the resobolo currently
		self.kappa = (1/2 + self.Delta/k_B/self.T_b).to(1)
		self.P_leg = self.P_opt * (1 + self.x) # Total power into the resobolo thermal link
		self.gamma_g = (self.K_leg * self.gamma_leg *
				self.T_b**self.gamma_leg / self.P_leg).to(1)
		self.G_b = (self.gamma_g * self.P_leg/self.T_b).to('pW/K') # Conductance of the resobolo
		self.P_b = self.G_b * self.T_b
		self.tau_b = (self.C_b/self.G_b).to('us')
		self.f_b = (1/self.tau_b).to(u.Hz) # roll off frequency for bolometer

		niceprint ("")
		niceprint ("dln n_qp/ d ln T_b", self.kappa)
		niceprint ("Conductance of Island ", self.G_b)
		niceprint ("Heat Capacity of Island ", self.C_b)
		niceprint ("Quasiparticle lifetime, tau_qp ", self.tau_qp)
		niceprint ("Equilibrium qp lifetime, tau_th ", self.tau_th)
		niceprint ("Bolometer Time constant", self.tau_b)
		niceprint ("Bolometer Roll off frequency", self.f_b)

		niceprint ("")
		niceprint ("The optical power is ", self.P_opt)
		niceprint ("The readout power is ", self.P_read)
		niceprint ("The generated power is ", self.P_g)
		niceprint ("Critical readout power ", self.P_crit)
		niceprint ("The island power, P_b ", self.P_b)

# Calculating the responsivities on resonance
# r = dx/dPopt
r = (0.5 * (chi_qp * beta/Q_i) * tau_qp/tau_th * kappa/P_b).to(1/u.pW)
r_f = (f_r * r).to(u.kHz/u.pW)

niceprint ("")
niceprint ("resobolo responsivity ignoring bolometer rolloff", r_f)

# Phonon NEP
S_ph = (4 * chi_ph * k_B * self.T_b**2 * G_b ).to(u.aW**2/u.Hz)

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
bolo_rolloff =  1/(1 + 1j * 2 * np.pi * nu * tau_b)
qp_rolloff =  1/(1 + 1j * 2 * np.pi * nu * tau_qp)
bolo_rolloff = np.abs(bolo_rolloff)
qp_rolloff = np.abs(qp_rolloff)

NEP_ph = NEP_ph * bolo_rolloff
NEP_amp = NEP_amp * ones
NEP_gr = NEP_gr * qp_rolloff * bolo_rolloff
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
