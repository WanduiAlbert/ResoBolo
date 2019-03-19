#! /usr/bin/env python3

# My version of the resonator code. Most of the calculations are adapted from the ipython notebook that I have in this folder.
# I initialize the resonator bolometer with some of its parameters so I can easily adjust it.

import numpy as np
from math import pi
import astropy.units as u
from astropy.constants import h,k_B, c, m_e, N_A
import matplotlib.pyplot as plt
from scipy.special import kn, iv
import pdb

np.set_printoptions(precision=4)
verbose = False
datatype = u.quantity.Quantity

def niceprint(*args):
	print(*("{0:0.4f}".format(a) if isinstance(a, datatype) or\
		isinstance(a, float) else a for a in args))

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

hbar = h/(2*np.pi)
dB = u.dB(1)
dBm = u.dB(u.mW)

# Resonator Parameters
gamma_s = 1
Q_int = 3e8
# f_g = 205.128 * u.MHz

T_amp = 5.22 * u.Kelvin
eta_read = 0.1
chi_ph = 0.658 # flink factor for phonon noise

# I'll keep these here for now but I'm more interested in doing calculations for a dark run.
nu_opt = 150 * u.GHz
dnu_opt = 0.25
eta_opt = 0.25
T_rj = 31* u.Kelvin #0.1 * u.Kelvin
N_pol = 1
P_opt = 14.00*u.pW #((Vdc/Rb)**2 * Rh).to(u.pW)

gamma_leg = 3.131 # conductivity index = beta + 1
K_leg = 403.705 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.329 * u.Kelvin
T_0 = 0.25 * u.Kelvin # Previously 0.23K Temperature of the thermal bath

# Material properties of the Aluminum superconductor
tau_max = 500 * u.microsecond
n_qp_star = 100 * 1/u.um**3
gamma = 1.35 * u.mJ/u.mol/u.Kelvin**2
density = 2.7 * u.g/u.cm**3
A_r = 26.98 * u.g/u.mol
rho = 1.15 * u.uOhm * u.cm

#T_b = 0.38 * u.Kelvin # Temperature of the bolometer
# Calculated parameters

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
# C_b = 0.9 * u.pJ/u.Kelvin

Delta = (1.764 * k_B * T_c).to('J')

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
L_k = (h * Rs/(2 * np.pi**2 * Delta) * N_sq).to('nH') # Kinetic inductance contribution
Z0 = 50 * u.Ohm # Characteristic impedance of the line

class OperatingPoint:
	"""
	Defines the operating point of a TKID. Used to model the noise terms
	for the bolometer.
	"""


	def __init__(self, P_opt=0*u.pW, f_r=305.8*u.MHz, T_0=250*u.mK,\
		P_read=-81*dBm):
		self.f_g = f_r #We will work exactly on resonance
		self.f_r = f_r
		self.T_0 = T_0 #Temperature of the thermal bath
		self.P_read = P_read.to(u.pW)
		self.P_opt = P_opt
		self.K_leg = 0
		self.gamma_leg = 1
		self.C_c1 = 0.3984 * u.pF
		self.C_c = self.C_c1/2

		print ("Resonator parameters set.")

	def calculate_noise(self):
		self.x = 0.5*self.P_read/self.P_opt
		self.P_total = 0.5*self.P_read + self.P_opt
		#pdb.set_trace()
		self.T_b= (((self.P_total/self.K_leg).si.value +
			self.T_0.si.value**self.gamma_leg)**(1./self.gamma_leg))*u.Kelvin
		niceprint("The temperature of the island", self.T_b)

		self.C_b = 0.25*(self.T_b/(350*u.mK))**2*u.pJ/u.Kelvin #Heat capacity


		self.L = L_g + L_k # total inductance
		self.alpha = (L_k/self.L).to(1)
		if verbose: niceprint ("The total inductance", self.L)

		# f_r = 329 * u.MHz
		self.omega_r = (2*np.pi*self.f_r).to(1/u.s)
		self.C = (1/(self.L * self.omega_r**2)).to(u.pF)

		# C_i = 27.98 * u.pF
		#C_c1 = 0.19 * u.pF

		# C = C_c1 + C_i
		self.C_i = self.C - self.C_c
		# omega_r = (1./np.sqrt(L*C)).to('1/s')
		# f_r = (omega_r/(2*np.pi)).to('MHz')

		self.eta = (h * self.f_g / (2*k_B * self.T_b)).to(1).value # Weird conversion because astropy
		# niceprint ("The parameter eta has a value ", eta)
		S_1 = ((2/np.pi)*np.sqrt(2*Delta/(np.pi*k_B*self.T_b))*np.sinh(self.eta * u.rad)*K_0(self.eta)).to(1)
		S_2 = (1 + np.sqrt(2*Delta/(np.pi*k_B*self.T_b)) * np.exp(-self.eta) * I_0(self.eta)).to(1)
		# niceprint (S_1, S_2)
		self.beta = S_2/S_1

		N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')
		Gamma_gen = (eta_read * self.P_read/Delta).to('1/s')
		self.n_th = (2*N_0 * np.sqrt(2*np.pi* k_B * self.T_b* Delta)*np.exp(-Delta/(k_B*self.T_b))).to('1/um^3')
		E_crit = 2 * N_0 * Delta**2 * V_sc
		niceprint ("Quasiparticle parameters calculated.\n")

		# Quality factors
		Gamma_th = ((self.n_th * V_sc/tau_max)* (1 + 0.5 * self.n_th/n_qp_star)).to('1/s')
		self.n_qp = (np.sqrt((self.n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star).to('1/um^3')
		self.tau_th = (tau_max/(1 + self.n_th/n_qp_star)).to('us')
		self.tau_qp = (tau_max/(1 + self.n_qp/n_qp_star)).to('us')
		self.Q_qp = ((2 * N_0 * Delta)/(self.alpha * gamma_s * S_1 * self.n_qp)).to(1)
		self.Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * self.T_b))/np.sinh(self.eta)/K_0(self.eta)
		self.Q_c = (2 * self.C_i/(self.omega_r * self.C_c**2 * Z0)).to(1)
		self.Q_i = 1./(1/self.Q_qp + 1./Q_int)
		self.Q_r  = 1./(1./self.Q_c + 1./self.Q_i)
		self.P_crit = (0.8 * (self.omega_r)* E_crit * self.Q_c/2/self.Q_r**3).to(u.pW)

		if verbose:
			niceprint ("Q factor from qp losses, Q_qp ", self.Q_qp)
			niceprint ("Resonant Frequency, f_r ", self.f_r)
			niceprint ("Internal Q factor, Q_i ", self.Q_i)
			niceprint ("Coupling Q factor, Q_c ", self.Q_c)
			niceprint ("Resonator Q factor, Q_r ", self.Q_r)
			niceprint ("Kinetic Inductance Fraction, alpha ", self.alpha)
			niceprint ("Beta ", self.beta)
			niceprint ("surface resistance, Rs", Rs)

		self.dx = ((self.f_g - self.f_r)/self.f_r).to(1)
		self.S_21 = 1 - (self.Q_r/self.Q_c)* 1./(1 + 2 * 1j * self.Q_r * self.dx)
		self.df_r = (self.f_r/(2*self.Q_r)).to('kHz')
		self.df_g = (self.f_g - self.f_r).to('kHz')

		self.chi_c = (4 * self.Q_r**2)/(self.Q_i * self.Q_c)
		self.chi_g = 1./(1 + (self.df_g/self.df_r)**2)
		self.chi_qp = self.Q_i/self.Q_qp
		self.P_g = (2/self.chi_c) * self.P_read

		if verbose:
			niceprint("")
			niceprint ("Resonator Bandwidth", self.df_r)
			niceprint ("Coupling efficiency", self.chi_c)
			niceprint ("Detuning efficiency", self.chi_g)
			niceprint ("Fraction of Q_i from qp losses", self.chi_qp)

		# Now we include the NEP estimates for the resobolo currently
		self.kappa = (1/2 + Delta/k_B/self.T_b).to(1)
		self.P_leg = self.P_total # Total power into the resobolo thermal link
		gamma_g = (self.K_leg * self.gamma_leg * self.T_b**self.gamma_leg / self.P_leg).to(1)
		#self.G_b = (gamma_g * self.P_leg/self.T_b).to('pW/K') # Conductance of the resobolo
		self.G_b = (self.gamma_leg*self.K_leg*self.T_b**(self.gamma_leg-1)).to('pW/K') # Conductance of the resobolo
		self.P_b = self.G_b * self.T_b
		self.tau_b = (self.C_b/self.G_b).to('ms')
		self.f_b = (1/2/pi/self.tau_b).to(u.Hz) # roll off frequency for bolometer

		if verbose:
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
		self.r = (0.5 * (self.chi_qp * self.beta/self.Q_i) *\
			self.tau_qp/self.tau_th * self.kappa/self.P_b).to(1/u.pW)
		self.r_f = (self.f_r * self.r).to(u.kHz/u.pW)

		niceprint ("")
		niceprint ("resobolo responsivity ignoring bolometer rolloff", self.r_f)

	def calculate_noise_spectra(self):

		# Phonon NEP
		self.S_ph = (4 * chi_ph * k_B * self.T_b**2 * self.G_b ).to(u.aW**2/u.Hz)

		self.NEP_ph = self.S_ph ** 0.5

		# Amplifier NEP
		self.S_amp = (k_B * T_amp/self.P_g).to(1/u.Hz)

		self.NEP_amp = (2 * self.S_amp**0.5/self.r/self.chi_c/self.Q_i).to(u.aW/u.Hz**0.5)

		# Shot NEP
		self.S_shot = 2 * self.n_qp * V_sc * (1/tau_max + 1/self.tau_qp)

		self.NEP_gr = (self.S_shot**0.5 * (self.tau_th/(self.n_qp * V_sc) * self.P_b/self.kappa)).to(u.aW/u.Hz**0.5)

		# Total NEP
		self.NEP_total = (self.NEP_gr**2 + self.NEP_amp**2 + self.NEP_ph**2)**0.5

		niceprint ("")
		niceprint ("Phonon NEP", self.NEP_ph)
		niceprint ("Amplifier NEP", self.NEP_amp)
		niceprint ("Shot NEP", self.NEP_gr)
		niceprint ("Total NEP", self.NEP_total)


	def plot_noise(self, ax, fmax=6):
		self.calculate_noise_spectra()
		self.nu = np.logspace(-1, fmax, 5000) * u.Hz
		self.ones = np.ones_like(self.nu.value)
		self.bolo_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_b).to(1))
		self.qp_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_qp).to(1))
		self.reso_rolloff =  1/(1 + 1j * (self.nu / self.df_r).to(1))
		self.bolo_rolloff = np.abs(self.bolo_rolloff)
		self.qp_rolloff = np.abs(self.qp_rolloff)
		self.reso_rolloff = np.abs(self.reso_rolloff)

		self.NEP_ph_spec = self.NEP_ph * self.bolo_rolloff
		self.NEP_amp_spec = self.NEP_amp * self.ones
		self.NEP_gr_spec = self.NEP_gr * self.qp_rolloff * self.reso_rolloff
		self.NEP_total_spec = (self.NEP_gr_spec**2 + self.NEP_amp_spec**2 +
				self.NEP_ph_spec**2)**0.5

		ax.loglog(self.nu, self.NEP_ph_spec, 'b', label='Phonon')
		ax.loglog(self.nu, self.NEP_amp_spec, 'r--', label='Amplifier')
		ax.loglog(self.nu, self.NEP_gr_spec, 'g', label='Gen-Recomb')
		ax.loglog(self.nu, self.NEP_total_spec, 'k', label='Total')


if __name__=="__main__":
	# We are using the heater pad to mimic the optical load on the bolometer
	Rh = 0.150 * u.Ohm
	frequencies = np.array([305.8, 318.4, 337.4])*u.MHz
	Rb = np.array([300, 404, 510]) * u.kOhm
	Pref = 10e-12
	Vmax = (Pref/Rh.value)**0.5*510e3
	Npowers = 10
	Vdc = np.linspace(0, Vmax, Npowers)*u.Volt
	P_opt = ((Rh/Rb[:, np.newaxis]**2)*Vdc**2).to(u.pW)
	print (P_opt.shape)
	gamma_leg = np.array([3.131, 2.918, 3.031])# conductivity index = beta + 1
	K_leg = np.array([403.705, 186.552, 141.179])* u.picoWatt
	Cc = np.array([0.3984, 0.3731, 0.3478])*u.pF
	T_c = 1.329 * u.Kelvin
	T_0 = 0.25 * u.Kelvin # Previously 0.23K Temperature of the thermal bath
	P_read = -81*dBm
	P_read -= 20*np.log10(4)*dB # Accounting for the 4 resonators being read at the
	# same time
	#P_opt = 10*u.pW
	NEPs = np.zeros((len(frequencies), Npowers))
	NEPs_ph = np.zeros((len(frequencies), Npowers))
	NEPs_amp = np.zeros((len(frequencies), Npowers))
	for ireso, reso in enumerate(frequencies):
		op = OperatingPoint(P_opt=0*u.pW, f_r=reso, T_0=T_0, P_read=P_read)
		op.C_c1 = Cc[ireso]
		op.gamma_leg = gamma_leg[ireso]
		op.K_leg = K_leg[ireso]/u.Kelvin**gamma_leg[ireso]
		for ipow in range(Npowers):
			#if ipow != 6: continue
			op.P_opt = P_opt[ireso, ipow]
			op.calculate_noise()
			op.calculate_noise_spectra()
			NEPs[ireso, ipow] = op.NEP_total.value
			NEPs_ph[ireso, ipow] = op.NEP_ph.value
			NEPs_amp[ireso, ipow] = op.NEP_amp.value
			fig, ax = plt.subplots(figsize=(10,10))
			op.plot_noise(ax)
			ax.set_title("%.1f MHz reso biased at %.1f pW"%(frequencies[ireso].value, P_opt[ireso, ipow].value))
			ax.set_xlabel(r'Frequency [Hz]')
			ax.set_ylabel(r'NEP [aW/rtHz]')
			ax.set_ylim([1, 1000])
			ax.set_xlim([0.1, 1e6])
			ax.grid(which='major', axis='both')
			ax.legend(loc='best')
			plt.savefig('Waffle_TKID_%.1fMHz_%.1fpW_Noise_Prediction.png'%(frequencies[ireso].value, P_opt[ireso, ipow].value))
			plt.close()
			#op.calculate_noise_spectra()
	#op.calculate_noise()
	#fig, ax = plt.subplots(figsize=(10,10))
	#op.plot_noise(ax)
	#ax.set_xlabel(r'Frequency [Hz]')
	#ax.set_ylabel(r'NEP [aW/rtHz]')
	#ax.set_ylim([1, 1000])
	#ax.set_xlim([0.1, 1e6])
	#ax.grid(which='major', axis='both')
	#ax.legend(loc='best')
	#plt.savefig('Waffle_TKID_305.8MHz_Noise_Prediction.png')
	exit()
	P_readout = (1*u.mW*10**(P_read.value/10)).to(u.pW)
	print (P_readout)
	for ireso, reso in enumerate(frequencies):
		fig, ax = plt.subplots(figsize=(10,10))
		ax.plot(P_opt[ireso].value, NEPs[ireso], label='total')
		#ax.plot(P_opt[ireso].value, NEPs_ph[ireso], 'k--', label='phonon')
		#ax.plot(P_opt[ireso].value, NEPs_amp[ireso], 'k-.', label='amplifier')
		ax.set_title("%.1f MHz reso"%(frequencies[ireso].value))
		ax.set_xlabel(r'Optical Power [pW]')
		ax.set_ylabel(r'NEP [aW/rtHz]')
		ax.legend(loc="lower right")
		#ax.set_ylim([1, 1000])
		#ax.set_xlim([0.1, 1e6])
		ax.grid(which='major', axis='both')
		plt.savefig('Waffle_TKID_%.1fMHz_NEP_vs_Power.png'%(frequencies[ireso].value))
		plt.close()
