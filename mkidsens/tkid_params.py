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
Q_int = 187519
# f_g = 205.128 * u.MHz

T_amp = 3 * u.Kelvin
eta_read = 0.1
chi_ph = 0.658 # flink factor for phonon noise

# I'll keep these here for now but I'm more interested in doing calculations for a dark run.
nu_opt = 150 * u.GHz
dnu_opt = 0.25
eta_opt = 0.25
T_rj = 31* u.Kelvin #0.1 * u.Kelvin
N_pol = 1

gamma_leg = 2.65 # conductivity index = beta + 1
K_leg = 120 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.32 * u.Kelvin
T_0 = 0.06 * u.Kelvin 

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


    def __init__(self, P_opt=0*u.pW, f_r=329*u.MHz, T_0=60*u.mK,\
        P_read=-110*dBm):
        self.f_g = f_r #We will work exactly on resonance
        self.f_r = f_r
        self.T_0 = T_0 #Temperature of the thermal bath
        self.P_read = P_read.to(u.pW)
        self.P_opt = P_opt

        print ("Resonator parameters set.")

        self.x = self.P_read/self.P_opt
        self.P_total = self.P_read + self.P_opt
        self.T_b= ((self.P_total/K_leg + self.T_0**gamma_leg)**(1./gamma_leg)).to('K')
        niceprint("The temperature of the island", self.T_b)

        self.C_b = 0.9*u.pJ/u.Kelvin #Heat capacity


        self.L = L_g + L_k # total inductance
        self.alpha = (L_k/self.L).to(1)
        niceprint ("The total inductance", self.L)
        
        # f_r = 329 * u.MHz
        self.omega_r = (2*np.pi*self.f_r).to(1/u.s)
        self.C = (1/(self.L * self.omega_r**2)).to(u.pF)

        # C_i = 27.98 * u.pF
        C_c1 = 0.19 * u.pF
        C_c = C_c1 # The same now but will change once I add capacitance to ground

        # C = C_c1 + C_i
        self.C_i = self.C - C_c
        # omega_r = (1./np.sqrt(L*C)).to('1/s')
        # f_r = (omega_r/(2*np.pi)).to('MHz')

        self.eta = (h * self.f_g / (2 * k_B * self.T_b)).to(1).value # Weird conversion because astropy
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
        self.Q_c = (2 * self.C_i/(self.omega_r * C_c**2 * Z0)).to(1)
        self.Q_i = 1./(1/self.Q_qp + 1./Q_int)
        self.Q_r  = 1./(1./self.Q_c + 1./self.Q_i)
        self.P_crit = (0.8 * (self.omega_r)* E_crit * self.Q_c/2/self.Q_r**3).to(u.pW)

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

        niceprint("")
        niceprint ("Resonator Bandwidth", self.df_r)
        niceprint ("Coupling efficiency", self.chi_c)
        niceprint ("Detuning efficiency", self.chi_g)
        niceprint ("Fraction of Q_i from qp losses", self.chi_qp)

        # Now we include the NEP estimates for the resobolo currently
        self.kappa = (1/2 + Delta/k_B/self.T_b).to(1)
        self.P_leg = self.P_opt * (1 + self.x) # Total power into the resobolo thermal link
        gamma_g = (K_leg * gamma_leg * self.T_b**gamma_leg / self.P_leg).to(1)
        self.G_b = (gamma_g * self.P_leg/self.T_b).to('pW/K') # Conductance of the resobolo
        self.P_b = self.G_b * self.T_b
        self.tau_b = (self.C_b/self.G_b).to('us')
        self.f_b = (1/self.tau_b).to(u.kHz) # roll off frequency for bolometer


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

    def calculate_noise(self):

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

        

    def plot_noise(self, ax):
        calculate_noise(self)
        nu = np.logspace(-1, 6, 5000) * u.Hz

        # I'll ignore the bolometer roll off factor for now
        #H = np.ones_like(nu.value)
        ones = np.ones_like(nu.value)
        H =  1/(1 + 1j * 2 * np.pi * nu * tau_b)
        H = np.abs(H)

        NEP_ph = self.NEP_ph * ones
        NEP_amp = self.NEP_amp * H
        NEP_gr = self.NEP_gr * H
        NEP_total = (NEP_gr**2 + NEP_amp**2 + NEP_ph**2)**0.5

        # fig, ax = plt.subplots(figsize=(10,10))
        ax.loglog(nu, NEP_ph, 'b', label='Phonon')
        ax.loglog(nu, NEP_amp, 'r--', label='Amplifier')
        ax.loglog(nu, NEP_gr, 'g', label='Gen-Recomb')
        ax.loglog(nu, NEP_total, 'k', label='Total')

        # ax.set_xlabel(r'Frequency [Hz]')
        # ax.set_ylabel(r'NEP [aW/rtHz]')
        # ax.set_ylim([1, 1000])
        # ax.set_xlim([0.1, 1e6])
        # ax.grid(which='major', axis='both')
        # ax.legend(loc='best')

        # plt.savefig(r'NEP.png')


if __name__=="__main__":
    # We are using the heater pad to mimic the optical load on the bolometer
    Rh = 0.106 * u.Ohm
    Rb = 854 * u.kOhm
    Vdc = 7.0 * u.Volt

    P_opt = ((Vdc/Rb)**2 * Rh).to(u.pW)
    op = OperatingPoint(P_opt=P_opt)
    op.calculate_noise()
    
