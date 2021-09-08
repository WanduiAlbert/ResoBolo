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
verbose = True
datatype = u.quantity.Quantity
useBoloMat = False
if useBoloMat:
    plot_dir = 'resobolo_plots/'
else:
    plot_dir = 'waffle_plots/'

def niceprint(*args):
    print(*("{0:0.4f}".format(a) if isinstance(a, datatype) or\
        isinstance(a, float) else a for a in args))

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

hbar = h/(2*np.pi)
dB = u.dB(1)
dBm = u.dB(u.mW)

T_debye = 433*u.K#K
Omega_debye = (k_B*T_debye).to(u.J)


# Resonator Parameters
gamma_s = 1
Q_int = 3e8
# f_g = 205.128 * u.MHz

T_amp = 5.22 * u.Kelvin
eta_read = 0.0001
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
#tau_max = 955.1 * u.microsecond
#tau_max = 500 * u.microsecond
#n_qp_star = 100 * 1/u.um**3
#tau_max = 688.1 * u.microsecond
#n_qp_star = 663 * 1/u.um**3
tau_max = 500.1 * u.microsecond
#n_qp_star = 100 * 1/u.um**3
#R = 2.1 *u.um**3 / u.s
n_qp_star = 518 * 1/u.um**3
#tau_max = 1./(R * n_qp_star)
R = 1./(tau_max * n_qp_star)
#tau_max = 488.1 * u.microsecond
#R *= 0.5
gamma = 1.35 * u.mJ/u.mol/u.Kelvin**2
density = 2.7 * u.g/u.cm**3
A_r = 26.98 * u.g/u.mol
rho = 1.15 * u.uOhm * u.cm
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')

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

dNV = np.log(2*Omega_debye/Delta)

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
#L_k = 0.851 * L_g
Z0 = 50 * u.Ohm # Characteristic impedance of the line
C_load = 7*u.pF

class OperatingPoint:
    """
    Defines the operating point of a TKID. Used to model the noise terms
    for the bolometer.
    """


    def __init__(self, P_opt=0*u.pW, f_r=305.8*u.MHz, f_g=305.8*u.MHz, T_0=250*u.mK,\
        P_g=-81*dBm, nu0 = 150*u.GHz, opt_bw = 0.25, is_optical=False):
        self.f_g = f_g #We will work exactly on resonance
        self.omega_g = (2*pi*self.f_g).to(1/u.s)
        self.f_r = f_r
        self.T_0 = T_0 #Temperature of the thermal bath
        self.P_g = P_g.to(u.pW)
        self.P_opt = P_opt
        self.K_leg = 0
        self.gamma_leg = 1
        self.C_c1 = 0.3984 * u.pF
        self.C_c = self.C_c1/2
        self.is_optical = is_optical
        self.nu0 = nu0
        self.opt_bw = opt_bw

        print ("Resonator parameters set.")

    def calculate_noise(self):
        self.x = 0.5*self.P_g/self.P_opt
        self.P_total = 0.5*self.P_g + self.P_opt
        #pdb.set_trace()
        self.T_b= (((self.P_total/self.K_leg).value +
            self.T_0.value**self.gamma_leg)**(1./self.gamma_leg))*u.Kelvin
        niceprint("The temperature of the island", self.T_b)

        self.C_b = (0.25*(self.T_b/(350*u.mK))**2*u.pJ/u.Kelvin).to(u.pJ/u.Kelvin) #Heat capacity

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

        Gamma_gen = (eta_read * self.P_g/Delta).to('1/s')
        self.n_th = (2*N_0 * np.sqrt(2*np.pi* k_B * self.T_b* Delta)*np.exp(-Delta/(k_B*self.T_b))).to('1/um^3')
        E_crit = 0.5 * N_0 * Delta**2 * V_sc/self.alpha**2
        #E_crit = 0.5/dNV * Delta**2/self.alpha**2
        niceprint ("Quasiparticle parameters calculated.\n")

        # Quality factors
        Gamma_th = ((self.n_th * V_sc/tau_max)* (1 + 0.5 * self.n_th/n_qp_star)).to('1/s')
        self.n_qp = (np.sqrt((self.n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star).to('1/um^3')
        self.tau_th = (1./(R*self.n_th)).to('us')
        self.tau_qp = (1./(R*self.n_qp)).to('us')

        #self.tau_th = (tau_max/(1 + self.n_th/n_qp_star)).to('us')
        #self.tau_qp = (tau_max/(1 + self.n_qp/n_qp_star)).to('us')
        if verbose:
            niceprint ("")
            niceprint ("Thermal quasiparticle density, n_th ", self.n_th)
            niceprint ("Equilibrium qp lifetime, tau_th ", self.tau_th)
            niceprint ("Total quasiparticle density, n_qp ", self.n_qp)
            niceprint ("Quasiparticle lifetime, tau_qp ", self.tau_qp)
            niceprint ("QP rolloff frequency, f3dB_qp ", (1./(2*pi*self.tau_qp)).to(u.kHz))
        niceprint ("Quasiparticle parameters calculated.\n")

        self.xMB = ((self.alpha * gamma_s * S_2 * self.n_qp)/(4 * N_0 * Delta)).to(1)
        self.Q_qp = ((2 * N_0 * Delta)/(self.alpha * gamma_s * S_1 * self.n_qp)).to(1)
        self.Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * self.T_b))/np.sinh(self.eta)/K_0(self.eta)
        #y = (self.omega_r*self.C_c*Z0).to(1)
        #Zin = Z0
        #Zout = 1./(1/Z0 + 1j*self.omega_r*C_load)
        #Gprime = (self.omega_r*self.C_c*Zin*y/(Zin + Zout)) - (1j*self.omega_r*self.C_c*Zin**2*y**2)/(Zin + Zout)**2
        #dQe = (Gprime/self.omega_r/self.C_i).to(1)
        #self.Q_e = self.Q_c*np.cos(self.phi_c)*np.exp(1j*self.phi_c)
        self.Q_c = (2 * self.C_i/(self.omega_r * self.C_c**2 * Z0)).to(1)
        #self.Q_c = 1./np.real(dQe)
        #self.Q_e = 1./dQe
        #self.phi_c = np.arctan2(np.imag(self.Q_e), np.real(self.Q_e))
        #self.Q_c = (2 * self.C_i/(self.omega_r * self.C_c**2 * Z0)).to(1)
        self.Q_i = 1./(1/self.Q_qp + 1./Q_int)
        self.Q_r  = 1./(1./self.Q_c + 1./self.Q_i)
        self.P_crit = (0.8 * (self.omega_r)* E_crit * self.Q_c/2/self.Q_r**3).to(u.pW)

        if verbose:
            niceprint ("")
            niceprint ("Q factor from qp losses, Q_qp ", self.Q_qp)
            niceprint ("Resonant Frequency, f_r ", self.f_r.to(u.MHz))
            niceprint ("Internal Q factor, Q_i ", self.Q_i)
            niceprint ("Coupling Q factor, Q_c ", self.Q_c)
            #niceprint ("Coupling phase angle, phi_c ", self.phi_c)
            niceprint ("Resonator Q factor, Q_r ", self.Q_r)
            niceprint ("Kinetic Inductance Fraction, alpha ", self.alpha)
            niceprint ("Beta ", self.beta)
            niceprint ("surface resistance, Rs", Rs)

        print (self.f_g, self.f_r)
        self.dx = ((self.f_g - self.f_r)/self.f_r).to(1)
        self.S_21 = 1 - (self.Q_r/self.Q_c)* 1./(1 + 2 * 1j * self.Q_r * self.dx)
        self.df_r = (self.f_r/(2*self.Q_r)).to('kHz')
        self.df_g = (self.f_g - self.f_r).to('kHz')

        self.chi_c = (4 * self.Q_r**2)/(self.Q_i * self.Q_c)
        self.chi_g = 1./(1 + (self.df_g/self.df_r)**2)
        self.chi_qp = self.Q_i/self.Q_qp
        self.P_diss = (self.chi_g*self.chi_c/2) * self.P_g
        self.E_stored = (self.P_diss*self.Q_i/self.omega_r).to(u.eV)
        self.N_ph = (self.E_stored/(hbar*self.omega_r)).to(1) # number of #microwave photons
        self.a_swenson = (2*self.Q_r**3/self.Q_c*self.P_g/self.omega_r/E_crit).to(1)

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

        self.a_phi = self.beta*self.kappa
        self.a_amp = self.kappa
        self.b_phi = 0
        self.b_amp = 0

        if verbose:
            niceprint ("")
            niceprint ("dln n_qp/ d ln T_b", self.kappa)
            niceprint ("Conductance of Island ", self.G_b)
            niceprint ("Heat Capacity of Island ", self.C_b)
            niceprint ("Bolometer Time constant", self.tau_b)
            niceprint ("Bolometer Roll off frequency", self.f_b)

            niceprint ("")
            niceprint ("The optical power is ", self.P_opt.to(u.pW))
            niceprint ("The readout power is ", self.P_g)
            niceprint ("The dissipated power is ", self.P_diss)
            niceprint ("The non-linearity parameter is ", self.a_swenson)
            niceprint ("E_crit/E_stored", (E_crit/self.E_stored).to(1))
            niceprint ("Critical readout power ", self.P_crit.to(dBm))
            niceprint ("The island power, P_b ", self.P_b)
            niceprint ("The Energy stored , E ", self.E_stored)
            print ("The Number of Microwave Photons, Nph %1.2e"%(self.N_ph))

        # Calculating the responsivities on resonance
        # r = dx/dPopt
        self.r = (self.xMB*self.kappa/((1 + self.xMB)*self.T_b*self.G_b)).to(1./u.pW)
        #self.r = (0.5 * (self.chi_qp * self.beta/self.Q_i) *\
        #    self.tau_qp/self.tau_th * self.kappa/self.P_b).to(1/u.pW)
        #pdb.set_trace()
        self.r_f = (self.f_r * self.r).to(u.kHz/u.pW)

        niceprint ("")
        niceprint ("resobolo responsivity ignoring bolometer rolloff", self.r_f)
        niceprint ("resobolo resonator rolloff", self.df_r)

    def plot_bolomat_noise(self, ax, fmin=-1, fmax=6):
        self.calculate_bolomat_noise(fmin, fmax)

        #ax.loglog(self.nu, resp, 'b', label='Phonon')

        ax.loglog(self.nu, self.bolomat_NEP_phonon, 'b', label='Phonon')
        ax.loglog(self.nu, self.bolomat_NEP_amp, 'r--', label='Amplifier')
        ax.loglog(self.nu, self.bolomat_NEP_jamp, 'g', label='Johnson Amplitude')
        ax.loglog(self.nu, self.bolomat_NEP_jph, 'g:', label='Johnson Phase')
        ax.loglog(self.nu, self.bolomat_NEP_phi, 'k', label='Total')

    def plot_bolomat_NEF(self, ax, fmin=-1, fmax=6):
        self.calculate_bolomat_noise(fmin, fmax)

        #ax.loglog(self.nu, resp, 'b', label='Phonon')

        ax.loglog(self.nu, self.bolomat_NEF_phonon, 'b', label='Phonon')
        ax.loglog(self.nu, self.bolomat_NEF_amp, 'r--', label='Amplifier')
        ax.loglog(self.nu, self.bolomat_NEF_jamp, 'g', label='Johnson Amplitude')
        ax.loglog(self.nu, self.bolomat_NEF_jph, 'g:', label='Johnson Phase')
        ax.loglog(self.nu, self.bolomat_NEF_phi, 'k', label='Total')

    def calculate_bolomat_noise(self, fmin=-1, fmax=6):

        self.nu = np.logspace(fmin, fmax, 5000) * u.Hz
        self.bolomat = np.zeros((self.nu.size,3,3), dtype=np.complex64)
        om = 2*pi*self.nu
        self.ones = np.ones_like(self.nu.value)
        self.bolomat[:,0,0] = (1j*2*om*self.Q_i/self.omega_g + self.Q_i/self.Q_r +\
            self.b_amp + 4*self.Q_i*self.Q_r*self.dx**2).to(1)
        self.bolomat[:,0,1] = (-1j*4*om*self.Q_i*self.Q_r*self.dx/self.omega_g -\
                2*self.Q_r*self.dx*self.b_amp).to(1)
        self.bolomat[:,0,2] = (2*self.a_amp/self.T_b).to(1./u.K).value
        self.bolomat[:,1,0] = (1j*4*om*self.Q_i*self.Q_r*self.dx/self.omega_g -\
                self.b_phi).to(1)
        self.bolomat[:,1,1] = (1j*2*om*self.Q_i/self.omega_g + self.Q_i/self.Q_r +\
                4*self.Q_i*self.Q_r*self.dx**2 +\
                2*self.Q_r*self.dx*self.b_phi).to(1)
        self.bolomat[:,1,2] = (-2*self.a_phi/self.T_b).to(1/u.K).value
        self.bolomat[:,2,0] = (-(1 + self.b_amp/2)*self.P_diss).to(u.W).value
        self.bolomat[:,2,1] = (self.dx*self.Q_r*(self.b_amp +\
                2)*self.P_diss).to(u.W).value
        self.bolomat[:,2,2] = (1j*om*self.C_b + self.G_b -\
                self.P_diss*self.a_amp/self.T_b).to(u.W/u.K).value
        self.invbolomat = np.linalg.inv(self.bolomat)

        #self.invbolomat[:,0,2] *= 1./u.pW
        #self.invbolomat[:,1,2] *= 1./u.pW
        #self.invbolomat[:,2,0] *= u.K
        #self.invbolomat[:,2,1] *= u.K
        #self.invbolomat[:,2,2] *= u.K/u.pW

        n = self.gamma_leg - 1
        chi_ph = (n+1)/(2*n+3) * ((self.T_0/self.T_b).to(1)**(2*n+3)-1)/\
                ((self.T_0/self.T_b).to(1)**(n+1) - 1)
        S_ph = (4 * chi_ph * k_B * self.T_b**2 * self.G_b ).to(u.W**2/u.Hz)
        A_ph = np.sqrt(S_ph).to(u.W/u.Hz**0.5).value
        self.N_ph = np.array([0,0,A_ph])

        self.y_ph = np.abs(np.dot(self.invbolomat, self.N_ph))[:, 1]/u.Hz**0.5
        #print (self.invbolomat[0, :, :])
        #exit()

        Me = np.sqrt(Delta/k_B/self.T_b).to(1)
        self.N_gamp = Me*np.array([
                np.sqrt(16*k_B*self.T_b/self.P_diss).to(1/u.Hz**0.5).value,
                0,
                -np.sqrt(4*k_B*self.T_b*self.P_diss).to(u.W/u.Hz**0.5).value])
        self.y_gamp = np.abs(np.dot(self.invbolomat, self.N_gamp))[:, 1]/u.Hz**0.5

        self.N_gphi = Me*np.array([0,
                1j*np.sqrt(16*k_B*self.T_b/self.P_diss).to(1/u.Hz**0.5).value,
                0])
        self.y_gphi = np.abs(np.dot(self.invbolomat, self.N_gphi))[:, 1]/u.Hz**0.5

        self.y_amp = np.sqrt(8*k_B*T_amp/self.P_diss).to(1/u.Hz**0.5)*self.ones


        # RMS phase noise
        self.y_phi_rms = np.sqrt(self.y_ph**2 + self.y_gamp**2 + self.y_gphi**2 + self.y_amp**2)

        #resp = 1./u.W
        self.bolomat_resp = np.abs(self.invbolomat[:,1,2]/u.W)
        self.bolomat_freqresp = (self.f_r/4/self.Q_r).to(u.Hz)


        self.bolomat_NEF_phi = np.abs(self.y_phi_rms*self.bolomat_freqresp).to(u.Hz/u.Hz**0.5).value
        self.bolomat_NEF_phonon = np.abs(self.y_ph*self.bolomat_freqresp).to(u.Hz/u.Hz**0.5).value
        self.bolomat_NEF_amp = np.abs(self.y_amp*self.bolomat_freqresp).to(u.Hz/u.Hz**0.5).value
        self.bolomat_NEF_jamp = np.abs(self.y_gamp*self.bolomat_freqresp).to(u.Hz/u.Hz**0.5).value
        self.bolomat_NEF_jph = np.abs(self.y_gphi*self.bolomat_freqresp).to(u.Hz/u.Hz**0.5).value

        self.bolomat_NEP_phi = np.abs(self.y_phi_rms/self.bolomat_resp).to(u.aW/u.Hz**0.5).value
        self.bolomat_NEP_phonon = np.abs(self.y_ph/self.bolomat_resp).to(u.aW/u.Hz**0.5).value
        self.bolomat_NEP_amp = np.abs(self.y_amp/self.bolomat_resp).to(u.aW/u.Hz**0.5).value
        self.bolomat_NEP_jamp = np.abs(self.y_gamp/self.bolomat_resp).to(u.aW/u.Hz**0.5).value
        self.bolomat_NEP_jph = np.abs(self.y_gphi/self.bolomat_resp).to(u.aW/u.Hz**0.5).value


    def calculate_noise_spectra(self):

        # Phonon NEP
        n = self.gamma_leg - 1
        chi_ph = (n+1)/(2*n+3) * ((self.T_0/self.T_b).to(1)**(2*n+3)-1)/\
                ((self.T_0/self.T_b).to(1)**(n+1) - 1)
        niceprint ("The bolometer flink factor", chi_ph)
        self.S_ph = (4 * chi_ph * k_B * self.T_b**2 * self.G_b ).to(u.aW**2/u.Hz)

        self.NEP_ph = self.S_ph ** 0.5
        self.NEF_ph = (self.r_f * self.NEP_ph).to(u.Hz**0.5)

        # Amplifier NEP
        self.S_amp = (self.Q_c/2/self.Q_r**2)**2*(k_B * T_amp/self.P_g).to(1/u.Hz)

        self.NEP_amp = (self.S_amp**0.5/self.r).to(u.aW/u.Hz**0.5)
        self.NEF_amp = (self.r_f * self.NEP_amp).to(u.Hz**0.5)
        #self.S_amp = (k_B * T_amp/self.P_g).to(1/u.Hz)

        #self.NEP_amp = (2 * self.S_amp**0.5/self.r/self.chi_c/self.Q_i).to(u.aW/u.Hz**0.5)

        # Shot NEP
        self.S_shot = 4 * self.n_qp * self.tau_qp / V_sc

        self.NEP_gr = (self.S_shot**0.5 * self.P_b/(self.n_qp * self.kappa)).to(u.aW/u.Hz**0.5)
        self.NEF_gr = (self.r_f * self.NEP_gr).to(u.Hz**0.5)
        #self.S_shot = self.n_qp * V_sc * (1/tau_max + 1/self.tau_qp)

        #self.NEP_gr = (2*self.S_shot**0.5 * self.P_b/(self.n_qp * self.kappa *
        #    V_sc) * self.tau_th).to(u.aW/u.Hz**0.5)
#        self.S_shot = 2 * self.n_qp * V_sc * (1/tau_max + 1/self.tau_qp)
#
#        self.NEP_gr = (self.S_shot**0.5 * (self.tau_th/(self.n_qp * V_sc) * self.P_b/self.kappa)).to(u.aW/u.Hz**0.5)

        if self.is_optical:
            # Photon NEP
            self.NEP_photon = np.sqrt(2*h*self.nu0*self.P_opt +
                    2*self.P_opt**2/(self.opt_bw*self.nu0)).to(u.aW/u.Hz**0.5)
            self.NEF_photon = (self.r_f * self.NEP_photon).to(u.Hz**0.5)

        # TLS NEP
        self.alpha_tls = 0.5
        self.beta_tls = 2
        # Made this prediction off of Jonas's review paper
        self.kappatls0 = 1.9e-19/u.Hz*u.Kelvin**self.beta_tls*u.Hz**0.5
        #self.kappatls0 = 1.9e-17/u.Hz*u.Kelvin**self.beta_tls*u.Hz**0.5
        self.nu_tls = 1 * u.Hz
        self.nu_ref = 1 * u.kHz
        # TLS spectrum at 1 Hz
        self.Stls = self.kappatls0/np.sqrt(self.N_ph)*self.T_b**(-self.beta_tls)* self.nu_tls**(-self.alpha_tls)
        self.NEF_tls = (self.Stls**0.5*self.f_r).to(u.Hz**0.5)
        self.NEP_tls = (self.Stls**0.5/self.r).to(u.aW/u.Hz**0.5)

        if self.is_optical:
            # Total NEP
            self.NEP_total = (self.NEP_photon**2 + self.NEP_gr**2 +
                    self.NEP_amp**2 + self.NEP_ph**2 + self.NEP_tls**2)**0.5
            self.NEF_total = (self.NEF_photon**2 + self.NEF_gr**2 +
                    self.NEF_amp**2 + self.NEF_ph**2 + self.NEF_tls**2)**0.5
        else:
            # Total NEP
            self.NEP_total = (self.NEP_gr**2 + self.NEP_amp**2 + self.NEP_ph**2 + self.NEP_tls**2)**0.5
            self.NEF_total = (self.NEF_gr**2 + self.NEF_amp**2 + self.NEF_ph**2 + self.NEF_tls**2)**0.5

        niceprint ("")
        if self.is_optical:
            niceprint ("Photon NEF", self.NEF_photon)
        niceprint ("Phonon NEF", self.NEF_ph)
        niceprint ("Amplifier NEF", self.NEF_amp)
        niceprint ("Shot NEF", self.NEF_gr)
        niceprint ("TLS NEF @ 1Hz", self.NEF_tls)
        niceprint ("Total NEF", self.NEF_total)
        niceprint ("")
        if self.is_optical:
            niceprint ("Photon NEP", self.NEP_photon)
        niceprint ("Phonon NEP", self.NEP_ph)
        niceprint ("Amplifier NEP", self.NEP_amp)
        niceprint ("Shot NEP", self.NEP_gr)
        niceprint ("TLS NEP @ 1Hz", self.NEP_tls)
        niceprint ("Total NEP", self.NEP_total)


    def plot_NEF(self, ax, fmin=-1, fmax=6):
        self.calculate_noise_spectra()
        self.nu = np.logspace(fmin, fmax, 5000) * u.Hz
        self.ones = np.ones_like(self.nu.value)
        self.bolo_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_b).to(1))
        self.qp_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_qp).to(1))
        self.reso_rolloff =  1/(1 + 1j * (self.nu / self.df_r).to(1))
        self.bolo_rolloff = np.sqrt(np.abs(self.bolo_rolloff))
        self.qp_rolloff = np.sqrt(np.abs(self.qp_rolloff))
        self.reso_rolloff = np.sqrt(np.abs(self.reso_rolloff))

        if self.is_optical:
            self.NEF_photon_spec = self.NEF_photon * self.bolo_rolloff
        self.NEF_ph_spec = self.NEF_ph * self.bolo_rolloff
        self.NEF_amp_spec = self.NEF_amp * self.ones
        self.NEF_gr_spec = self.NEF_gr * self.qp_rolloff * self.reso_rolloff
        self.NEF_tls_spec = self.NEF_tls * (self.nu/self.nu_tls)**(-0.25)
        self.NEF_total_spec = (self.NEF_gr_spec**2 + self.NEF_amp_spec**2 +
                self.NEF_ph_spec**2 + self.NEF_tls_spec**2)**0.5

        if self.is_optical:
            ax.loglog(self.nu, self.NEF_photon_spec, 'b-.', label='Photon')
        ax.loglog(self.nu, self.NEF_ph_spec, 'b', label='Phonon')
        ax.loglog(self.nu, self.NEF_amp_spec, 'r--', label='Amplifier')
        ax.loglog(self.nu, self.NEF_gr_spec, 'g', label='Gen-Recomb')
        ax.loglog(self.nu, self.NEF_tls_spec, 'm-.', label='TLS')
        ax.loglog(self.nu, self.NEF_total_spec, 'k', label='Total')

    def plot_noise(self, ax, fmin=-1, fmax=6):
        self.calculate_noise_spectra()
        self.nu = np.logspace(fmin, fmax, 5000) * u.Hz
        self.ones = np.ones_like(self.nu.value)
        self.bolo_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_b).to(1))
        self.qp_rolloff =  1/(1 + 1j * 2 * np.pi * (self.nu * self.tau_qp).to(1))
        self.reso_rolloff =  1/(1 + 1j * (self.nu / self.df_r).to(1))
        self.bolo_rolloff = np.sqrt(np.abs(self.bolo_rolloff))
        self.qp_rolloff = np.sqrt(np.abs(self.qp_rolloff))
        self.reso_rolloff = np.sqrt(np.abs(self.reso_rolloff))

        if self.is_optical:
            self.NEP_photon_spec = self.NEP_photon * self.bolo_rolloff
        self.NEP_ph_spec = self.NEP_ph * self.bolo_rolloff
        self.NEP_amp_spec = self.NEP_amp * self.ones
        self.NEP_gr_spec = self.NEP_gr * self.qp_rolloff * self.reso_rolloff
        self.NEP_tls_spec = self.NEP_tls * (self.nu/self.nu_tls)**(-0.25)
        self.NEP_total_spec = (self.NEP_gr_spec**2 + self.NEP_amp_spec**2 +
                self.NEP_ph_spec**2 + self.NEP_tls_spec**2)**0.5

        if self.is_optical:
            ax.loglog(self.nu, self.NEP_photon_spec, 'b-.', label='Photon')
        ax.loglog(self.nu, self.NEP_ph_spec, 'b', label='Phonon')
        ax.loglog(self.nu, self.NEP_amp_spec, 'r--', label='Amplifier')
        ax.loglog(self.nu, self.NEP_gr_spec, 'g', label='Gen-Recomb')
        ax.loglog(self.nu, self.NEP_tls_spec, 'm-.', label='TLS')
        ax.loglog(self.nu, self.NEP_total_spec, 'k', label='Total')


if __name__=="__main__":
    # We are using the heater pad to mimic the optical load on the bolometer
    P_opt = np.linspace(0, 20, 21)*u.pW + 1e-3*u.pW
    Npowers = P_opt.size

    Rh = 0.150 * u.Ohm
    frequencies = np.array([305.8, 318.4, 337.4])*u.MHz
    Rb = np.array([300, 404, 510]) * u.kOhm
    Pref = 10e-12
    gamma_leg = np.array([2.962, 2.754, 2.862])# conductivity index = beta + 1
    K_leg = np.array([352, 165, 122])* u.picoWatt
    Cc = np.array([0.3984, 0.3731, 0.3478])*u.pF

    #frequencies = np.arange(400, 800, 50)*u.MHz
    #gamma_leg = np.ones_like(frequencies)*2.862# conductivity index = beta + 1
    #K_leg = np.ones_like(frequencies) * 122 * u.picoWatt
    #Cc = np.ones_like(frequencies) * 0.3478 * u.pF

    #Vmax = (Pref/Rh.value)**0.5*510e3
    #Vdc = np.linspace(0, Vmax, Npowers)*u.Volt
    #P_opt = ((Rh/Rb[:, np.newaxis]**2)*Vdc**2).to(u.pW)
    #P_opt = 4.7 * u.pW
    #print (P_opt)
    T_c = 1.284 * u.Kelvin
    T_0 = 0.25 * u.Kelvin
    P_g = -100*dBm
    #P_g -= 20*np.log10(4)*dB # Accounting for the 4 resonators being read at the
    # same time
    #P_opt = 10*u.pW
    NEPs = np.zeros((len(frequencies), Npowers))
    NEFs = np.zeros((len(frequencies), Npowers))
    NEPs_ph = np.zeros((len(frequencies), Npowers))
    NEPs_amp = np.zeros((len(frequencies), Npowers))
    for ireso, reso in enumerate(frequencies):
        op = OperatingPoint(P_opt=P_opt, f_r=reso, f_g=reso, T_0=T_0, P_g=P_g)
        op.C_c1 = Cc[ireso]
        op.gamma_leg = gamma_leg[ireso]
        op.K_leg = K_leg[ireso]/u.Kelvin**gamma_leg[ireso]
        for ipow in range(Npowers):
            #if ipow > 0: continue
            op.P_opt = P_opt[ipow]
            op.is_optical = False
            op.calculate_noise()
            op.calculate_noise_spectra()
            fig, ax = plt.subplots(figsize=(10,10))
            if useBoloMat:
                op.plot_bolomat_noise(ax,fmin=-1, fmax=6)
                NEPs[ireso, ipow] = op.bolomat_NEP_phi[0]
                NEFs[ireso, ipow] = op.bolomat_NEF_phi[0]
                NEPs_ph[ireso, ipow] = op.bolomat_NEP_phonon[0]
                NEPs_amp[ireso, ipow] = op.bolomat_NEP_amp[0]
            else:
                op.plot_noise(ax,fmin=-1, fmax=6)
                NEPs[ireso, ipow] = op.NEP_total.value
                NEFs[ireso, ipow] = op.NEF_total.value
                NEPs_ph[ireso, ipow] = op.NEP_ph.value
                NEPs_amp[ireso, ipow] = op.NEP_amp.value
            ax.set_title("%.1f MHz reso biased at %.1f pW"%(frequencies[ireso].value, P_opt[ipow].value))
            ax.set_xlabel(r'Frequency [Hz]')
            ax.set_ylabel(r'NEP [aW/rtHz]')
            ax.set_ylim([1, 1000])
            ax.set_xlim([0.1, 1e6])
            ax.grid(which='major', axis='both')
            ax.legend(loc='best')
            plt.savefig(plot_dir + 'Waffle_TKID_%.1fMHz_%.1fpW_NEP.png'%(frequencies[ireso].value, P_opt[ipow].value))
            plt.close()
            fig, ax = plt.subplots(figsize=(10,10))
            if useBoloMat:
                op.plot_bolomat_NEF(ax,fmin=-1, fmax=6)
            else:
                op.plot_NEF(ax,fmin=-1, fmax=6)
            ax.set_title("%.1f MHz reso biased at %.1f pW"%(frequencies[ireso].value, P_opt[ipow].value))
            ax.set_xlabel(r'Frequency [Hz]')
            ax.set_ylabel(r'NEF [Hz/rtHz]')
            ax.set_ylim([0.01, 100])
            ax.set_xlim([0.1, 1e6])
            ax.grid(which='major', axis='both')
            ax.legend(loc='best')
            plt.savefig(plot_dir + 'Waffle_TKID_%.1fMHz_%.1fpW_NEF.png'%(frequencies[ireso].value, P_opt[ipow].value))
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
    P_readout = (1*u.mW*10**(P_g.value/10)).to(u.pW)
    print (P_readout)
    for ireso, reso in enumerate(frequencies):
        fig, ax = plt.subplots(figsize=(10,10))
        ax.plot(P_opt.value, NEPs[ireso], 'bo', label='total', ls='None', ms=12)
        ax.plot(P_opt.value, NEPs_ph[ireso], 'kd', label='phonon', ls='None', ms=12)
        ax.plot(P_opt.value, NEPs_amp[ireso],'rs', label='amplifier', ls='None', ms=12)
        ax.set_title("%.1f MHz resonator"%(frequencies[ireso].value))
        ax.set_xlabel(r'Loading [pW]')
        ax.set_ylabel(r'NEP [aW/rtHz]')
        ax.legend(loc="upper left")
        #ax.set_ylim([1, 1000])
        #ax.set_xlim([0.1, 1e6])
        ax.grid(which='major', axis='both')
        plt.savefig(plot_dir + 'Waffle_TKID_%.1fMHz_NEP_vs_loading.png'%(frequencies[ireso].value))
        plt.close()

        fig, ax = plt.subplots(figsize=(10,10))
        ax.plot(P_opt.value, NEFs[ireso], 'bo', ls='None', ms=12)
        ax.set_title("%.1f MHz resonator"%(frequencies[ireso].value))
        ax.set_xlabel(r'Loading [pW]')
        ax.set_ylabel(r'NEF [Hz/rtHz]')
        #ax.set_ylim([1, 1000])
        #ax.set_xlim([0.1, 1e6])
        ax.grid(which='major', axis='both')
        plt.savefig(plot_dir + 'Waffle_TKID_%.1fMHz_NEF_vs_loading.png'%(frequencies[ireso].value))
        plt.close()
