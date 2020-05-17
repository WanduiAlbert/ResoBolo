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
Q_int = 3e7
# f_g = 205.128 * u.MHz
f_g = 305.8 * u.MHz
T_amp = 5.22 * u.Kelvin
eta_read = 0.0
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

gamma_leg = 2.130 # conductivity index = beta + 1
K_leg = 403.705 * u.picoWatt/u.Kelvin**gamma_leg
T_c = 1.329 * u.Kelvin
T_0 = 0.25 * u.Kelvin # Previously 0.23K Temperature of the thermal bath

# Material properties of the Aluminum superconductors
#tau_max = 955.1 * u.microsecond
tau_max = 488.1 * u.microsecond
n_qp_star = 763 * 1/u.um**3
gamma = 1.35 * u.mJ/u.mol/u.Kelvin**2
density = 2.7 * u.g/u.cm**3
A_r = 26.98 * u.g/u.mol
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')

rho = 1.15 * u.uOhm * u.cm
L_g = 6.1 * u.nH # From simulations
Z0 = 50 * u.Ohm # Characteristic impedance of the line
C_c1 = 0.3984 * u.pF * 2
C_c = C_c1/2.

class TKIDBolometer:

    def __init__(self, fr, f_g, Tc,  Cc, T_0,
            P_opt, P_read, K_leg, gamma_leg, L_g=6.1, **kwargs):
        self.C_b = 0.14 * u.pJ/u.Kelvin # From TKID paper
        self.f_r = fr * u.MHz
        self.f_g = f_g * u.MHz
        self.Tc = Tc * u.Kelvin
        self.L_g = L_g * u.nH
        self.Cc1 = Cc * u.pF
        self.C_c = self.Cc1/2
        self.T_0 = T_0 * u.Kelvin #Bath temperature
        self.P_opt = P_opt * u.pW
        self.P_read = (P_read * dBm).to(u.pW)
        self.P_total = self.P_opt + 0.5*self.P_read #kind of approximate
        self.gamma_leg = gamma_leg #n+1
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
        self.Delta = (1.764 * k_B * self.Tc).to('J')

    def set_operating_point(self):
        self.T_b= (((self.P_total)/self.K_leg +
            self.T_0**self.gamma_leg)**(1./self.gamma_leg)).to('K')
        niceprint("The temperature of the island", self.T_b)
        self.l = self.N_sq * self.w_trace
        self.A = self.t * self.w_trace # cross-sectional area
        self.V_sc = self.l * self.A

        self.Rs = (rho/self.t).to('Ohm') # surface resistance in ohms/sq
        niceprint ("The surface resistance", self.Rs)
        self.L_k_psq = (h * self.Rs/(2 * np.pi**2 * self.Delta)).to(u.pH)
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
        self.f_qp = (1./2/np.pi/self.tau_qp).to(u.kHz) # qp roll off frequency
        self.x_MB = - (self.alpha * self.S_2 * self.n_qp / (4 * N_0 *
                self.Delta)).to(1)
        self.Q_qp = ((2 * N_0 * self.Delta)/(self.alpha * gamma_s * self.S_1 * self.n_qp)).to(1)
        self.Q_sigma = (np.pi/4)*np.exp(self.Delta/(k_B * self.T_b))/np.sinh(self.eta)/K_0(self.eta)
        self.Q_i = 1./(1/self.Q_qp + 1./Q_int)
        self.Q_c = self.Q_i #(2 * self.C_i/(self.omega_r * self.C_c**2 * Z0)).to(1)
        self.Q_r  = 1./(1./self.Q_c + 1./self.Q_i)
        self.f_qr = (1./2/np.pi/self.tau_qp).to(u.Hz) # qp roll off frequency
        # Overwrite to simulate more realistic devices
        #Q_i = 48494
        #Q_c = 155298
        #Q_r = 1./(1./Q_c + 1./Q_i)
        self.P_crit = (0.8 * (2*np.pi*self.f_r)* self.E_crit *
                self.Q_c/2/self.Q_r**3).to(u.pW)

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
        self.S_21 = 1 - (self.Q_r/self.Q_c)* 1./(1 + 2j * self.Q_r * self.dx)
        self.df_r = (self.f_r/(2*self.Q_r)).to('kHz')
        self.df_g = (self.f_g - self.f_r).to('kHz')

        self.chi_c = (4 * self.Q_r**2)/(self.Q_i * self.Q_c)
        self.chi_g = 1./(1 + (self.df_g/self.df_r)**2)
        self.chi_qp = self.Q_i/self.Q_qp
        self.P_g = (0.5 * self.chi_c * self.chi_g) * self.P_read

        niceprint("")
        niceprint ("Coupling efficiency", self.chi_c)
        niceprint ("Detuning efficiency", self.chi_g)
        niceprint ("Fraction of Q_i from qp losses", self.chi_qp)

        # Now we include the NEP estimates for the resobolo currently
        self.kappa = (1/2 + self.Delta/k_B/self.T_b).to(1)
        self.P_leg = self.P_total # Total power into the resobolo thermal link
        self.gamma_g = (self.K_leg * self.gamma_leg *
                self.T_b**self.gamma_leg / self.P_leg).to(1)
        self.G_b = (self.K_leg*self.gamma_leg*self.T_b**(self.gamma_leg-1)).to('pW/K') # Conductance of the resobolo
        self.P_b = self.G_b * self.T_b
        self.tau_b = (self.C_b/self.G_b).to('us')
        self.f_b = (1./2/np.pi/self.tau_b).to(u.Hz) # roll off frequency for bolometer

        niceprint ("")
        niceprint ("dln n_qp/ d ln T_b", self.kappa)
        niceprint ("Conductance of Island ", self.G_b)
        niceprint ("Heat Capacity of Island ", self.C_b)
        niceprint ("Quasiparticle lifetime, tau_qp ", self.tau_qp)
        niceprint ("Equilibrium qp lifetime, tau_th ", self.tau_th)
        niceprint ("Bolometer Time constant", self.tau_b)
        niceprint ("Resonator Bandwidth", self.df_r)
        niceprint ("QP Roll off frequency", self.f_qp)
        niceprint ("Bolometer Roll off frequency", self.f_b)

        niceprint ("")
        niceprint ("The optical power is ", self.P_opt)
        niceprint ("The readout power is ", self.P_read)
        niceprint ("The dissipated power is ", self.P_g)
        niceprint ("Critical readout power ", self.P_crit)
        niceprint ("The island power, P_b ", self.P_b)

    def get_responsivity(self):
        self.r = (np.abs(self.x_MB) * self.kappa / (self.T_b * self.G_b)).to(1/u.pW)
        self.r_f = (self.f_r * self.r).to(u.kHz/u.pW)
        niceprint ("")
        print ("resobolo responsivity ignoring bolometer rolloff", self.r)
        niceprint ("resobolo responsivity ignoring bolometer rolloff", self.r_f)
        return self.r_f

    def get_noise(self):
        self.get_responsivity()
        # Phonon self.NEP
        n = self.gamma_leg - 1
        chi_ph = (n+1)/(2*n+3) * ((self.T_0/self.T_b).to(1)**(2*n+3)-1)/\
                ((self.T_0/self.T_b).to(1)**(n+1) - 1)
        niceprint ("The bolometer flink factor ", chi_ph)
        self.S_ph = (4 * chi_ph * k_B * self.T_b**2 * self.G_b ).to(u.aW**2/u.Hz)
        self.NEP_ph = self.S_ph ** 0.5

        # Amplifier self.NEP
        self.S_amp = (k_B * T_amp/self.P_read).to(1/u.Hz)
        # self.NEP_amp = (2 * S_amp**0.5/r/chi_c/Q_i).to(u.aW/u.Hz**0.5)
        #self.NEP_amp = (2 * self.S_amp**0.5 * self.Q_c/self.Q_r**2/self.r).to(u.aW/u.Hz**0.5)
        self.NEP_amp = (self.S_amp**0.5 * self.Q_c/2/self.Q_r**2/self.r).to(u.aW/u.Hz**0.5)

        # Shot self.NEP
        #self.S_shot = 2 * self.n_qp * self.V_sc * (1/tau_max + 1/self.tau_qp)
        #self.NEP_gr = (self.S_shot**0.5 * (self.tau_th/(self.n_qp * self.V_sc) *
        #    self.P_b/self.kappa)).to(u.aW/u.Hz**0.5)
        self.NEP_gr = (2*(self.tau_qp/(self.n_qp * self.V_sc))**0.5 *
            self.P_b/self.kappa).to(u.aW/u.Hz**0.5)

        # Photon NEP
        self.NEP_photon = np.sqrt(2*self.P_opt*h*nu_opt*dnu_opt +
                2*self.P_opt**2/(dnu_opt*nu_opt)).to(u.aW/u.Hz**0.5)


        # Total self.NEP
        self.NEP_total = (self.NEP_gr**2 + self.NEP_amp**2 + self.NEP_ph**2)**0.5

        niceprint ("")
        niceprint ("Phonon NEP", self.NEP_ph)
        niceprint ("Amplifier NEP", self.NEP_amp)
        niceprint ("GR NEP", self.NEP_gr)
        niceprint ("Photon NEP", self.NEP_photon)
        niceprint ("Total NEP", self.NEP_total)


        return self.NEP_ph, self.NEP_amp, self.NEP_gr, self.NEP_total

    def get_noiseindBcHz(self):
        responsivity = self.r
        nep_dark = self.NEP_total
        nep_light = self.NEP_photon
        #nep = np.sqrt(nep_dark**2 + nep_light**2)
        nep = nep_dark
        chifactor = 2*self.Q_r**2/self.Q_c
        noise = ((chifactor*responsivity*nep)**2).to(1./u.Hz)
        self.perHznoise = noise
        self.dBcperHznoise = 10*np.log10(noise.value)/u.Hz
        self.ampdBcperHznoise = 10*np.log10(self.S_amp.value)/u.Hz
        niceprint ("")
        print ("Noise in 1/Hz", self.perHznoise)
        print ("Amplifier Noise in 1/Hz", self.S_amp)
        print ("Amplifier Noise in dBc/Hz", self.ampdBcperHznoise)
        niceprint ("Noise in dBc/Hz", self.dBcperHznoise)
        return self.perHznoise, self.dBcperHznoise

    def get_noise_spectra(self):
        self.get_noise()
        self.nu = np.logspace(-1, 6, 5000) * u.Hz
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

        fig, ax = plt.subplots(figsize=(10,10))
        ax.loglog(self.nu, self.NEP_ph_spec, 'b', label='Phonon')
        ax.loglog(self.nu, self.NEP_amp_spec, 'r--', label='Amplifier')
        ax.loglog(self.nu, self.NEP_gr_spec, 'g', label='Gen-Recomb')
        ax.loglog(self.nu, self.NEP_total_spec, 'k', label='Total')

        ax.set_xlabel(r'Frequency [Hz]')
        ax.set_ylabel(r'NEP [aW/rtHz]')
        ax.set_ylim([1, 1000])
        ax.set_xlim([0.1, 1e6])
        ax.grid(which='major', axis='both')
        ax.legend(loc='best')

        plt.savefig(r'NEP_spectra.png')


if __name__=="__main__":
    f0 = 400
    df = 2.3
    Nreso = 128
    frs = f0 + df*np.arange(Nreso)
    P_opt = 5.00 #pW
    L_g = 6.1
    C_c1 = 0.2878
    P_read = -90 #dBm
    #P_opt = 10.80 #pW
    #P_opt = 13.90 #pW

    gamma_leg = 3#2.975
    K_leg = 122.9
    #K_leg = 164.9
    #K_leg = 352.5


    T_c = 1.278
    T_0 = 0.25

    noiseperhz = np.zeros(Nreso)
    dBcnoiseperhz = np.zeros(Nreso)
    NEP_phs = np.zeros(Nreso)
    NEP_amps = np.zeros(Nreso)
    NEP_grs = np.zeros(Nreso)
    NEP_totals = np.zeros(Nreso)

    for ireso in range(Nreso):
        fg = fr = frs[ireso]
        print ("Working on resonator ", fr)
        fr = fg
        bolo = TKIDBolometer(fr, fg, T_c, C_c1, T_0,
                P_opt, P_read, K_leg, gamma_leg)
        bolo.set_operating_point()
        NEP_ph, NEP_amp, NEP_gr, NEP_tot = bolo.get_noise()
        perhz, dbcperhz = bolo.get_noiseindBcHz()
        NEP_phs[ireso] = NEP_ph.value
        NEP_amps[ireso] = NEP_amp.value
        NEP_grs[ireso] = NEP_gr.value
        NEP_totals[ireso] = NEP_tot.value
        NEP_phs[ireso] = NEP_ph.value
        noiseperhz[ireso] = perhz.value
        dBcnoiseperhz[ireso] = dbcperhz.value
    #bolo.get_noise_spectra()

    plt.figure()
    plt.plot(frs, dBcnoiseperhz)
    plt.xlabel('Resonator Frequency [MHz]')
    plt.ylabel('Noise [dBc/Hz]')
    plt.title('TKID Dynamic Range with No Sky Loading')
    plt.grid()
    plt.savefig('tkid_dynamicrange_vs_freq_dark.png')
    plt.show()

    plt.figure()
    plt.plot(frs, NEP_phs, 'b', label='Phonon')
    plt.plot(frs, NEP_amps, 'r--', label='Amplifier')
    plt.plot(frs, NEP_grs, 'g', label='Gen-Recomb')
    plt.plot(frs, NEP_totals, 'k', label='Total')
    plt.xlabel('Resonator Frequency [MHz]')
    plt.ylabel(r'NEP [aW/rtHz]')
    t.legend(loc='best')
    plt.grid()
    plt.savefig('tkid_nep_vs_freq.png')
    pHlt.show()



