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
T_amp = 5.22 * Kelvin
eta_read = 0.01
chi_ph = 0.658 # flink factor for phonon noise
# Material properties of the Aluminum superconductor
tau_max = 0.1809 * microsecond
n_qp_star = 518 * 1/um**3
R = 1./(n_qp_star * tau_max)
# Physical properties of the superconductor + resonator capacitor
t = 0.05 * um
w_trace = 1 * um #width of inductor trace
s_trace = 1 * um #spacing between inductor meanders
N_sq = 16200  # Number of squares
l = N_sq * w_trace # Total length of the inductor (approximately. More exact is +21um + ...)
A = t * w_trace # cross-sectional area
V_sc = l * A

P_opt = 14.00*pW #((Vdc/Rb)**2 * Rh).to(u.pW)

n = 2.000 # conductivity index = beta + 1
band = 270
if band == 150:
    K_leg = 122.9 * pW/Kelvin**(n+1)  #150 GHz case
elif band == 220:
    K_leg = 275.2 * pW/Kelvin**(n+1)   #220 GHz case
elif band == 270:
    K_leg = 352.5 * pW/Kelvin**(n+1)   #270 GHz case

T_c = 1.278 * Kelvin
Tstart = 0.2
T_0 = 0.23 * Kelvin # Previously 0.23K Temperature of the thermal bath
alphak = 0.378

Pg_dBm = -90
Pg = 1e-3 * 10**(Pg_dBm/10)

gamma = 1.35 * mJ/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g/mol
rho = 1.255 * uOhm * cm

N_0 = (3 * (gamma * (density/A_r)))/(2*pi**2 * k_B**2)


#plt.figure(321, figsize=(10,10))
#plt.figure(213, figsize=(10,10))
Delta = 1.764 * k_B * T_c

varyQc = False
varyfreq = True
varyTemp = False
analyzeTempandNumber = False

def get_islandtemperature(Tbath, P):
    return (P/K + Tbath**(n+1))**(1./(n+1))

def get_islandpower(T, Tbath):
    return K*(T**(n+1) - Tbath**(n+1))

def get_rjtemp(P, nu0, bw):
    return P/(k*nu0*bw)

def get_beta(T, nu):

    delta = 3.5*k*Tc/2
    q = (h*nu/2/k/T)
    nqp = 2*N0*np.sqrt(2*pi*delta*k*T)*np.exp(-delta/(k*T))
    S1 = 2/pi * np.sqrt(2*delta/(pi*k*T))*np.sinh(q)*K0(q)
    S2 = 1 + np.sqrt(2*delta/(pi*k*T))*np.exp(-q)*I0(q)
    return S2/S1

def get_xQrMB(T, f):
    eta = h * f / (2*k_B * T)
    S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k_B*T))*np.sinh(eta)*K_0(eta)
    S_2 = 1 + np.sqrt(2*Delta/(pi*k_B*T)) * np.exp(-eta) * I_0(eta)

    n_th = 2*N_0 * np.sqrt(2*pi* k_B * T* Delta)*np.exp(-Delta/(k_B*T))
    n_qp = n_th
    Q_qp = (2 * N_0 * Delta)/(alphak * S_1 * n_qp)
    Q_i = 1./(1/Q_qp + 1./Q_int)
    Q_c = Q_i
    Q_r  = 1./(1./Q_c + 1./Q_i)
    x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)
    return x, Q_r

def get_xQMB(T, f):

    eta = h * f / (2*k_B * T)
    S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k_B*T))*np.sinh(eta)*K_0(eta)
    S_2 = 1 + np.sqrt(2*Delta/(pi*k_B*T)) * np.exp(-eta) * I_0(eta)

    n_th = 2*N_0 * np.sqrt(2*pi* k_B * T* Delta)*np.exp(-Delta/(k_B*T))
    n_qp = n_th
    Q_qp = (2 * N_0 * Delta)/(alphak * S_1 * n_qp)
    Q_i = 1./(1/Q_qp + 1./Q_int)
    x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)
    return x, 1./Q_i

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


if varyQc:
    if band == 150:
        f_r = 324*MHz
    elif band == 220:
        f_r = 475*MHz
    elif band == 270:
        f_r = 583*MHz
    f_g = f_r
    omega_r = 2*pi*f_r
    Qc = np.logspace(3,5,2000)
    T_op = 380e-3
    T = T_op
    kappa = 0.5 + Delta/(k_B*T)
    G = K_leg * (n+1)*T**n

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
    #Q_c = 22400
    Q_c = Qc
    Q_i = 1./(1/Q_qp + 1./Q_int)
    Q_r  = 1./(1./Q_c + 1./Q_i)
    chi_c = 4*Q_r**2/(Q_c*Q_i)
    print ("Qi ", Q_i)
    print ("Qc ", Q_c)
    print ("chi_c", chi_c)
    chi_g = 1
    x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)
    P_diss = (chi_g*chi_c/2) * Pg
    Pc_tls = 1e-3*10**(-95.3/10)
    Ec_tls = (Pc_tls*Q_i/omega_r)
    Nc_tls = (Ec_tls/(hbar*omega_r))
    E_stored = (P_diss*Q_i/omega_r)
    N_ph = (E_stored/(hbar*omega_r)) # number of #microwave photons
    #responsivity
    S = f_r * x * kappa / (G * T)

    NEP_ph = np.sqrt(4*chi_ph*k_B*T**2*G)
    #NEP_amp = (G*T/(kappa*x))*(2/Q_i)*np.sqrt(k_B*T_amp/Pg)
    NEP_amp = (f_r/S)*(Q_c/(2*Q_r**2))*np.sqrt(k_B*T_amp/Pg)
    NEP_gr = (2*G*T/n_qp/kappa)/np.sqrt(R*V_sc)
    # TLS NEP
    alpha_tls = 0.5
    beta_tls = 2
    kappatls0 = 3.2e-16/Hz*Kelvin**beta_tls*Hz**0.5/np.sqrt(Nc_tls)
    #print (kappatls0)
    nu_tls = 1 * Hz
    # TLS spectrum at 1 Hz
    Stls = kappatls0/np.sqrt(1 + N_ph/Nc_tls)*T**(-beta_tls)* nu_tls**(-alpha_tls)
    #print (Stls)
    NEP_tls = (Stls**0.5*f_r/S)
    ones = np.ones_like(Qc)
    if band == 150:
        N_photon = 43#aW/rtHz
    elif band == 220:
        N_photon = 82#aW/rtHz
    elif band == 270:
        N_photon = 92#aW/rtHz
    NEP_total = np.sqrt(NEP_ph**2 + NEP_amp**2 + NEP_gr**2 + NEP_tls**2)

    fig, ax = plt.subplots(figsize=(10,10))
    ax.semilogx(Qc, NEP_total/aW, label='Total Noise')
    ax.semilogx(Qc, NEP_ph/aW*ones, label='Phonon Noise')
    ax.semilogx(Qc, NEP_amp/aW, label='Amplifier Noise')
    ax.axvline(Q_i, color='k', ls='--', lw=2)
    ax.axhline(N_photon, color='k', ls='-.', lw=2)
    ax.grid()
    ax.legend(loc='upper left', title='NEP')
    ax.set_xlim(left=Qc[0], right=Qc[-1])
    ax.set_ylim(bottom=0, top=100)
    ax.set_xlabel('Qc')
    ax.set_title('%d GHz band, NEP(T=%d mK, f=1 Hz)'%(band, T_op*1e3))
    ax.set_ylabel('NEP [aW/$\sqrt{\mathrm{Hz}}$]')
    plt.savefig('Qc_dependence_of_NEP_%dGHzband.png'%(band))
    plt.show()
elif varyfreq:
    ximax = 0.009
    f_r = np.r_[300:1300:2000j]*MHz
    f_g = f_r
    omega_r = 2*pi*f_r
    Qc = 10000
    T_op = 380e-3
    T = T_op
    kappa = 0.5 + Delta/(k_B*T)
    G = K_leg * (n+1)*T**n

    eta = h * f_g / (2*k_B * T)
    S_1 = (2/pi)*np.sqrt(2*Delta/(pi*k_B*T))*np.sinh(eta)*K_0(eta)
    S_2 = 1 + np.sqrt(2*Delta/(pi*k_B*T)) * np.exp(-eta) * I_0(eta)
    beta = S_2/S_1

    Gamma_gen = 0
    n_th = 2*N_0 * np.sqrt(2*pi* k_B * T* Delta)*np.exp(-Delta/(k_B*T))

    # Quality factors
    n_qp = np.sqrt((n_th + n_qp_star)**2 + (2*Gamma_gen*n_qp_star*tau_max/V_sc)) - n_qp_star
    tau_qp = tau_max/(1 + n_qp/n_qp_star)
    Q_qp = (2 * N_0 * Delta)/(alphak * S_1 * n_qp)
    Q_sigma = (np.pi/4)*np.exp(Delta/(k_B * T))/np.sinh(eta)/K_0(eta)
    #Q_c = 22400
    Q_c = Qc
    Q_i = 1./(1/Q_qp + 1./Q_int)
    Q_c = Q_i
    Q_r  = 1./(1./Q_c + 1./Q_i)
    chi_c = 4*Q_r**2/(Q_c*Q_i)
    print ("Qi ", Q_i)
    print ("Qc ", Q_c)
    print ("chi_c", chi_c)
    chi_g = 1
    x = (alphak * n_qp * S_2)/(4 * N_0 * Delta)

    geom_delta = x/beta/ximax**0.5
    arith_deltaf = f_r*x*(1-x)/beta/ximax**0.5
    fr_actual = f_r*(1 + x)
    resobw = fr_actual/(2*Q_r)
    targetspacing = resobw/ximax**0.5
    ymax = resobw.max()/kHz
    ymin = resobw.min()/kHz
    mask150 = np.ones_like(f_r, dtype=bool)
    mask150[f_r/MHz < 400] = False
    mask150[f_r/MHz > 2*400] = False
    mask220 = np.ones_like(f_r, dtype=bool)
    mask220[f_r/MHz < 475] = False
    mask220[f_r/MHz > 2*475] = False
    mask270 = np.ones_like(f_r, dtype=bool)
    mask270[f_r/MHz < 583] = False
    mask270[f_r/MHz > 2*583] = False

    plt.figure()
    plt.plot(f_r/MHz, resobw/kHz, 'k')
    plt.fill_between(f_r/MHz, ymin, ymax, mask150, color='r', alpha=0.2,
    label='150 GHz')
    plt.fill_between(f_r/MHz, ymin, ymax, mask220, color='b', alpha=0.2,
    label='220 GHz')
    plt.fill_between(f_r/MHz, ymin, ymax, mask270, color='g', alpha=0.2,
    label='270 GHz')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Resonator BW [kHz]')
    plt.title("Resonator BW @ %d mK"%(T*1e3))
    plt.legend(loc='lower right', title='Octave of BW')
    plt.grid()
    plt.axis('tight')
    plt.savefig("resobw_vs_resonator_frequency_%dmK.png"%(T*1e3))
    #plt.show()
    #exit()

    ymax = targetspacing.max()/MHz
    ymin = targetspacing.min()/MHz
    plt.figure()
    plt.plot(f_r/MHz, targetspacing/MHz)
    plt.fill_between(f_r/MHz, ymin, ymax, mask150, color='r', alpha=0.2,
    label='150 GHz')
    plt.fill_between(f_r/MHz, ymin, ymax, mask220, color='b', alpha=0.2,
    label='220 GHz')
    plt.fill_between(f_r/MHz, ymin, ymax, mask270, color='g', alpha=0.2,
    label='270 GHz')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Resonator spacing [MHz]')
    plt.title("Target spacing for < 1%% crosstalk @ %d mK"%(T*1e3))
    plt.legend(loc='lower right', title='Octave of BW')
    plt.grid()
    plt.axis('tight')
    plt.savefig("targetspacing_vs_resonator_frequency_%dmK.png"%(T*1e3))
    #plt.show()

    plt.figure()
    plt.plot(f_r/MHz, fr_actual/MHz)
    plt.xlabel('Frequency [MHz]')
    plt.xlabel('Frequency at 380 mK[MHz]')
    plt.grid()
    plt.axis('tight')
    plt.savefig('freqshifted_vs_resonator_frequency_380mK.png')


    plt.figure()
    plt.plot(f_r/MHz, beta)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('beta')
    plt.grid()
    plt.axis('tight')
    plt.savefig('beta_vs_resonator_frequency_%dmK.png'%(T*1e3))

    plt.figure()
    plt.plot(f_r/MHz, geom_delta)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('minimum geometric delta')
    plt.grid()
    plt.axis('tight')
    plt.savefig('min_geom_delta_vs_resonator_frequency_%dmK.png'%(T*1e3))

    plt.figure()
    plt.plot(f_r/MHz, arith_deltaf/kHz)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('minimum arithmetic delta f [kHz]')
    plt.grid()
    plt.axis('tight')
    plt.savefig('min_arith_delta_vs_resonator_frequency_%dmK.png'%(T*1e3))
    #plt.show()
    plt.close('all')




    nu0s = [150, 220, 270]
    f0s = [400., 475., 583.]
    deltas = [0.0055, 0.0004, 0.0003]
    dfs = [2.3, 0.28, 0.24]
    #dfs = [3, 3, 3]
    Ns = [128, 1680, 2400]
    T_ops = [380, 380, 380]

    for i in range(3):
        # 150 GHz case
        N = Ns[i]
        index = np.arange(N)
        f0 = f0s[i]
        delta = deltas[i]
        df = dfs[i]
        print ("delta ", delta)
        print ("df in MHz ", df)
        f_geom = f0*(1+delta)**index*MHz
        f_arith = (f0 + df*index)*MHz

        T_op = T_ops[i]*1e-3
        T = T_op

        f_equallw = np.zeros_like(f_geom)
        f_equallw[0] = f0s[i]*MHz
        for istep in range(1, N):
            x, Q_r = get_xQrMB(T, f_equallw[istep-1])
            delta_fr = f_equallw[istep-1]*(1-x)/(2*Q_r)
            print (istep, delta_fr/MHz)
            f_equallw[istep] = f_equallw[istep-1] + delta_fr/ximax**0.5

        #f_r = f_arith*MHz
        f_r = f_equallw
        omega_r = 2*pi*f_r
        f_g = f_r
        Qc = 10000

        x_geom, Qr_geom = get_xQrMB(T, f_geom)
        x_arith, Qr_arith = get_xQrMB(T, f_arith)
        x_equallw, Qr_equallw = get_xQrMB(T, f_equallw)

        f_geom *= (1-x_geom)
        f_arith *= (1-x_arith)
        f_equallw *= (1-x_equallw)
        dfr_geom = f_geom/(2*Qr_geom)
        dfr_arith = f_arith/(2*Qr_arith)
        dfr_equallw = f_equallw/(2*Qr_equallw)

        crosstalk_geom = 1./(1 + np.diff(f_geom)**2/dfr_geom[:-1]**2)
        crosstalk_arith = 1./(1 + np.diff(f_arith)**2/dfr_arith[:-1]**2)
        crosstalk_equallw = 1./(1 + np.diff(f_equallw)**2/dfr_equallw[:-1]**2)

        Nlw_geom = np.diff(f_geom)/dfr_geom[:-1]
        Nlw_arith = np.diff(f_arith)/dfr_arith[:-1]
        Nlw_equallw = np.diff(f_equallw)/dfr_equallw[:-1]

        plt.figure()
        plt.plot(index, f_geom/MHz, color='k', marker='o', ls='None',
                label='geom')
        plt.plot(index, f_arith/MHz, color='b', marker='s', ls='None',
                label='arith')
        plt.plot(index, f_equallw/MHz, color='r', marker='d', ls='None',
                label='equallw')
        plt.grid()
        plt.xlabel('Index')
        ymin, ymax = plt.ylim()
        if ymax > 1200:
            plt.ylim(bottom=f0-50, top=1200)
        loc='upper left'
        if i==1:
            loc="upper right"
        elif i==2:
            loc="lower right"
        plt.legend(loc=loc)
        plt.ylabel('Frequency [MHz]')
        plt.title('%d GHz frequency schedule'%nu0s[i])
        plt.savefig('%dGHz_freqvsindex.png'%nu0s[i])
        plt.show()

        plt.figure()
        plt.plot(index, Qr_geom, color='k', marker='o', ls='None',
                label='geom')
        plt.plot(index, Qr_arith, color='b', marker='s', ls='None',
                label='arith')
        plt.plot(index, Qr_equallw, color='r', marker='d', ls='None',
                label='equallw')
        plt.grid()
        plt.xlabel('Index')
        loc='lower left'
        if i==2:
            loc="upper right"
        plt.legend(loc=loc)
        plt.ylabel('Qr')
        plt.title('%d GHz frequency schedule'%nu0s[i])
        plt.savefig('%dGHz_Qrvsindex.png'%nu0s[i])
        plt.show()

        plt.figure()
        plt.plot(index[1:], Nlw_geom, color='k', marker='o', ls='None',
                label='geom')
        plt.plot(index[1:], Nlw_arith, color='b', marker='s', ls='None',
                label='arith')
        plt.plot(index[1:], Nlw_equallw, color='r', marker='d', ls='None',
                label='equallw')
        plt.axhline(10, color='k', ls='--')
        plt.grid()
        plt.xlabel('Index')
        loc='upper right'
        if i==1:
            loc="center right"
        elif i==2:
            loc="center right"
        plt.legend(loc=loc)
        plt.ylabel('Number of linewidths')
        plt.title('%d GHz frequency schedule'%nu0s[i])
        plt.savefig('%dGHz_Nlwvsindex.png'%nu0s[i])
        plt.show()

        plt.figure()
        plt.plot(index[1:], crosstalk_geom, color='k', marker='o', ls='None',
                label='geom')
        plt.plot(index[1:], crosstalk_arith, color='b', marker='s', ls='None',
                label='arith')
        plt.plot(index[1:], crosstalk_equallw, color='r', marker='d', ls='None',
                label='equallw')
        plt.axhline(ximax, color='r', ls='--')
        plt.grid()
        plt.xlabel('Index')
        loc='center right'
        if i==1:
            loc="upper left"
        elif i==2:
            loc="upper left"
        plt.legend(loc=loc)
        plt.ylabel('Crosstalk')
        plt.title('%d GHz frequency schedule'%nu0s[i])
        plt.savefig('%dGHz_crosstalkvsindex.png'%nu0s[i])
        plt.show()
        continue

        plt.figure()
        plt.plot(index, f_geom, color='k', marker='o', ls='None', label='geometric spacing')
        plt.plot(index, f_arith,color='b', marker='s', ls ='None', label='arithmetic spacing')
        plt.grid()
        plt.xlabel('Index')
        plt.ylabel('Frequency [MHz]')
        plt.title('%d GHz frequency schedule'%nu0s[i])
        plt.legend(loc='best')
        plt.show()




elif analyzeTempandNumber:

    ximax = 0.01
    nu0s = [150, 220, 270]
    f0s = [400., 475., 583.]
    deltas = [0.0055, 0.0004, 0.0003]
    dfs = [2.3, 0.28, 0.24]
    #dfs = [3, 3, 3]
    Ns = [128, 1680, 2400]

    for iband in range(3):
        BW = f0s[iband]*MHz
        temps = np.r_[250:450:1000j]*1e-3
        bounds = []
        N = np.arange(100, 2500,10)
        for nreso in N:
            f0 = BW
            df = BW/nreso
            f = f0 + df*np.arange(nreso)
            X, Y = np.meshgrid(temps, f)

            x, dQi = get_xQMB(X, Y)
            Qi = 1./dQi
            Qc = 2.0*Qi
            Qr = Qi*Qc/(Qi+Qc)
            Nlw = 2*Qr*df/f[:, np.newaxis]
            print (nreso, np.max(Nlw),np.min(Nlw))

            valid = Nlw >= 1./np.sqrt(ximax)
            row = valid[-1, :]
            bounds.append(row)

            #plt.figure()
            #plt.pcolormesh(X/1e-3, Y/MHz, valid, cmap='gray')
            #plt.xlabel('Temperature [mK]')
            #plt.ylabel('Frequency [MHz]')
            #plt.show()


        print ("\n")
        bounds = np.asarray(bounds)
        print (bounds.shape)
        X, Y = np.meshgrid(temps, N)
        plt.figure()
        plt.pcolormesh(X/1e-3, Y, bounds, cmap='gray')
        plt.grid()
        plt.axhline(Ns[iband] , color='k', ls='-.', lw=2, label="%d"%Ns[iband])
        plt.axhline(Ns[iband]//2, color='b', ls='-', lw=2,
                label="%d"%(Ns[iband]//2))
        plt.axhline(Ns[iband]//3, color='r', ls='-.', lw=2,
                label="%d"%(Ns[iband]//3))
        plt.axhline(Ns[iband]//4, color='g', ls='--', lw=2,
                label="%d"%(Ns[iband]//4))
        plt.axhline(Ns[iband]//5, color='m', ls=':', lw=2,
                label="%d"%(Ns[iband]//5))
        plt.axvline(380, color='r', ls='--', lw=2)
        plt.legend(loc='upper right', title="Nreso")
        plt.xlabel('Temperature [mK]')
        plt.ylabel('Number of resonators')
        plt.title('%d GHz frequency schedule'%nu0s[iband])
        plt.savefig('%dGHz_numberofresospacked_vs_temperature.png'%nu0s[iband])
        plt.show()
        #exit()





