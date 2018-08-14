#! /usr/bin/env python3

import numpy as np
from math import pi
from scipy.constants import epsilon_0, mu_0,c
import matplotlib.pyplot as plt

# Unit conversions
Z0 = 50 #Ohms
mm = 1e-3
um = 1e-6
pF = 1e-12
nH = 1e-9
pH = 1e-12
MHz = 1e6



def S21(nu, Zc, v_ph, l):
    """
    S21 (nu, Zc, v_ph, l)
    returns the complex forward transmission coefficient given:
        nu - frequency of interest
        Zc - impedance of transmission line in Ohms. Zc = sqrt(L/C)
        v_ph - phase velocity of waves on line v_ph = 1/sqrt(L*C) in m/s
        l - electrical length of the transmission line in m
    """
    beta = 2*pi*nu/v_ph
    return (2*Z0*Zc)/((Zc + Z0)**2 * np.exp(1j*beta*l) + (Zc- Z0)**2*np.exp(-1j*beta*l))

def S11(nu, Zc, v_ph, l):
    """
    S11 (nu, Zc, v_ph, l)
    returns the complex reflection coefficient given:
        nu - frequency of interest
        Zc - impedance of transmission line in Ohms. Zc = sqrt(L/C)
        v_ph - phase velocity of waves on line v_ph = 1/sqrt(L*C) in m/s
        l - electrical length of the transmission line in m
    """
    beta = 2*pi*nu/v_ph
    return (1j*(Zc**2- Z0**2)*np.tan(beta*l))/(2*Zc*Z0 + 1j*(Zc**2- Z0**2)*np.tan(beta*l))


if __name__=="__main__":
    # Do a calculation for the case of a simple microstrip line with the
    # following parameters:
    # w = 2 um
    # d = 300nm
    # dielectric: SiO2 with er=3.9
    # surface kinetic inductance Lk = 0.158pH/sq
    er = 3.9
    d = 0.3*um
    w = 2*um
    Lk = 0.158*pH
    Cs = epsilon_0*(1 + er)/2*w/d
    Ls = mu_0 * d/w + Lk*w

    Zc = (Ls/Cs)**0.5
    v_ph = 1/(Ls*Cs)**0.5
    #er_eff = 3.35
    #v_ph = c/er_eff**0.5
    #Zc = 29
    print ("Characteristic impedance: {0:3.1f} Ohms".format(Zc))
    print ("Wave velocity: {0:3.1f} c".format(v_ph/c))

    nu = np.r_[50:1000:1000j]*MHz
    l = 22*mm

    S = S21(nu, Zc, v_ph, l)

    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(nu/MHz, 20*np.log10(np.abs(S)))
    ax.axis('tight')
    ax.grid(which='both')
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|S21|')
    plt.show()


