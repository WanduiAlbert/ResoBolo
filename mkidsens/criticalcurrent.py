#! /usr/bin/env python3

import numpy as np
import astropy.units as u

from astropy.constants import c, eps0, mu0, h, k_B

from math import pi

hbar = h/2/pi

# These are the calculations of the critical current for the new TKID transmission lines. since the widths are vastly smaller
# we expect that there is only so much current we can handle before non-linearities become important. My goal here is to calculate
# the critical current and readout power.

t = 400 * u.nm
w = 2 * u.um

l = 20*u.mm

Tc = 9.3 * u.K # Nb critical temperature

Delta = 3.5/2 * k_B * Tc

Ls = 0.158 * u.pH # surface kinetic inductance of Nb. I'll use this to back out the resistivity of Nb.

rho_n = (Ls * (0.15 * u.um) * pi* Delta/hbar).to('uOhm cm')

print ("The resistivity of Niobium is {0:1.3f}".format(rho_n))

gamma = 7.53 * u.mJ /u.mol/u.K**2 # From the paper: https://journals.aps.org/pr/abstract/10.1103/PhysRev.127.1501

A_r = 92.906 * u.g/u.mol
density = 8.57 * u.g/u.cm**3

N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k_B**2)).to('1/(J um^3)')

J_star = np.sqrt(pi*N_0*Delta**3/(hbar*rho_n)).to('A/um^2')

I_c = (0.42 * J_star * w * t).to('mA')

print ("The critical current of the line is {0:1.4f}".format(I_c))

# For a 50 Ohm source

R = 50 * u.Ohm

P = (I_c**2 * R).to('mW')

print ("The critical readout power is {0:1.4f}".format(P))

dBm  = u.dB(u.mW)

P_dB = 10 * np.log10((P/(1*u.mW)).to(1).value) * dBm

print ("In dBm this critical power is {0:1.4f}".format(P_dB))