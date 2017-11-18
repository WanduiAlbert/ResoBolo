import numpy as np
from astropy.constants import c, eps0, mu0
import astropy.units as u
import pandas as pd

# Calculations of the coupling capacitance for the resobolos

# Define some useful dimensions
W_r = 1001.5 * u.um # Width of the resonator capacitor
W_c = 76.5 * u.um # width of the coupling capacitor

w = 2 * u.um # width of the finger
g = 2 * u.um # gap between fingers

Cs = 1.1297e-4 * u.pF # Surface capacitance in pF/sq. From sonnet simulations

def N_sq(L, W):
	return (L*W/w/(w + g)).to(1)

def C(L, W):
	return Cs * N_sq(L, W)

# Determined from calculations + simulations of the resobolo at T=0.38K
Qi = 11492.87
C1 = 0.19 * u.pF
Inductance = 9.99 * u.nH # Both surface and kinetic
Z0 = 50 * u.Ohm

# resonant frequency
def f_r(L, W):
	return (1/((Inductance * C(L, W))**0.5 * 2 * np.pi)).to('MHz')

# compute the optimal coupling capacitance as a function of the dimensions
# of the resonator capacitor
def C_c(L, W):
	Ci = C(L, W)
	Cc = (4*Inductance*Ci**3/(Z0**2 * Qi**2))**0.25
	return (2 * Cc).to('pF') # I want C1 = C2 = 2 * Cc

def coupling_l(L, W):
	C1 = C_c(L, W) # The value of the coupling capacitance either to ground or to feedline

	return (C1/Cs).to(1) *(w * (w + g)/W_c).to('um')


if __name__ == "__main__":
	L = np.array([1346, 1462, 1598, 1750, 1930]) * u.um # Lengths of the 5 resobolos

	C_r = C(L, W_r)
	f_rs = f_r(L, W_r)

	C1s = C_c(L, W_r)
	L_c = coupling_l(L, W_r)

	# Let's organize this data into a table
	summary_data = np.array([L, C_r, f_rs, C1s, L_c]).T

	df = pd.DataFrame(summary_data, columns=['Length of C_r [um]', 'C_r [pF]', 'f_r [MHz]', 'C_1 [pF]', 'Length of C_1 [um]'], index=[1,2,3,4,5])

	df.to_csv('reso_cap_params.csv')