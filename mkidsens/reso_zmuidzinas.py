
import numpy as np
import pint

u = pint.UnitRegistry()
k = u.boltzmann_constant
hbar = u.hbar

# Aluminum
N0 = 1.72e10 / ( u.electron_volt * u.micrometer**3 )
Tc = 1.4 * u.K
nstar = 100 / u.micrometer**3
Ft = 3.5

Delta = 1.76 * k * Tc

# Design
Pg = 1.0 * u.picowatt
Po = 8.0 * u.picowatt
V = 1000 * u.micrometer**3
Tamp = 3.0 * u.K

Hnu = 1.0

dfg = 0.0
dfr = 1.0
chic = 4 * Qr**2 / (Qc * Qi)
chig = 1.0 / (1.0 + (dfg/dfr)**2)
chiread = 0.5 * chic * chig
Pread = chiread * Pg

Pleg = Pread + Po

gammaleg = 3.0

Gammareadgen = etaread * Pread / Delta

T0 = 0.27
Tb = 0.4
#Pleg = Kleg * (T**gammaleg - T0**gammaleg)
Kleg = Pleg / (T**gammaleg - T0**gammaleg)
#gammaG = (Tb/Pleg) * dPleg / dTb
gammaG = (Tb/Pleg) * Kleg * (gammaleg-1)*T**gammaleg
Gb = gammaG * Pleg / Tb

chiph = 1.0
Sph = 2 * chiph

nth = 2 * N0 * np.sqrt(2*pi*k*Tb*Delta) * np.exp(-Delta/(k*Tb))
R = 1.0 / (nstar * taumax)
nqp = np.sqrt((nth + nstar)**2 + 2 * Gammareadgen / (V*R)) - nstar

tauqp = taumax / (1 + nqp / nqpstar)
tauth = taumax / (1 + nth / nqpstar)

kappa = 0.5 * (1 + Ft*Tc/Tb)

NEPamp2 = Hnu * (Pb * tauth * 4 / (kappa * tauqp * chic * chiqp * beta))**2 * k*Tamp / Pg
NEPshot2 = 2*Hnu * (Pb*tauth / (kappa * Nqp)) * Sshot
