
import math
from math import pi
import pint
import numpy as np
from scipy import special

u = pint.UnitRegistry()

''' Sensitivity calculation from McCarrick et al '''

# physical constants
k = 1.38e-23 * u.joule / u.kelvin
h = 6.63e-34 * u.joule * u.second

# material properties
Tc = 1.4 * u.kelvin
N0 = 1.72e10 / (u.electron_volt * u.micrometer**3)
taumax = 0.0005 * u.second
nstar = 100 / (u.micrometer**3)
delta0 = 1.76 * k * Tc		# gap energy
etapb = 0.57				# pair breaking efficiency

# optical design
nu = 270 * u.gigahertz
bw = 0.27
etadewar = 0.4
Tsky = 18. * u.kelvin
Tdewar = 10. * u.kelvin
Po = k * (Tsky + Tdewar) * etadewar * nu * bw
etahorn = 0.7

# resonator design
T = 0.3 * u.kelvin		# operating temperature
fro = 0.165 * u.gigahertz		# readout frequency
d = 400.0 * u.micrometer	# waveguide radius
w = 2.0 * u.micrometer		# trace width
t = 0.02 * u.micrometer		# trace thickness
L = 4.1 * u.nanohenry		# total inductance
Lg = 0.8 * u.nanohenry		# geometric inductance
Tn = 4.0 * u.kelvin			# amplifier noise temperature
Qi = 1e6					# Intrinsic Q of unloaded resonator
Cc = 0.3 * u.picofarad		# Coupling capacitor capacitance
Z0 = 50.0 * u.ohms			# Readout microstrip impedance
Pa = 10 * u.picowatt		# Readout tone power
chic = 0.5					# Fraction of readout power that goes into MKID

# resonator derived properties
l = d * 2 * 2				# trace length
Lk = L - Lg					# kinetic inductance
alpha = Lk / L				# kinetic inductance fraction
C = 1.0 / ((2*pi*fro)**2 * L)
Qc = (C / (pi * fro * Cc**2 * Z0)).to('dimensionless')	# coupling Q
V = w * t * l
chia = 0.0		# Efficiency of readout power at breaking cooper pairs
gammao = Po * etapb * etahorn / delta0
gammaa = chia * Pa * chic / delta0
gamma = gammao + gammaa
nqpth = (2 * N0 * np.sqrt(2. * pi * k * T * delta0) * np.exp(-delta0 / (k * T))).to(1./u.micrometer**3)
nqp = np.sqrt((nqpth + nstar)**2 + 2 * gamma * nstar * taumax / V) - nstar

nqpo = (np.sqrt(gammao * nstar * taumax / V)).to('micrometer^-3')
nqpa = (np.sqrt(gammaa * nstar * taumax / V)).to('micrometer^-3')

print "nqpth: ",nqpth
print "nqpo: ",nqpo
print "nqpa: ",nqpa
print "nqp: ",nqp

hf2kt = (h * fro / (2 * k * T)).to('dimensionless')
y = np.sqrt(2*delta0 / (pi * k * T))
S1 = (2./pi) * np.sqrt(2 * delta0 / (pi * k * T)) * np.sinh(hf2kt) * special.k0(hf2kt)
S2 = 1.0 +  np.sqrt(2 * delta0 / (pi * k * T)) * np.exp(-hf2kt) * special.i0(hf2kt)

beta = S2 / S1

tauqp = taumax / np.sqrt(1. + 2 * gamma * taumax / (nstar * V))

Qqp = 2 * N0 * delta0 / (alpha * S1 * nqp)
Qr = 1. / ((1. / Qqp) + (1. / Qi) + (1./Qc))

print "Qqp: ",Qqp.to('dimensionless')
print "Qc: ",Qc
print "Qr: ",Qr

b32 = (etapb*etahorn*taumax / (delta0 * V * nstar)).to('seconds / electron_volt')
dxdnqp = (alpha * S2 / ( 4 * N0 * delta0)).to(u.micrometer**3)
dnqpdPo = (nstar * b32 / (2 * math.sqrt(1 + Po * b32))).to('1 / (picowatt * um^3)')
dxdPo = (dxdnqp * dnqpdPo).to('1/ picowatt')
dfdPo = (fro * dxdPo).to('gigahertz / picowatt')

efamp = (np.sqrt(4 * k * Tn / Pa) * (Qc/Qr**2) * fro).to(u.hertz / u.hertz**0.5)
print "efamp: ",efamp

nepamp = (efamp / dfdPo).to('attowatt / hertz ** 0.5')

nepamp2 = (nqp * V * delta0 / (beta * etahorn * etapb * tauqp))**2 * (8 * k * Tn / (chic * Pa))
nepamp2 = nepamp2.to(u.attowatt**2 / u.hertz)
nepamp_alt = np.sqrt(nepamp2)

print "nepamp McCarrick: ",nepamp
print "nepamp McKenney: ",nepamp_alt

nepshot2 = (2 * h * nu * Po / (etapb * etahorn)).to(u.attowatt**2 / u.hertz)
nepshot = np.sqrt(nepshot2)
nepwave2 = (2*Po*Po/(nu * bw)).to(u.attowatt**2 / u.hertz)
nepwave = np.sqrt(nepwave2)
nepphot = np.sqrt(nepshot2 + nepwave2)

nepgr2 = (2 * Po * delta0 / (etapb * etahorn)).to(u.attowatt**2 / u.hertz)
nepgr =  np.sqrt(nepgr2)

print "amplifier noise: ",nepamp
print "g-r noise: ",nepgr
print "photon noise: ",nepphot

neptotal = np.sqrt(nepamp**2 + nepphot**2 + nepgr**2)
print "total nep: ",neptotal


def dpdt_fun(fs,t):
	u =  h * fs / k
	eu = np.exp(u/t)
	dpdt = h * fs * eu * u / (t * (eu-1.))**2
	return dpdt
dpdtcmb = dpdt_fun(nu,2.725*u.kelvin)
dpdtrj = dpdt_fun(nu,300*u.kelvin)

dtrjdtcmb = dpdtrj / dpdtcmb

