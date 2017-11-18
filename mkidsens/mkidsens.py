
from math import pi
import pint
import numpy as np
from scipy import special

u = pint.UnitRegistry()

''' Sensitivity calculation from McKenney et al '''

# physical constants
k = 1.38e-23 * u.joule / u.kelvin
h = 6.63e-34 * u.joule * u.second

# material properties
Tc = 2.5 * u.kelvin
N0 = 8.9e9 / (u.electron_volt * u.micrometer**3)
taumax = 0.00004 * u.second
nstar = 100 / (u.micrometer**3)
delta0 = 1.76 * k * Tc		# gap energy
etapb = 0.57				# pair breaking efficiency

# optical design
nu = 270 * u.gigahertz
bw = 0.27
etadewar = 0.4
Tsky = 3. * u.kelvin
Tdewar = 5. * u.kelvin
Po = k * (Tsky + Tdewar) * etadewar * nu * bw
etahorn = 0.75
print "eff: ",etadewar * etahorn

# resonator design
T = 0.3 * u.kelvin		# operating temperature
fro = 1.0 * u.gigahertz		# readout frequency
d = 400.0 * u.micrometer	# waveguide radius
w = 2.0 * u.micrometer		# trace width
t = 0.02 * u.micrometer		# trace thickness
L = 4.1 * u.nanohenry		# total inductance
Lg = 0.8 * u.nanohenry		# geometric inductance
Tn = 4.0 * u.kelvin			# amplifier noise temperature
Qi = 1e6
Cc = 0.12 * u.picofarad
Z0 = 50.0 * u.ohms
Pa = 3 * u.picowatt
chic = 0.5

# resonator derived properties
l = d * 2 * 2				# trace length
Lk = L - Lg					# kinetic inductance
alpha = Lk / L				# kinetic inductance fraction
C = 1.0 / ((2*pi*fro)**2 * L)
Qc = (C / (pi * fro * Cc**2 * Z0)).to('dimensionless')	# coupling Q
V = w * t * l
gammao = (Po * etapb * etahorn / delta0).to('gigahertz')
gammaa = (Pa * chic / delta0).to('gigahertz')
gamma = gammao + gammaa
print "gammaa: ",gammaa
print "gammao: ",gammao
nqpth = (2 * N0 * np.sqrt(2. * pi * k * T * delta0) * np.exp(-delta0 / (k * T))).to(1./u.micrometer**3)
nqp = np.sqrt((nqpth + nstar)**2 + 2 * gamma * nstar * taumax / V) - nstar
print "nqp: ",nqp

hf2kt = (h * fro / (2 * k * T)).to('dimensionless')
y = np.sqrt(2*delta0 / (pi * k * T))
S1 = (2./pi) * np.sqrt(2 * delta0 / (pi * k * T)) * np.sinh(hf2kt) * special.k0(hf2kt)
S2 = 1.0 +  np.sqrt(2 * delta0 / (pi * k * T)) * np.exp(-hf2kt) * special.i0(hf2kt)

beta = S2 / S1

tauqp = taumax / (1. + nqp / nstar)

nepamp2 = (nqp * V * delta0 / (beta * etahorn * etapb * tauqp))**2 * (8 * k * Tn / (chic * Pa))
nepamp2 = nepamp2.to(u.attowatt**2 / u.hertz)
nepamp = np.sqrt(nepamp2)

nepshot2 = (2 * h * nu * Po / (etapb * etahorn)).to(u.attowatt**2 / u.hertz)
nepshot = np.sqrt(nepshot2)
nepwave2 = (2*Po*Po/(nu * bw)).to(u.attowatt**2 / u.hertz)
nepwave = np.sqrt(nepwave2)
nepphot = np.sqrt(nepshot2 + nepwave2)

nepgr2 = (2 * Po * delta0 / (etapb * etahorn)).to(u.attowatt**2 / u.hertz)
nepgr =  np.sqrt(nepgr2)

Qqp = (2 * N0 * delta0 / (alpha * nqp * S1)).to('dimensionless')
Qr = 1.0 / (1./Qqp + 1./Qc)
print "Qc: ",Qc
print "Qqp: ",Qqp
print "Qr: ",Qr

Ec = 2 * N0 * delta0**2 * V
pbif = (0.8 * Qc * (2*pi*fro) * Ec / (2 * Qr**3)).to('picowatt')
print "Pbifurcation: ",pbif

print "amplifier noise: ",nepamp
print "g-r noise: ",nepgr
print "photon noise: ",nepphot

neptotal = np.sqrt(nepamp**2 + nepphot**2 + nepgr**2)
print "total nep: ",neptotal


def dpdt_fun(fs,bw,t):
	x =  h * fs / k
	eu = np.exp(x/t)
	dpdt = bw * h * fs * eu * x / (t * (eu-1.))**2
	dpdt = dpdt.to(u.picowatt / u.kelvin)
	return dpdt

dpdtcmb = etadewar * dpdt_fun(nu,nu*bw,2.725*u.kelvin)
dpdtrj = etadewar * dpdt_fun(nu,nu*bw,300*u.kelvin)

print "dpdtrj / dpdtcmb",dpdtrj/dpdtcmb

ukrts = u.microkelvin * u.second**0.5
netcmb = (neptotal / (2**0.5*dpdtcmb)).to(ukrts)
netrj = (neptotal / (2**0.5*dpdtrj)).to(ukrts)

print "netcmb: ",netcmb
print "netrj: ",netrj



