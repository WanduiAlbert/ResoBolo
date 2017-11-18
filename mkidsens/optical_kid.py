#!/usr/bin/env python

'''
Generation-recombination noise prediction for a lumped element kinetic inductance detector with no optical load
Source: Zmuidzinas - Superconducting Microresonators for Bolometer Readout
'''

from math import pi
import pint
import numpy as np
from scipy import special
import matplotlib.pyplot as pl

u = pint.UnitRegistry()

# physical constants
k = u.boltzmann_constant
h = u.planck_constant

# material properties
Tc = 1.4 * u.kelvin
N0 = 1.72e10 / (u.electron_volt * u.micrometer**3)
taumax = 0.0005 * u.second
nstar = 100 / (u.micrometer**3)
delta0 = 1.76 * k * Tc
rho = 3.0 * u.microohm * u.cm
t = 0.04*u.micrometer

'''
Tc = 1.4 * u.kelvin
N0 = 8.9e9 / (u.electron_volt * u.micrometer**3)
taumax = 0.0001 * u.second
nstar = 100 / (u.micrometer**3)
delta0 = 1.76 * k * Tc
rho = 100.0 * u.microohm * u.cm
t = 0.25 * u.micrometer
'''

Rs = rho / t
Lk = (h * Rs / (2 * pi*pi * delta0)).to('picohenry')
print "Lk: ",Lk,"/ square"

#T = np.linspace(0.1,0.7,100) * u.kelvin	
Tbath = 0.27*u.kelvin

# optical design
nu = 270 * u.gigahertz
bw = 0.27
etadewar = 0.4
Tdewar = 31 * u.kelvin
#Tdewar = 0 * u.kelvin
etahorn = 0.75

# resonator design
fro = 232.0 * u.megahertz	# readout frequency
#Lk = 2.896 * u.nanohenry
Lk = 3.0 * u.nanohenry
Lg = 4.0 * u.nanohenry
#Lg = (1.69 + 1.4) * u.nanohenry
L = Lg + Lk # total inductance
Tamp = 3.0 * u.kelvin			# amplifier noise temperature
Qiloss = 3e5				# Intrinsic Q of unloaded resonator
Cc = 0.4 * u.picofarad		# Coupling capacitor capacitance
Z0 = 50.0 * u.ohms			# Readout microstrip impedance
V = 240 * u.micrometer**3
#V = 3.5 * (5/6.5) *2* 800 * .5 * u.micrometer**3
print "V: ",V
Pg = 1.0 * u.picowatt			# Readout tone power
filmgamma = 1.0		# thin film
chic = 1.0
niter = 5
Lk = L - Lg							# kinetic inductance
alphak = Lk / L						# kinetic inductance fraction
print "alphak: ",alphak
C = (1.0 / ((2*pi*fro)**2 * L)).to('picofarad')		# resonator capacitance
Qc = (C / (pi * fro * Cc**2 * Z0)).to('dimensionless')	# coupling Q
nqpth = (2 * N0 * np.sqrt(2. * pi * k * Tbath * delta0) * np.exp(-delta0 / (k * Tbath))).to(1./u.micrometer**3)
etao = 0.7				# Efficiency of optical photons to break pairs
chiqp = 1.0
etaa = 0.0				# Efficiency of readout power to break pairs

print "IDC capacitor: ",C

def update_resonator():
	global Po,Tdewar,etadewar,nu,bw
	global gammao,etao,etahorn,delta0
	global gammaa,chia,Pg,Pa,chic
	global nqp,tauqp,taue,xMB
	global s21_gr_psd_re, s21_gr_psd_im, s21_amp
	global Qi,Qiqp,Qr
	Po = k * Tdewar * etadewar * nu * bw
	gammao = Po * etao * etahorn / delta0

	for i in range(niter):
		chia = 0.5*chic
		Pa = chia * Pg
		gammaa = (etaa * Pa / delta0).to('Hz')
		gamma = gammao + gammaa

		nqp = np.sqrt((nqpth + nstar)**2 + 2 * gamma * nstar * taumax / V) - nstar

		hf2kt = (h * fro / (2 * k * Tbath)).to('dimensionless')
		S1 = (2./pi) * np.sqrt(2 * delta0 / (pi * k * Tbath)) * np.sinh(hf2kt) * special.k0(hf2kt)

		Qiqp = (2 * N0 * delta0 / (alphak * S1 * nqp * filmgamma)).to('')
		Qi = 1. / ((1./Qiqp) + (1./Qiloss))
		Qr = 1. / ((1./Qi) + (1./Qc))
		chic = 4*Qr*Qr/(Qi*Qc)
	
	'''
	print Pa
	print gammaa
	print nqpth
	print nqp
	print Qc
	print Qi
	'''

	S2 = 1.0 +  np.sqrt(2 * delta0 / (pi * k * Tbath)) * np.exp(-hf2kt) * special.i0(hf2kt)
	beta = S2 / S1

	xMB = -alphak * S2 * nqp / (4 * N0 * delta0)

	tauqp = taumax / (1. + nqp/nstar)
	nqpo = (np.sqrt(gammao * nstar * taumax / V)).to('micrometer^-3')
	nqpa = (np.sqrt(gammaa * nstar * taumax / V)).to('micrometer^-3')
	taue = Qr / (pi * fro)

	tauth = taumax / (1. + nqpth/nstar)
	gammath = 0.5 * nqpth * V * (1.0/taumax + 1.0/tauth)
	gammarec = 0.5 * nqp * V * (1.0/taumax + 1.0/tauqp)

	chiqp = Qi/Qiqp
	shot_psd = 2*nqp*V*(1/taumax + 1/tauqp)
	Nqp_psd = tauqp**2 * shot_psd
	s21_gr_psd_re = (((chic*chiqp/4)*(1/(nqp*V)))**2 * Nqp_psd).to('1/Hz')
	s21_gr_psd_im = s21_gr_psd_re*beta**2

	s21_amp = (k * Tamp / (2 * Pg)).to('1/Hz')

	nep_ogen2 = (4 * gammao * delta0**2 / etao**2).to('W^2/Hz')
	nep_elgen2 = (4 * gammaa * delta0**2 / etao**2).to('W^2/Hz')
	nep_thgen2 = (4 * gammath * delta0**2 / etao**2).to('W^2/Hz')
	nep_rec2 = (4 * gammarec * delta0**2 / etao**2 ).to('W^2/Hz')
	nep_gr2 = nep_rec2 + nep_thgen2 + nep_elgen2

	nep_elgen = np.sqrt(nep_elgen2).to('aW/Hz**0.5')
	nep_thgen = np.sqrt(nep_thgen2).to('aW/Hz**0.5')
	nep_rec = np.sqrt(nep_rec2).to('aW/Hz**0.5')
	nep_gr = np.sqrt(nep_gr2).to('aW/Hz**0.5')

def db(x):
	return 10*np.log10(x*u.Hz)/u.Hz

Tdewar0 = Tdewar
Tdewar1 = Tdewar0 + 0.1 * u.kelvin
Tdewar = Tdewar0
update_resonator()
nqp0 = nqp
xMB0 = xMB
Tdewar = Tdewar1
update_resonator()
nqp1 = nqp
xMB1 = xMB

dxMBdT = (xMB1 - xMB0) / (Tdewar1 - Tdewar0)
dxMBdT = dxMBdT.to('1/kelvin')
print "dxMBdT:",dxMBdT

'''
Tdewars = np.linspace(0,40,100)*u.K
nqps = []
for x in Tdewars:
	Tdewar = x
	update_resonator()
	nqps.append(nqp.m)
nqps = np.array(nqps)

pl.plot(Tdewars.m,nqps)
pl.show()
exit()
'''

print "Qc: ",Qc
print "Qi: ",Qi
print "Qr: ",Qr
print "df: ",6.0*fro/Qr
print "chic: ",chic

dnqpdT = (nqp1 - nqp0) / (Tdewar1 - Tdewar0)
print dnqpdT

s21_gr_psd_re = db(s21_gr_psd_re).m
s21_gr_psd_im = db(s21_gr_psd_im).m
s21_amp = db(s21_amp).m

f3dbqp = (1.0 / (2 * pi * tauqp)).to('hertz')
f3dbe = (1.0 / (2*pi*taue)).to('hertz')

print 'Tbath: ',Tbath
print 'f3dbqp: ',f3dbqp
print 'f3dbe: ',f3dbe
print 's21 gr psd re: ',s21_gr_psd_re
print 's21 gr psd im: ',s21_gr_psd_im

'''
pl.figure()
pl.plot(Tbath,f3dbqp,label='quasiparticle')
pl.plot(Tbath,f3dbe,label='electrical')
pl.xlabel('Temperature (K)')
pl.ylabel('f3db (Hz)')
pl.grid()
pl.gca().set_yscale('log')
pl.title('Quasiparticle time constant')
pl.legend()

pl.figure()
pl.plot(Tbath,chic)
pl.xlabel('Temperature (K)')
pl.ylabel('chi_c')
pl.grid()
pl.title('Readout coupling efficiency')

pl.figure()
pl.plot(Tbath,s21_gr_psd_re,label='Amplitude readout')
pl.plot(Tbath,s21_gr_psd_im,label='Phase readout')
pl.axhline(s21_amp,color='black',label='Cold amp noise (%.1fK)'%Tamp.m)
pl.xlabel('Temperature (K)')
pl.ylabel('PSD (dBc/Hz)')
pl.title('KID noise vs temperature at 1pW bias')
pl.grid()
pl.legend(loc='upper right')
pl.ylim(ymin=-150,ymax=-60)
pl.show()

'''

