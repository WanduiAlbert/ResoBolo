#!/usr/bin/env python

'''
Generation-recombination noise prediction for a lumped element kinetic inductance detector with no optical load
Source: Zmuidzinas - Superconducting Microresonators for Bolometer Readout
'''

import sys
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

# Compute KID properties for a range of bath temperatures
T = np.linspace(0.1,0.7,100) * u.kelvin	
#T = 0.5*u.kelvin

# resonator design
fro = 232.0 * u.megahertz	# readout frequency
L = 9.0 * u.nanohenry		# total inductance
Lg = (2.35+2.75) * u.nanohenry		# geometric inductance
Tamp = 4.0 * u.kelvin			# amplifier noise temperature
Qiloss = 3e5				# Intrinsic Q of unloaded resonator
Cc = 0.51 * u.picofarad		# Coupling capacitor capacitance
Z0 = 50.0 * u.ohms			# Readout microstrip impedance
V = 240 * u.micrometer**3
Pg = 1 * u.picowatt			# Readout tone power


# resonator derived properties
Lk = L - Lg							# kinetic inductance
alphak = Lk / L						# kinetic inductance fraction
C = (1.0 / ((2*pi*fro)**2 * L)).to('pF')		# resonator capacitance
Qc = (C / (pi * fro * Cc**2 * Z0)).to('dimensionless')	# coupling Q
nqp = (2 * N0 * np.sqrt(2. * pi * k * T * delta0) * np.exp(-delta0 / (k * T))).to(1./u.micrometer**3)

print "fro: ",fro
print "C: ",C
print "L: ",L
print "coupling Q: ",Qc

q = (h * fro / (2 * k * T)).to('dimensionless')
v = np.sqrt(2*delta0 / (pi * k * T))
S1 = (2./pi) * v * np.sinh(q) * special.k0(q)
S2 = 1.0 +  v * np.exp(-q) * special.i0(q)
beta = S2 / S1

tauqp = taumax / (1. + nqp/nstar)

filmgamma = 1.0		# thin film
Qiqp = (2 * N0 * delta0 / (alphak * S1 * nqp * filmgamma)).to('')
Qi = 1. / ((1./Qiqp) + (1./Qiloss))
Qr = 1. / ((1./Qi) + (1./Qc))

taue = Qr / (pi * fro)

# Note typo in Zmuidzinas Eq. 9 missing square
chic = 4 * Qr**2 / (Qc * Qi)

Nqp = nqp * V
chiqp = Qi/Qiqp
shot_psd = 2*Nqp*(1/taumax + 1/tauqp)
Nqp_psd = tauqp**2 * shot_psd
s21_gr_psd_re = (((chic*chiqp/4)*(1/Nqp))**2 * Nqp_psd).to('1/Hz')
s21_gr_psd_im = s21_gr_psd_re*beta**2
s21 = 1.0 - Qr / Qc

s21_amp = (k * Tamp / (2 * Pg)).to('1/Hz')

def db(x):
	return 10*np.log10(x*u.Hz)/u.Hz

s21_gr_psd_re = db(s21_gr_psd_re).m
s21_gr_psd_im = db(s21_gr_psd_im).m
s21_amp = db(s21_amp).m
s21 = db(s21/u.Hz).m

f3dbqp = (1.0 / (2 * pi * tauqp)).to('hertz')
f3dbe = (1.0 / (2*pi*taue)).to('hertz')

#print T
#print 'f3dbqp: ',f3dbqp
#print 'f3dbe: ',f3dbe
#print 's21 gr psd re: ',s21_gr_psd_re
#print 's21 gr psd im: ',s21_gr_psd_im

'''
pl.figure()
pl.plot(T,f3dbqp,label='quasiparticle')
pl.plot(T,f3dbe,label='electrical')
pl.xlabel('Temperature (K)')
pl.ylabel('f3db (Hz)')
pl.grid()
pl.gca().set_yscale('log')
pl.title('Quasiparticle time constant')
pl.legend()

pl.figure()
pl.plot(T,chic)
pl.xlabel('Temperature (K)')
pl.ylabel('chi_c')
pl.grid()
pl.title('Readout coupling efficiency')
'''


# Comparing to real data

fn = sys.argv[1]
#Tm,_,sig_on,noisere_on,noiseim_on,sig_off,noisere_off,noiseim_off = np.loadtxt(fn).T
Tm,carrpower,noiseim_on = np.loadtxt(fn).T

refpower = carrpower[0]
Tm = Tm[1:]
carrpower = carrpower[1:]
noiseim_on = noiseim_on[1:]


def todb(x):
	return 10*np.log10(x)
def fromdb(x):
	return 10**(0.1*x)

pred_noisefloor = todb(fromdb(s21_gr_psd_im) + fromdb(s21_amp))

s21_meas = 0.5*(carrpower - refpower)

# s21_meas is the carrier amplitude
# to get the carrier power to noise ratio, multiply carrier amplitude by two
snr_on_meas = -carrpower + noiseim_on
snr_on_pred = -2*s21 + pred_noisefloor

#s21_meas = sig_on - sig_off + 3		# The on signal was biased at 1/2 the power of the off signal
#s21_meas = np.interp(T,Tm,s21_meas)

pl.figure()
pl.scatter(Tm,s21_meas,label='Measured')
pl.plot(T,s21,label='Predicted')
pl.xlabel('Temperature (K)')
pl.ylabel('|S21| (dB)')
pl.grid()
pl.title('measured S21 vs predicted S21')
pl.legend(loc='lower right')
pl.savefig('measvspreds21.png')

pl.figure()
pl.scatter(Tm,snr_on_meas,color='b',label='meas noise')
#pl.scatter(Tm,sig_off-noisere_off,color='g',label='meas off res SNR')
pl.plot(T,snr_on_pred,'b-',label='pred noise')
#pl.plot(T,-s21_amp+np.zeros_like(T),'g-',label='pred off res SNR')
pl.grid()
pl.xlabel('Temperature (K)')
pl.ylabel('Noise to carrier ratio (dBc/Hz)')
pl.title('Measured noise of 232MHz dark resonator')
pl.legend(loc='upper right')
pl.savefig('measvsprednoise.png')
pl.show()
exit()


pl.figure()
pl.plot(T,s21_gr_psd_re-s21,label='Amplitude readout')
pl.plot(T,s21_gr_psd_im-s21,label='Phase readout')
pl.plot(T,s21_amp-s21,color='black',label='Cold amp noise (%.1fK)'%Tamp.m)
pl.xlabel('Temperature (K)')
pl.ylabel('PSD (dBc/Hz)')
pl.title('KID noise vs temperature at 1pW bias')
pl.grid()
pl.legend(loc='upper right')
pl.ylim(ymin=-150,ymax=-60)
pl.show()


