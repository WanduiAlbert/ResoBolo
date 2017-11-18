
from math import pi
import pint
import numpy as np
from scipy import special


u = pint.UnitRegistry()

h = u.planck_constant
k = u.boltzmann_constant


def dpdt_fun(fs,bw,t):
	x =  h * fs / k
	eu = np.exp(x/t)
	dpdt = bw * h * fs * eu * x / (t * (eu-1.))**2
	dpdt = dpdt.to(u.picowatt / u.kelvin)
	return dpdt

class OperatingPoint():
	''' DC operating point of a resonator bolometer '''
	def __init__(self):
		# Material properties
		Tc = 2.0 * u.kelvin
		#N0 = 1.7e10 / (u.electron_volt * u.micrometer**3)
		N0 = 0.89e10 / (u.electron_volt * u.micrometer**3)
		taumax = 0.0001 * u.second
		nstar = 100 / (u.micrometer**3)

		#Tc = 1.8 * u.kelvin
		#N0 = 8.7e9 / (u.electron_volt * u.micrometer**3)
		#taumax = 0.000114 * u.second
		#nstar = 175 / (u.micrometer**3)

		delta0 = (1.76 * k * Tc).to('joule')
		print 'delta0: ',delta0.to('eV')

		#rho = 2.5*u.microohm * u.centimeter
		rho = 100*u.microohm * u.centimeter
		t = 0.05*u.micrometer
		w = 1.0*u.micrometer
		A = t*w
		Z0 = 50.*u.ohm
		
		p = 2 * pi * N0*delta0**3 * A**2 * Z0 / (h * rho)
		print "critical readout power",p.to('watt')

		island_height = 150*u.micrometer
		island_width = 300*u.micrometer
		filling_fraction = 0.4
		area = island_height * island_width
		l = filling_fraction * area / w
		V = l * t * w
		#Rs = 0.5*u.ohm/u.micrometer
		Rs = (rho /  ( t * w )).to('ohm/micrometer')
		Rsheet = (Rs * w).to('ohm')
		print "Rsheet: ",Rsheet
		print "Rs: ",Rs
		Ls = (h * Rs / (2*pi*pi * delta0)).to('picohenry / micrometer')
		Lk = l * Ls
		L = Lk + 8.0*u.nanohenry
		L = L.to('nanohenry')
		#L = 200*u.nanohenry
		print "L: ",L

		# Optical properties
		nu = 220 * u.gigahertz
		bw = 0.25
		eta = 0.3					# Front of dewar to dissipation resistor efficiency
		Tsky = 31.* u.kelvin		# Total optical load referenced to front of dewar (rj)
		#Tsky = 43. * u.kelvin
		rjtop = (k * eta * nu * bw).to('picowatt/kelvin')
		Po = (Tsky * rjtop).to('picowatt')
		print "Po: ",Po

		#Psat = 1.5*Po
		#Pr = Psat - Po				# Readout power
		Pr = 1.0*u.picowatt
		Psat = Pr + Po
		Pe = 0.5 * Pr				# Fraction of readout power dissipated in bolometer

		Pt = Po + Pe				# Total power dissipated on bolometer
		print "Psat: ",Pt

		To = 0.50 * u.kelvin		# Desired operating temperature of bolometer
		Tb = 0.25 * u.kelvin		# Bolometer bath temperature
		thermal_beta = 1.6					# Thermal conductivity exponent of legs

		kc = Pt / (To**(thermal_beta+1.) - Tb**(thermal_beta+1.))

		G = kc * (thermal_beta + 1.) * To**thermal_beta		# bolometer G at operating temperature
		Gb = kc * (thermal_beta+1.) * Tb**thermal_beta
		print "Gc: ",G

		Z0 = 50 * u.ohm				# Readout line impedance
		#fro = 0.3 * u.gigahertz		# KID resonator frequency
		#C = 1./((2*pi*fro)**2 * L)	# KID capacitance
		#C = C.to('picofarad')

		C = 16.9*u.picofarad
		fro = (1.0/(2.*pi*np.sqrt(L*C))).to('megahertz')
		print "fro: ",fro

		print "C: ",C
		alpha = Lk/L					# Kinetic inductance fraction
		print "alpha: ",alpha.to('')
		Tn = 3.0 * u.kelvin			# Noise temperature of readout amplifier
		Qiloss = 3e5					# Intrinsic resonator Q

		# Thermal quasiparticle density
		nqp = (2 * N0 * np.sqrt(2. * pi * k * To * delta0) * np.exp(-delta0 / (k * To))).to(1./u.micrometer**3)
		print "nqp: ",nqp
		# Change in nqp density due to operating temperature change
		dnqpdt = nqp * (k * To + 2 * delta0) / (2 * k * To**2)

		# change in nqp from optical power
		dnqpdp = dnqpdt / G
		# change in nqp from sky temperature
		dnqpdts = dnqpdp * rjtop

		tauqp = taumax / (1. + nqp / nstar)
		print "tauqp: ",tauqp.to('microseconds')

		gammarec = 0.5 * nqp * V * (1./taumax + 1./tauqp)
		print "gammarec: ",gammarec*1e-9

		hf2kt = (h * fro / (2 * k * To)).to('')
		y = np.sqrt(2*delta0 / (pi * k * To))
		S1 = ((2./pi) * np.sqrt(2 * delta0 / (pi * k * To)) * np.sinh(hf2kt) * special.k0(hf2kt)).to('')
		S2 = (1.0 +  np.sqrt(2 * delta0 / (pi * k * To)) * np.exp(-hf2kt) * special.i0(hf2kt)).to('')
		beta = S2/S1
		self.beta = beta
		print "beta: ",beta

		# Fractional shift in resonant frequency from quasiparticles
		x = (alpha * S2 * nqp / (4 * N0 * delta0)).to('')
		dxdnqp = (alpha * S2 / (4 * N0 * delta0)).to('micrometer**3')
		# Responsivity - frequency shift due to change in sky temperature, ignoring change in S2
		dxdts = (dxdnqp * dnqpdts).to('1/K')

		# Absorption due to quasiparticles
		filmgamma = 1.0
		Qiqp = (2 * N0 * delta0 / (alpha * nqp * S1 * filmgamma)).to('')
		Qi = 1. / (1./Qiqp + 1./Qiloss)

		# Optimal Qc
		Qc = Qi
		Cc = np.sqrt(C / (pi *fro * Qc * Z0)).to('picofarad')
		# Fixed Cc
		#Cc = 0.2 * u.picofarad
		Qc = (C / (pi * fro * Cc**2 * Z0)).to('')	# Coupling Q

		# Effective Q of resonator
		Qr = 1. / (1./Qi + 1./Qc)

		print "Qi: ",Qi
		print "Qr: ",Qr
		print "Cc: ",Cc

		s21 = 1 - Qr/Qc
		ds21dx =  0.5 * 0.5 * Qi

		# Amplifier noise
		efamp = (np.sqrt(4*k*Tn/Pr) * (Qc/Qr**2) * fro).to('hertz/hertz**0.5')
		resp = dxdts * fro
		print "responsivity: ",resp.to('kilohertz/kelvin')
		netamp = 0.5**0.5*(efamp / resp).to('microkelvin * second**0.5')

		dpdtcmb = dpdt_fun(nu,nu*bw,2.725*u.kelvin)
		dpdtrj = dpdt_fun(nu,nu*bw,300*u.kelvin)
		cmbtorj = dpdtrj/dpdtcmb

		nepshot2 = (2 * h * nu * Po).to(u.attowatt**2 / u.hertz)
		nepshot = np.sqrt(nepshot2)
		nepwave2 = (2*Po*Po/(nu * bw)).to(u.attowatt**2 / u.hertz)
		nepwave = np.sqrt(nepwave2)
		nepphot = np.sqrt(nepshot2 + nepwave2)
		netphot = 0.5**0.5*(nepphot / rjtop).to('microkelvin * second**0.5')

		tau_eff = 1. / (1./taumax + 1./tauqp)
		kappa = 0.5 * (1 + 3.5 * Tc/To)

		x = Pr/Po
		xcr = (V * nstar * delta0/(Po*taumax))*(1 + nqp/nstar)**2
		psi = (1+x/xcr)*(1+x)**2	

		gammaG = (To/Pt)*G
		Nqp = nqp * V
		nepgr2 = (2*gammaG*(1+Pr/Po)/kappa)**2*(1+tauqp/taumax)*(tauqp/Nqp)*Po**2
		nepgr = (nepgr2**0.5).to('attowatt/hertz**0.5')

		netgr = 0.5**0.5*(nepgr/rjtop).to('microkelvin*second**0.5')

		nepamp2 = (2*gammaG/(kappa*beta))**2 * (psi/x)*(k*Tn/(h*nu)) * (2*h*nu*Po)
		nepamp = np.sqrt(nepamp2).to('attowatt/hertz**0.5')
		netamp = 0.5**0.5*(nepamp/rjtop).to('microkelvin*second**0.5')

		#gamma = 0.5
		gammaphonon = ((1 + thermal_beta) / (2*thermal_beta + 3)) * ((Tb/To)**(3.+2.*thermal_beta)-1.)/((Tb/To)**(thermal_beta+1.)-1)
		print "gamma phonon: ",gammaphonon
		nepphon2 = (4 * k * To**2 * gammaphonon * G).to('attowatt**2 / hertz')
		nepphon = nepphon2 ** 0.5
		netphon = 0.5**0.5*(nepphon / rjtop).to('microkelvin * second**0.5')


		self.Qc = Qc
		self.Qiqp = Qiqp
		self.Qr = Qr
		self.x = x.to('')
		self.fro = fro.to('hertz')
		self.efamp = efamp
		self.cmbtorj = cmbtorj
		self.netphot = netphot
		self.netphon = netphon
		self.netgr = netgr
		self.nqp = nqp
		self.dnqpdt = dnqpdt
		self.dxdts = dxdts
		self.rjtop = rjtop
		self.nepamp = nepamp
		self.netamp = netamp

o = OperatingPoint()

print '###'
print 'device parameters'
print "cmbtorj: ",o.cmbtorj
print "Qiqp: ",o.Qiqp
print "Qr: ",o.Qr

resp = o.dxdts * o.fro
print "responsivity: ",resp
print 'nqp: ',o.nqp

netphot = o.netphot
netamp = 0.5**0.5*(o.efamp / resp).to('microkelvin * second**0.5')
nepamp = (netamp * 2**0.5 * o.rjtop).to('W/Hz^0.5')

netphon = o.netphon
netgr = o.netgr

net = (netphot**2 + netamp**2 + netphon**2 + netgr**2)**0.5
netcmb = net * o.cmbtorj

s = o.rjtop * 2**0.5
nepgr = (netgr*s).to('aW/Hz^0.5')
nepphon = (netphon*s).to('aW/Hz^0.5')
nepamp = (netamp*s).to('aW/Hz^0.5')
nepphot = (netphot*s).to('aW/Hz^0.5')
nep = (net*s).to('aW/Hz^0.5')

print ''
print 'NEP'
print 'gr:\t',nepgr
print 'phonon:\t',nepphon
print 'amp:\t',nepamp
print 'amp:\t',o.nepamp
print 'photon:\t',nepphot
print 'total:\t',nep


print ''
print 'NET'
print 'gr:\t',netgr
print 'phonon:\t',netphon
print 'amp:\t',netamp
print 'amp:\t',o.netamp
print 'photon:\t',netphot
print 'netrj:\t',net
print ''
print 'netcmb:\t',netcmb
print ''



