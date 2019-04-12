#! /usr/bin/env python3
"""
Script to analyze a set of simulations of capacitors for the optically coupled
TKID devices. The code takes in sonnet simulation results expressed in the Y
parameters and computes the Capacitance and parasitic Inductance for the
different geometries.


"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi
from scipy.optimize import curve_fit, minimize
import os
import glob
import pdb

nH = 1e-9
pF = 1e-12
MHz = 1e6
Z0 = 50

datadir = '../numerical_sims/'
plotdir = 'waffle_cap_fig/'

def load_data(fn, nports=2, paramtype='Y'):
	dataset = np.loadtxt(fn, skiprows=9, delimiter=',')
	p = paramtype
	if nports==2:
		dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
			(p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
		p11 = dataset[:, 1] + 1j*dataset[:, 2]
		p12 = dataset[:, 3] + 1j*dataset[:, 4]
		p21 = dataset[:, 5] + 1j*dataset[:, 6]
		p22 = dataset[:, 7] + 1j*dataset[:, 8]
		tup = list(zip(dataset[:, 0], p11, p12, p21, p22))
		return np.array(tup, dtype=dtypes)
	elif nports==3:
		dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
			(p + "12", np.complex128), (p + "13", np.complex128), (p + "21", np.complex128),\
			(p + "22", np.complex128), (p + "23", np.complex128), (p +"31", np.complex128),\
			(p + "32", np.complex128), (p  +"33", np.complex128)])
		p11 = dataset[:, 1] + 1j*dataset[:, 2]
		p12 = dataset[:, 3] + 1j*dataset[:, 4]
		p13 = dataset[:, 5] + 1j*dataset[:, 6]
		p21 = dataset[:, 7] + 1j*dataset[:, 8]
		p22 = dataset[:, 9] + 1j*dataset[:,10]
		p23 = dataset[:,11] + 1j*dataset[:,12]
		p31 = dataset[:,13] + 1j*dataset[:,14]
		p32 = dataset[:,15] + 1j*dataset[:,16]
		p33 = dataset[:,17] + 1j*dataset[:,18]
		tup = lipt(zip(dataset[:, 0], p11, p12, p13, p21, p22, p23, p31, p32, p33))
		return np.array(tup, dtype=dtypes)

def admittance_model(x, C, L, R):
	w = 2*pi*x
	return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
	#return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

def tline_model(x, Zc, w0, eps):
	w = 2*pi*x
	return np.abs(-1/Zc/(eps*np.cos(w/w0) + 1j*np.sin(w/w0)))

def chisq(theta, x, y):
	return np.sum((y - admittance_model(x, *theta))**2)

def get_S21(Y):
	Y0 = 1/Z0
	DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
	return -2 * Y['Y21']*Y0/DY

def get_cap_params(fn):
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] * MHz
	nY21 = -Ydata['Y21']
	Y2gnd = np.imag(Ydata['Y11'] + Ydata['Y21'])

	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, np.abs(nY21), 'b')
	#ax.plot(f/MHz, Ydata['Y21'].real, 'r', label='Y21 real')
	#ax.plot(f/MHz, Ydata['Y21'].imag, 'b', label='Y21 imag')
	#ax.plot(f/MHz, Ydata['Y11'].real, 'r--', label='Y11 real')
	#ax.plot(f/MHz, Ydata['Y11'].imag, 'b--', label='Y11 imag')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|-Y21| [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	#ax.legend(loc='best')
	#plt.show()
	#exit()
	plt.savefig(plotdir + savename + "Y21_full.png")

	#S21 = np.abs(get_S21(Ydata))
	#plt.plot(f/1e6, S21)
	#plt.show()
	#exit()
	wpeak = 2*pi*f[np.argmax(nY21.imag)]
	C_est = nY21.imag[0]/(2*pi*f[0])
	L_est = 1/(wpeak**2*C_est)
	R_est = 1e-5
	C1_est = 2.0*pF
	nY21 = nY21[f/MHz < 500]
	f = f[f/MHz < 500]


	p0 = [C_est, L_est, R_est]
	#pdb.set_trace()
	popt, pcov = curve_fit(admittance_model, f, np.log(np.abs(nY21)), p0=p0, method='lm')
	#popt, pcov = curve_fit(tline_model, f, np.log(np.abs(nY21)), p0=p0, method='lm')
	#result = minimize(chisq, p0, args=(f, np.abs(nY21)),
	#		method='Nelder-Mead')
	C_fit, L_fit, R_fit = popt
	#C_fit, L_fit, R_fit, C1_fit = result['x']

	y_est = np.exp(admittance_model(f, C_est, L_est, R_est))

	y_fit = np.exp(admittance_model(f, C_fit, L_fit, R_fit))
	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, np.abs(nY21), 'b', label='Simulation')
	#ax.plot(f/MHz, y_est, 'g', label='Guess')
	ax.plot(f/MHz, y_fit, 'k--',
				label="Fit C = %1.3fpF L = %1.3fnH R=%1.3e Ohms "%(C_fit/pF,
					L_fit/nH, R_fit))
	#ax.plot(f/MHz, y_est, 'r-',
		#		label="Guess: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_est/pF,
		#			L_est/nH, R_est))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|-Y21| [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y21.png")

	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(f/MHz, np.abs(nY21) - y_fit,
				label="Fit C = %1.3fpF L = %1.3fnH R=%1.3e Ohms "%(C_fit/pF,
					L_fit/nH, R_fit))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Fit Residuals [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y21_residuals.png")
	plt.close('all')

	plt.close('all')
	return C_fit, L_fit, R_fit

if __name__=="__main__":
	L_al_g =  5.998* nH
	filenames = glob.glob(datadir + "Waffle_cap_*MHz.csv")
	expected_freqs = np.array([306.21, 321.64, 336.67, 352.06, 367.00]) * MHz
	filenames.sort(key=lambda x:int(x.split('/')[-1].split('.')[0].split('_')[-1][:3]))
	Nfingers = np.array([336.5,  365.5,  399.5,  437.5, 482.5])
	caps = []
	inds = []
	Rs =  []
	for fn in filenames:
		print (fn)
		C, L, R = get_cap_params(fn)
		caps.append(C)
		inds.append(L)
		Rs.append(np.abs(R))
	caps = np.array(caps)
	inds = np.array(inds)
	Rs = np.array(Rs)

	print (caps/pF)
	print (inds/nH)

	Lk = 4.059 *nH
	L_al = Lk + L_al_g
	alphak = Lk/(L_al + inds)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, caps/pF, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig(plotdir + 'cap_vs_design_freq.png')

	def cap_logfit(N, A, B):
		return A*np.log(N) + B
	def cap_powfit(N, A, B):
		return A*N**0.5 + B
	def cap_harmonicfit(N, A, B):
		return A*(1. - 1./N) + B

	cap_p = np.polyfit(1./Nfingers, caps/pF, 1)
	cap_l = np.polyfit(Nfingers, caps/pF, 1)
	p0 = [cap_p[0], cap_p[1]]
	#popt, _ = curve_fit(cap_logfit, Nfingers, caps/pF, p0=p0, method='lm')
	#popt2, _ = curve_fit(cap_powfit, Nfingers, caps/pF, p0=p0, method='lm')
	#popt3, _ = curve_fit(cap_harmonicfit, Nfingers, caps/pF, p0=p0, method='lm')
	x = np.linspace(250, 550, 100)
	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers, caps/pF, 'ko')
	#ax.plot(x, np.polyval(cap_p, 1./x), 'k--',
	#	label='Harmonic Fit: C (in pF) = {0:1.3f}/N + {1:1.3f}'.format(*cap_p))
	ax.plot(x, np.polyval(cap_l, x), 'r--',
		label='Linear Fit: C (in pF) = {0:1.3f}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig(plotdir + 'cap_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, inds/nH, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	plt.savefig(plotdir + 'ind_vs_design_freq.png')

	ind_p = np.polyfit(Nfingers, inds/nH, 2)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers, inds/nH, 'ko')
	ax.plot(x, np.polyval(ind_p, x), 'k--',
		label='Quadratic Fit:a={0:1.3e} b={1:1.3f} c={2:1.3f}'.format(*ind_p))
	ax.grid(which='both')
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	plt.savefig(plotdir + 'ind_vs_nfingers.png')


	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, alphak, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'alphak')
	plt.savefig(plotdir + 'alphak_vs_design_freq.png')

	alpha_p = np.polyfit(Nfingers[:-1], alphak[:-1], 1)
	alpha_m, alpha_b = alpha_p
	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers, alphak, 'ko')
	ax.plot(x, np.polyval(alpha_p, x), 'k-',
			label="m={0:1.3e}, b={1:1.3f}".format(alpha_m, alpha_b))
	ax.legend()
	ax.grid(which='both')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'alphak')
	plt.savefig(plotdir + 'alphak_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, Rs*1e5, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Resistance (x10^-5) [Ohms]')
	plt.savefig(plotdir + 'Rs_vs_design_freq.png')

	plt.close('all')
	L_tot = inds + L_al
	freqs = 1./np.sqrt(L_tot * caps)/2/pi
	print (freqs/MHz)

	ap = np.polyfit(freqs/MHz, alphak, 2)
	aa, ab,ac = ap
	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(freqs/MHz, alphak, 'ko')
	ax.plot(freqs/MHz, np.polyval(ap, freqs/MHz), 'k-',
			label='a={0:1.3e}, b={1:1.3e}, c={2:1.3e}'.format(*ap))
	ax.grid(which='both')
	ax.legend(loc='upper left')
	ax.set_xlabel(r'True Frequency [MHz]')
	ax.set_ylabel(r'alphak')
	plt.savefig(plotdir + 'alphak_vs_true_freq.png')

	p = np.polyfit(expected_freqs/MHz, freqs/MHz, 1)
	print (p)
	m, b = p
	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(expected_freqs/MHz, freqs/MHz)
	ax.plot(expected_freqs/MHz, np.polyval(p, expected_freqs/MHz), 'k',
			label='m=%1.3f\nb=%1.3f'%(m,b))
	ax.grid(which='both')
	ax.set_ylabel(r'True Frequency [MHz]')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.legend(loc='upper left')
	plt.savefig(plotdir + 'design_vs_actual_freq.png')
