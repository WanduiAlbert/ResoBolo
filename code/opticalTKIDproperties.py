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

nH = 1e-9
pF = 1e-12
MHz = 1e6

datadir = '../numerical_sims/'

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
		tup = lipt(zip(datapet[:, 0], p11, p12, p13, p21, p22, p23, p31, p32, p33))
		return np.array(tup, dtype=dtypes)

def admittance_model(x, C, L, C2):
	w = 2*pi*x
	return np.abs(1j*w*C2 + 1./(1./(1j*w*C) + 1j*w*L + 1e-5))

def chisq(theta, x, y):
	return np.sum((y - admittance_model(x, *theta))**2)

def get_cap_params(fn):
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] * MHz
	nY21 = -Ydata['Y21']

	wpeak = 2*pi*f[np.argmax(nY21.imag)]
	C_est = nY21.imag[0]/(2*pi*f[0])
	L_est = 1/(wpeak**2*C_est)
	C2_est = 0
	R_est = (1./nY21[0]).real

	mask = (f < (wpeak/2/pi - 100*MHz))
	#mask = (f > 0)
	p0 = [C_est, L_est, C2_est]
	#popt, pcov = curve_fit(admittance_model, f[mask], nY21[mask], method='lm')
	result = minimize(chisq, p0, args=(f[mask], np.abs(nY21[mask])),
			method='Nelder-Mead')
	C_fit, L_fit, C2_fit = result["x"]
	#C_fit, L_fit, R_fit = result["x"]
	#R_fit = np.abs(R_fit)
	#print (C_fit/pF, L_fit/nH, R_fit)

	y_est = admittance_model(f, C_est, L_est, C2_est)
	#print (wpeak/2/pi/MHz)
	#print (C_est/pF, L_est/nH, R_est)
	#fig, ax =plt.subplots(figsize=(10,10))
	#ax.semilogy(f/MHz, np.abs(nY21), 'b', label='Simulation')
	#ax.semilogy(f/MHz, y_est, 'k--',
	#		label="est C = %1.3fpF L = %1.3fnH R=%1.3f Ohms"%(C_est/pF,
	#			L_est/nH, R_est))
	#ax.legend()
	#ax.set_xlabel('Frequency [MHz]')
	#ax.set_ylabel('|-Y21| [1/Ohms]')
	#ax.grid()
	#ax.axis('tight')
	#plt.savefig(savename + "Y21.png")

	#exit()
	y_fit = admittance_model(f, C_fit, L_fit, C2_fit)
	fig, ax =plt.subplots(figsize=(10,10))
	ax.semilogy(f/MHz, np.abs(nY21), 'b', label='Simulation')
	ax.semilogy(f/MHz, y_fit, 'k--',
			label="Fit C = %1.3fpF L = %1.3fnH C2=%1.3f Ohms"%(C_fit/pF,
				L_fit/nH, C2_fit/pF))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|-Y21| [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(savename + "Y21.png")

	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f[mask]/MHz, np.abs(nY21)[mask] - y_fit[mask],
			'b', label="Fit C = %1.3fpF L = %1.3fnH C2=%1.3f Ohms"%(C_fit/pF,
				L_fit/nH, C2_fit/pF))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('-Y21 [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(savename + "Y21_residuals.png")

	return C_fit, L_fit, C2_fit

if __name__=="__main__":
	L_al = 10 * nH
	filenames = glob.glob(datadir + "*_with_boundary.csv")
	expected_freqs = np.array([300, 305, 310, 315]) * MHz
	caps = []
	inds = []
	caps2 =  []
	for fn in filenames:
		print (fn)
		C, L, R = get_cap_params(fn)
		caps.append(C)
		inds.append(L)
		caps2.append(R)
	caps = np.array(caps)
	inds = np.array(inds)
	caps2 = np.array(caps2)

	print (caps/pF)
	print (inds/nH)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, caps/pF, 'k-x')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig('cap_vs_design_freq.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, inds/nH, 'k-x')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	plt.savefig('ind_vs_design_freq.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, caps2, 'k-x')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'caps2istance [Ohms]')
	plt.savefig('caps2_vs_design_freq.png')

	plt.close('all')
	L_tot = inds + L_al
	freqs = 1./np.sqrt(L_tot * caps)/2/pi
	print (freqs/MHz)


	p = np.polyfit(freqs/MHz, expected_freqs/MHz, 1)
	print (p)
	m, b = p
	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(freqs/MHz, expected_freqs/MHz)
	ax.plot(freqs/MHz, np.polyval(p, freqs/MHz), 'k',
			label='m=%1.3f\nb=%1.3f'%(m,b))
	ax.grid(which='both')
	ax.set_xlabel(r'Actual Frequency [MHz]')
	ax.set_ylabel(r'Design Frequency [MHz]')
	ax.legend(loc='upper left')
	plt.savefig('design_vs_actual_freq.png')
