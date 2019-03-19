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

def real_of_complex(z):
	return np.hstack([z.real, z.imag])

def complex_of_real(y):
	r, i = np.split(y, 2)
	return r + 1j*i

def admittance_model(x, C, L, R):
	w = 2*pi*x
	return real_of_complex(1/(R + 1j*w*L + 1./(1j*w*C)))

def tline_model(x, Zc, w0, eps):
	w = 2*pi*x
	return np.abs(-1/Zc/(eps*np.cos(w/w0) + 1j*np.sin(w/w0)))

def chisq(theta, x, y):
	return np.sum((y - inductance_model(x, *theta))**2)

def get_S21(Y):
	Y0 = 1/Z0
	DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
	return -2 * Y['Y21']*Y0/DY

def get_cap_params(fn):
	print (fn)
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] * MHz
	nY21 = -Ydata['Y21']
	Y2gnd = np.imag(Ydata['Y11'] + Ydata['Y21'])
	nY21 = nY21[f/MHz < 800]
	f = f[f/MHz < 800]

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

	wpeak = 2*pi*f[np.argmax(nY21.imag)]
	C_est = nY21.imag[0]/(2*pi*f[0])
	L_est = 1/(wpeak**2*C_est)
	R_est = 1e-5
	C1_est = 2.0*pF
	nY21 = nY21[f/MHz < 500]
	f = f[f/MHz < 500]

	p0 = [C_est, L_est, R_est]
	popt, pcov = curve_fit(admittance_model, f, real_of_complex(nY21), p0=p0, method='lm')
	C_fit, L_fit, R_fit = popt

	y_est = complex_of_real(admittance_model(f, *p0))

	y_fit = complex_of_real(admittance_model(f, *popt))
	print (popt)
	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, np.abs(nY21), 'b', label='Simulation')
	ax.plot(f/MHz, np.abs(y_fit), 'r-',
				label="Fit: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_fit/pF,
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
	ax.scatter(f/MHz, np.abs(nY21) - np.abs(y_fit),
				label="Guess: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_fit/pF,
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
	fn = datadir + "Waffle_cap_329MHz.csv"
	C, L, R = get_cap_params(fn)
	print (C/pF, L/nH, R)
