#! /usr/bin/env python3
"""
Script to analyze the coupling between the GP cutouts and the transmission line
for the optical TKID devices.

01/03/2020

"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi
from scipy.optimize import curve_fit, minimize
from scipy.constants import c
import os
import glob
import pdb

nH = 1e-9
pF = 1e-12
MHz = 1e6
GHz = 1e9
Z0 = 50
mm = 1e-3
um = 1e-6

l = 8000*um

datadir = '../numerical_sims/'
plotdir = 'GPcoup_fig/'

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
	w = 2*pi*x/GHz
	return real_of_complex(-1/Zc/(eps*np.cos(w/w0) + 1j*np.sin(w/w0)))

def chisq(theta, x, y):
	return np.sum((y - inductance_model(x, *theta))**2)

def get_S21(Y):
	Y0 = 1/Z0
	DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
	return -2 * Y['Y21']*Y0/DY

def get_tline_params(fn):
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] * MHz
	Y21 = Ydata['Y21']

	#fig, ax =plt.subplots(figsize=(10,10))
	#ax.plot(f/MHz, np.abs(Y21), 'b')
	#ax.plot(f/MHz, Ydata['Y21'].real, 'r', label='Y21 real')
	#ax.plot(f/MHz, Ydata['Y21'].imag, 'b', label='Y21 imag')
	#ax.plot(f/MHz, Ydata['Y11'].real, 'r--', label='Y11 real')
	#ax.plot(f/MHz, Ydata['Y11'].imag, 'b--', label='Y11 imag')
	#ax.set_xlabel('Frequency [MHz]')
	#ax.set_ylabel('|-Y21| [1/Ohms]')
	#ax.grid()
	#ax.axis('tight')
	#ax.legend(loc='best')
	#plt.show()
	#exit()
	#plt.savefig(plotdir + savename + "Y21_full.png")

	er_est = 9.2
	w0_est = c/l/er_est**0.5/GHz

	p0 = [Z0, w0_est, 0]
	popt, pcov = curve_fit(tline_model, f, real_of_complex(Y21), p0=p0, method='lm')
	Zc_fit, w0_fit, alpha_fit = popt

	y_est = complex_of_real(tline_model(f, *p0))

	y_fit = complex_of_real(tline_model(f, *popt))
	print (popt)
	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, np.abs(Y21), 'b', label='Simulation')
	ax.plot(f/MHz, np.abs(y_fit), 'r-',
				label="Fit: Z0 = %1.3fOhms \nw0 = %1.3f GHz \nalpha=%1.3e "%(Zc_fit,
					w0_fit, alpha_fit))
	#ax.plot(f/MHz, y_est, 'r-',
	#			label="Guess: Z0 = %1.3fOhms w0 = %1.3e s^-1 alpha=%1.3e "%(Z0,
	#				w0_est, 0))
	ax.legend(loc='upper right')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|-Y21| [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y21.png")

	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(f/MHz, np.abs(Y21 - y_fit),
				label="Fit: Z0 = %1.3fOhms\nw0 = %1.3f GHz\nalpha=%1.3e "%(Zc_fit,
					w0_fit, alpha_fit))
	ax.legend(loc='upper right')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Fit Residuals [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y21_residuals.png")
	plt.close('all')

	plt.close('all')
	return Zc_fit, w0_fit, alpha_fit

if __name__=="__main__":
	filenames = glob.glob(datadir + "50Ohm_feedline_with_capacitor_cutout_*umborder.csv")
	filenames.sort(key=lambda x:int(x.split('/')[-1].split('.')[0].split('_')[-1][:-8]))
	gp_margin = np.array([50.,100.,150,200,250,300])

	Zcs = np.zeros_like(gp_margin)
	w0s = np.zeros_like(gp_margin)
	alphas = np.zeros_like(gp_margin)
	for i, fn in enumerate(filenames):
		Zc, w0, alpha = get_tline_params(fn)
		Zcs[i] = Zc
		w0s[i] = w0
		alphas[i] = alpha
	er_eff = (c/(w0s*GHz*l))**2
	att = np.abs(alphas/l)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(gp_margin, Zcs, 'bo')
	ax.set_xlabel('GP margin [um]')
	ax.set_ylabel('Zc [Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + "Zc_vs_gpmargin.png")
	plt.close()

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(gp_margin, er_eff, 'bo')
	ax.set_xlabel('GP margin [um]')
	ax.set_ylabel('epsilon_eff')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + "epseff_vs_gpmargin.png")
	plt.close()

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(gp_margin, att*1e6, 'bo')
	ax.set_xlabel('GP margin [um]')
	ax.set_ylabel('Attenuation [1/m] (x1e-6)')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + "attenuation_vs_gpmargin.png")
	plt.close()

