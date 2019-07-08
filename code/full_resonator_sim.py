#! /usr/bin/env python3
"""
06/07/2019
--------------------------------------------------------------------------------
Full resonator sim.

I generated a netlist of the optical TKID resonator without the filter in place
and this script is meant to explore the predicted resonances and other
properties.

"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi
from scipy.optimize import curve_fit, minimize
from scipy import interpolate
import os
import glob
import pdb

nH = 1e-9
pF = 1e-12
MHz = 1e6
Z0 = 50

datadir = '../numerical_sims/'
plotdir = ''

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

def Y_to_Z(Y):
	Y11 = Y['Y11']
	Y12 = Y['Y12']
	Y21 = Y['Y21']
	Y22 = Y['Y22']
	dtypes = np.dtype([("frequency", np.float64), ("Z11", np.complex128),\
		("Z12", np.complex128), ("Z21", np.complex128), ("Z22", np.complex128)])
	detY = Y11*Y22 - Y12*Y21
	Z11 = Y22/detY
	Z12 = -Y12/detY
	Z21 = -Y21/detY
	Z22 = Y11/detY

	tup = list(zip(Y["frequency"], Z11, Z12, Z21, Z22))
	return np.array(tup, dtype=dtypes)

def admittance_model(x, C, L, R):
	w = 2*pi*x
	return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
	#return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

def tline_model(x, Zc, w0, eps):
	w = 2*pi*x
	return np.abs(-1/Zc/(eps*np.cos(w/w0) - 1j*np.sin(w/w0)))

# Including a capacitive coupling term on the line
def tline_model_x(x, Zc, w0, eps, L):
	w = 2*pi*x
	#Series capacitance model
	#C1 = L
	#num = w**2*Zc*C1**2
	#denom_a = -2j*C1*Zc*w + eps*(1 - C1**2*Zc**2*w**2)
	#denom_b = 2*C1*Zc*eps*w - 1j*(1 - C1**2*Zc**2*w**2)
	#Series inductance model
	num = -Zc
	denom_a = Zc**2*eps + 2j*L*Zc*w - L**2*eps*w**2
	denom_b = 1j*(Zc**2 - L**2*w**2) - 2*L*Zc*eps*w
	print (Zc, w0, eps, L)
	return real_of_complex(num/(denom_a*np.cos(w/w0) + denom_b*np.sin(w/w0)))

def chisq(theta, x, y):
	#return np.sum((y - admittance_model(x, *theta))**2)
	return np.sum((y - tline_model(x, *theta))**2)

def get_S21(Y):
	Y0 = 1/Z0
	DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
	return -2 * Y['Y21']*Y0/DY

def real_of_complex(z):
	return np.hstack([z.real, z.imag])

def complex_of_real(y):
	r, i = np.split(y, 2)
	return r + 1j*i

def capacitor_model(f, C):
	return -1/(2*pi*f*C)

def overall_series_model(f, C, L):
	w = 2*pi*f
	return w*L - 1./(w*C)

def overall_parallel_model(f, C):
	w = 2*pi*f
	return w*C #- 1./(w*L)

def resonances_in_Z(f, Zdata):
	Zin_1 = Zdata['Z11'] - (Zdata['Z21']*Zdata['Z12'])/(Zdata['Z22'] + Z0)
	Zin_2 = Zdata['Z22'] - (Zdata['Z21']*Zdata['Z12'])/(Zdata['Z11'] + Z0)

	# Need to mask out all the places with possible resonances
	mask1 = (f/MHz > 230) &  (f/MHz < 330)
	mask2 = (f/MHz > 1100) & (f/MHz < 1300)
	mask3 = (f/MHz > 2000) & (f/MHz < 2300)
	mask4 = (f/MHz > 2700) & (f/MHz < 3000)
	mask5 = (f/MHz > 3400) & (f/MHz < 3800)

	mask = ~(mask1 | mask2 | mask3 | mask4 | mask5)

	# I want to fit the masked Zin,2 to a capacitive impedance and then subtract
	# that component out.
	p0 = [0.2*pF, 0.01*nH]
	popt, pcov = curve_fit(overall_series_model, f[mask], Zin_2[mask].imag, p0,
			method='lm')
	C, L = popt
	print (C/pF)
	print (L/nH)
	Zfit = overall_series_model(f, *popt)

	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, Zin_2.real, 'r', label='real')
	ax.plot(f/MHz, Zin_2.imag, 'b', label='imag')
	ax.plot(f/MHz, Zfit, 'k--', label='imag')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Z in,2 [Ohms]')
	ax.grid()
	ax.axis('tight')
	ax.legend(loc='best')
	plt.show()
	#exit()
	#plt.savefig(plotdir + savename + "Z21_fullreso.png")
	#plt.close()

	# I can use interpolate to find the roots and derivatives of Zin,2
	spl = interpolate.splrep(f/MHz, Zin_2.imag - Zfit, s=0)
	spl_real = interpolate.splrep(f/MHz, Zin_2.real, s=0)
	roots = interpolate.sproot(spl, mest=11)
	resofreqs = roots#[1::2]
	ks = interpolate.splev(resofreqs, spl, der=1) # find the slopes
	Rs = interpolate.splev(resofreqs, spl_real, der=0) # real impedances

	print (ks)
	Qs = ks*resofreqs*MHz/(2*Rs)
	print (resofreqs)
	print (Qs)

	return resofreqs, Qs

def resonances_in_Y(f, Ydata):
	Y0 = 1./Z0
	Yin_1 = Ydata['Y11'] - (Ydata['Y21']*Ydata['Y12'])/(Ydata['Y22'] + Y0)
	Yin_2 = Ydata['Y22'] - (Ydata['Y21']*Ydata['Y12'])/(Ydata['Y11'] + Y0)

	# Need to mask out all the places with possible resonances
	mask1 = (f/MHz > 230) &  (f/MHz < 330)
	mask2 = (f/MHz > 1100) & (f/MHz < 1300)
	mask3 = (f/MHz > 2000) & (f/MHz < 2300)
	mask4 = (f/MHz > 2700) & (f/MHz < 3000)
	#mask5 = (f/MHz > 3400) & (f/MHz < 3800)
	mask5 = (f/MHz > 3400)

	mask = ~(mask1 | mask2 | mask3 | mask4 | mask5)

	# I want to fit the masked Yin,2 to a capacitive impedance and then subtract
	# that component out.
	#p0 = 0.2*pF#, 0.01*nH]
	#popt, pcov = curve_fit(overall_parallel_model, f[mask], Yin_2[mask].imag, p0,
	#		method='lm')
	#C, = popt
	#print (C/pF)
	##print (L/nH)
	#Yfit = overall_parallel_model(f, *popt)

	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, Yin_2.real, 'r', label='real')
	ax.plot(f/MHz, Yin_2.imag, 'b', label='imag')
	#ax.plot(f/MHz, Yfit, 'k--', label='imag')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Y in,2 [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	ax.legend(loc='best')
	plt.show()
	exit()
	#plt.savefig(plotdir + savename + "Y21_fullreso.png")
	#plt.close()

	# I can use interpolate to find the roots and derivatives of Yin,2
	spl = interpolate.splrep(f/MHz, Yin_2.imag - Yfit, s=0)
	spl_real = interpolate.splrep(f/MHz, Yin_2.real, s=0)
	roots = interpolate.sproot(spl, mest=11)
	resofreqs = roots#[1::2]
	ks = interpolate.splev(resofreqs, spl, der=1) # find the slopes
	Rs = interpolate.splev(resofreqs, spl_real, der=0) # real impedances

	print (ks)
	Qs = ks*resofreqs*MHz/(2*Rs)
	print (resofreqs)
	print (Qs)

	return resofreqs, Qs

def get_reso_params(fn):
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] * MHz
	nY21 = Ydata['Y21']
	Y2gnd = np.imag(Ydata['Y11'] + Ydata['Y21'])

	# Construct the Z parameters from the Y parameters
	Zdata = Y_to_Z(Ydata)
	y_frs, y_Qs = resonances_in_Y(f, Ydata)
	#z_frs, z_Qs = resonances_in_Z(f, Zdata)

	return z_frs, z_Qs


if __name__=="__main__":
	fn = datadir + 'capacitor_and_inductor_sim_with_coupcaps_noloss.csv'
	params = get_reso_params(fn)
