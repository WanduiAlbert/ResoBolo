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
from scipy.constants import c
import os
import glob
import pdb

nH = 1e-9
pF = 1e-12
MHz = 1e6
Z0 = 50

datadir = '../numerical_sims/'
plotdir = 'coupcap_fig/'

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
		tup = list(zip(dataset[:, 0], p11, p12, p13, p21, p22, p23, p31, p32, p33))
		return np.array(tup, dtype=dtypes)

def admittance_model_toGND(x, Cp, C, L1, L2, M, R):
	w = 2*pi*x
	Cp  *= pF
	C  *= pF
	L1 *= nH
	L2 *= nH
	M  *= nH
	beta = -4*C*w + w**3*C**2*(L1+L2 - 2*M)
	a = 4 - 2*w**2*C*(L1+L2) + w**4*C**2*(L1*L2-M**2)
	alpha = w**2*C**2*R
	b = 2*w*C*R
	return w*Cp + (beta*a - alpha*b)/(a**2 + b**2)
	return w*Cp + beta/a

def admittance_model_simple(x, C, L, R):
	w = 2*pi*x
	C  *= pF
	L *= nH
	beta = -C*w
	a = 1 - w**2*C*L
	alpha = w**2*C**2*R
	b = 2*w*C*R
	return (beta*a - alpha*b)/(a**2 + b**2)


def admittance_model(x, C, L1, L2, M, R):
	w = 2*pi*x
	C  *= pF
	L1 *= nH
	L2 *= nH
	M  *= nH
	beta = -4*C*w + w**3*C**2*(L1+L2 - 2*M)
	a = 4 - 2*w**2*C*(L1+L2) + w**4*C**2*(L1*L2-M**2)
	alpha = w**2*C**2*R
	b = 2*w*C*R
	return (beta*a - alpha*b)/(a**2 + b**2)
	return beta/a


#def admittance_model(x, C, L, R):
#	w = 2*pi*x
#	return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
#	#return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

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
	Ydata = load_data(fn, nports=3)
	Y0 = 1./50
	n = 8.93**0.5
	l = 2656e-6
	f = Ydata['frequency'] * MHz
	beta = 2*pi*f*n/c
	Yload = (Y0*(Ydata['Y31'] - 1j*Y0*np.tan(beta*l))/(Y0 -
			1j*Ydata['Y31']*np.tan(beta*l))).imag
	Y31 = Ydata['Y31'].imag
	Y2gnd = np.imag(Ydata['Y11'] + Ydata['Y31'])
	Yodd = Ydata['Y11'] - Ydata['Y31']

	L_est = 0
	R_est = 1e-5
	C_est = 0.2#*pF
	#nY31 = nY31[f/MHz < 500]
	#f = f[f/MHz < 500]
	#print (C_est)
	#print (L_est)
	#exit()

	p0 = [C_est, L_est, R_est]
	popt, pcov = curve_fit(admittance_model_simple, f[f/MHz < 500], Yload[f/MHz < 500], p0=p0, method='lm')
	C_fit, L_fit, R_fit = popt
	y_est = admittance_model_simple(f, *p0)

	y_fit = admittance_model_simple(f, *popt)
	fig, ax =plt.subplots(figsize=(10,10))
	ax.plot(f/MHz, Yload, 'b', label='Simulation')
	#ax.plot(f/MHz, y_est, 'g', label='Guess')
	ax.plot(f/MHz, y_fit, 'k--',
				label="Fit C = %1.3fpF L = %1.3fnH "%(C_fit, L_fit))
	#ax.plot(f/MHz, y_est, 'r-',
	#			label="Fit C = %1.3fpF L = %1.3fnH "%(C_fit, L_fit))
	ax.legend(loc='upper left')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Yload [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y31.png")

	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(f/MHz, Yload - y_fit,
				label="Fit C = %1.3fpF L = %1.3fnH "%(C_fit, L_fit))
	ax.legend(loc='upper left')
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('Fit Residuals [1/Ohms]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y31_residuals.png")
	plt.close('all')

	return C_fit, L_fit, 0, 0



if __name__=="__main__":
	L_al = 10 * nH
	filenames = glob.glob(datadir + "coupcap_and_feedline_*.csv")
	expected_freqs = np.array([300, 305, 310, 315, 320, 325, 330, 335, 340,
			345, 350, 360, 370, 380, 390, 400, 410]) * MHz
	#filenames.sort(key=lambda x:int(x.split('/')[-1].split('.')[0].split('_')[1][:3]))
	Nfingers = 502 - np.array([0, 17, 32, 47, 61, 74, 87, 100, 111, 122, 133,
		153, 172, 189, 205, 219, 233])
	#mask = np.ones_like(Nfingers, dtype=bool)
	mask = Nfingers % 2 == 1
	Cs = []
	Ls = []
	Cp1s = []
	Cp2s = []
	for i, fn in enumerate(filenames):
		C, L, Cp1, Cp2 = get_cap_params(fn)
		Cs.append(C)
		Ls.append(L)
		Cp1s.append(Cp1)
		Cp2s.append(Cp2)
		print (fn)
		print (C)
		print (L)
	exit()
	Cs = np.array(Cs)
	L1s = np.array(L1s)
	L2s = np.array(L2s)
	Ms = np.array(Ms)
	Cp1s = np.array(Cp1s)
	Cp2s = np.array(Cp2s)



	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], Cs[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig(plotdir + 'C_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], L1s[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'L1 [nH]')
	plt.savefig(plotdir + 'L1_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], L2s[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'L2 [nH]')
	plt.savefig(plotdir + 'L2_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], a[mask], 'ko')
	ax.plot(Nfingers[mask], b[mask], 'bo')
	ax.plot(Nfingers[mask], c[mask], 'ro')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'a,b,c [nH]')
	plt.savefig(plotdir + 'abc_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], Ms[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'M [nH]')
	plt.savefig(plotdir + 'M_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], Cp1s[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Cp1 [pF]')
	plt.savefig(plotdir + 'Cp1_vs_nfingers.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers[mask], Cp2s[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Cp2 [pF]')
	plt.savefig(plotdir + 'Cp2_vs_nfingers.png')

	x = np.linspace(np.min(L1s), np.max(L1s), 1000)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(L1s[mask], L2s[mask], 'ko')
	ax.plot(x, x, 'r')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'L1 [nH]')
	ax.set_ylabel(r'L2 [nH]')
	plt.savefig(plotdir + 'L1_vs_L2.png')

	x = np.linspace(np.min(Cp1s), np.max(Cp1s), 1000)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(x, x, 'r')
	ax.plot(Cp1s[mask], Cp2s[mask], 'ko')
	#ax.plot(x, np.polyval(cap_l, x), 'r--',
	#	label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Cp1 [pF]')
	ax.set_ylabel(r'Cp2 [pF]')
	plt.savefig(plotdir + 'Cp1_vs_Cp2.png')



	exit()

	L_tot = L_al  + inds
	#alphak = Lk/(L_al + inds)
	fr = 1./(2*np.pi*np.sqrt(L_tot*caps))
	Cc = 0.1883 * pF # Using a number from actual calculations
	Qc =caps/(np.pi*fr*Cc**2*Z0)*(L_tot/L_al)**2
	
	print (Qc)
	exit()
	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(expected_freqs/MHz, caps/pF, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig(plotdir + 'Cp2_vs_design_freq.png')

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
		label='Linear Fit: C (in pF) = {0:1.3e}*N + {1:1.3f}'.format(*cap_l))
	#ax.plot(x, cap_logfit(x, *popt), 'k--', label='log')
	#ax.plot(x, cap_powfit(x, *popt2), 'r--', label='sqrt')
	#ax.plot(x, cap_harmonicfit(x, *popt3), 'k--', label='harmonic')
	ax.grid(which='both')
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Capacitance [pF]')
	plt.savefig(plotdir + 'Cp2_vs_nfingers.png')
	exit()

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
