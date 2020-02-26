#! /usr/bin/env python3
"""

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
plotdir = 'parallel_platecap/'

def load_data(fn, nports=1, paramtype='Y'):
	p = paramtype
	if nports==1:
		dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
	elif nports==2:
		dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
			(p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
	loopcondition = True
	with open(fn) as fp:
		for i in range(8): fp.readline()
		gpwidth = float(fp.readline().split(" ")[-1])
		dist2capedge = float(fp.readline().split(" ")[-1])
		distfromedge = gpwidth - dist2capedge
		#print (distfromedge)
		dataset = []
		freqs = []
		while loopcondition:
			line = fp.readline()
			if line is None: break
			#print (line)
			#print (line.startswith(" "))
			if line.startswith(" "):
				tup = list(zip(np.array(freqs), np.array(dataset)))
				yield distfromedge, np.array(tup, dtype=dtypes)
				dataset = []
				for i in range(7): fp.readline()
				gpwidth = float(fp.readline().split(" ")[-1])
				dist2capedge = float(fp.readline().split(" ")[-1])
				distfromedge = gpwidth - dist2capedge
				continue
			elif line == "": break
			if nports==1:
				freq, re, im = list(map(lambda x: float(x), line.split(",")))
				freqs.append([freq])
				dataset.append([re + 1j*im])
			elif nports==2:
				freq,\
				re11,im11,\
				re12, im12,\
				re21, im21,\
				re22, im22 = list(map(
					lambda x: float(x), fp.readline().split(",")))
				freqs.append([freq])
				dataset.append([
						re11 + 1j*im11,\
						re12 + 1j*im12,\
						re21 + 1j*im21,\
						re22 + 1j*im22])
	return

def admittance_model(x, C, L, R):
	w = 2*pi*x
	return np.abs(1./(R + 1j*w*L*nH + 1./(1j*w*C*pF)))
	#return np.abs(1./R + 1./(1j*w*L) + 1j*w*C)
	#return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
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
	return np.sum((y - admittance_model(x, *theta))**2)
	#return np.sum((y - tline_model(x, *theta))**2)

def get_S21(Y):
	Y0 = 1/Z0
	DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
	return -2 * Y['Y21']*Y0/DY

def real_of_complex(z):
	return np.hstack([z.real, z.imag])

def complex_of_real(y):
	r, i = np.split(y, 2)
	return r + 1j*i

def get_cap_params(dist, Ydata):
	savename = "dist{0:3.1f}um_fromedge".format(dist)
	f = Ydata['frequency'] * MHz
	Y11 = Ydata['Y11']

	#C_est = np.mean(np.gradient(Y11.imag))/np.mean(2*pi*f)/pF
	#L_est = 0.0#*nH
	C_est = 3.0#*pF
	L_est = 0.3#*nH
	R_est = 1e-15
	#print (C_est/pF, L_est/nH, R_est)
	#print (C_est, L_est, R_est)


	p0 = [C_est, L_est, R_est]
	popt, pcov = curve_fit(admittance_model, f, np.abs(Y11), p0=p0, method='lm')
	C_fit, L_fit, R_fit = popt
	y_fit = admittance_model(f, *popt)
	y_est = admittance_model(f, *p0)

	fig, ax =plt.subplots(figsize=(10,10))
	#ax.plot(f/MHz, y_est, 'go', label='Guess')
	ax.plot(f/MHz, np.imag(Y11), 'b', label='Simulation')
	ax.plot(f/MHz, y_fit, 'r-',
				label="Fit: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_fit,
					L_fit, R_fit))
	#ax.plot(f/MHz, np.abs(y_fit), 'k--',
	#			label="Fit Zc = %1.3fOhms v0 = %1.3fMHz eps=%1.3e "%(Zc_fit,\
	#					0.5*w0_fit/MHz, eps_fit))
	#ax.plot(f/MHz, np.abs(y_est), 'g--',
	#			label="Guess: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_est,
	#				L_est, R_est))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|Y11| [MOhs]')
	ax.grid()
	ax.axis('tight')
	#plt.show()
	#exit()
	plt.savefig(plotdir + savename + "Y11.png")

	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(f/MHz, np.abs(Y11) - np.abs(y_fit),
				label="Fit: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_fit,
					L_fit, R_fit))
	#ax.scatter(f/MHz, np.abs(Y11) - np.abs(y_fit),
	#			label="Fit Zc = %1.3fOhms v0 = %1.3fMHz eps=%1.3e "%(Zc_fit,\
	#					0.5*w0_fit/MHz, eps_fit))
	ax.legend()
	ax.set_xlabel('Frequency [MHz]')
	ax.set_ylabel('|Y11| [MOhs]')
	ax.grid()
	ax.axis('tight')
	plt.savefig(plotdir + savename + "Y11_residuals.png")
	#plt.show()
	plt.close('all')

	plt.close('all')
	return C_fit*pF, L_fit*nH, R_fit
	#return Zc_fit, w0_fit, eps_fit

if __name__=="__main__":
	L_al = 10 * nH
	filename = datadir + "cap2gnd_10umfromedge_new.csv"
	dataloader =  load_data(filename, nports=1, paramtype='Y')
	edgedist = []
	Cs = []
	Ls = []
	Rs =  []
	for dist, ydata in dataloader:
		print (dist)
		C, L, R = get_cap_params(dist, ydata)
		edgedist.append(dist)
		Cs.append(C)
		Ls.append(L)
		Rs.append(R)

	edgedist = np.array(edgedist)
	Cs = np.array(Cs)
	Ls = np.array(Ls)
	Rs = np.array(Rs)
	print (Cs/pF)
	print (Ls/1e-9)
	#exit()

	#pC = np.polyfit(np.log(disttoGNDplane)[mask], Cs[mask]/pF, 1)
	#pL = np.polyfit(np.log(disttoGNDplane)[mask], Ls[mask]/nH, 1)
	#labelC = "m = %1.5f pF/um, b = %1.5f pF"%(pC[0], pC[1])
	#labelL = "m = %1.5f nH/um, b = %1.5f nH"%(pL[0], pL[1])

	#h = np.linspace(10, 5000, 100)
	#Cfit = np.polyval(pC, np.log(h))
	#Lfit = np.polyval(pL, np.log(h))

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(edgedist, Cs/pF, 'bo', ms=12)
	#ax.plot(h, Cfit, 'k', label=labelC)
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Distance to GND plane [um]')
	ax.set_ylabel(r'Capacitance [pF]')
	#ax.set_ylabel(r'Characteristic Impedance [Ohms]')
	plt.savefig(plotdir + 'C_vs_design_size.png')

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(edgedist, Ls/nH, 'bo', ms=12)
	#ax.plot(h, Lfit, 'k', label=labelL)
	ax.grid(which='both')
	#ax.legend(loc='upper left')
	ax.set_xlabel(r'Distance to GND plane [um]')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	#ax.set_ylabel(r'Characteristic Impedance [Ohms]')
	plt.savefig(plotdir + 'L_vs_design_size.png')
	plt.show()
	exit()
	def cap_logfit(N, A, B):
		return A*np.log(N) + B
	def cap_powfit(N, A, B):
		return A*N**0.5 + B
	def cap_harmonicfit(N, A, B):
		return A*(1. - 1./N) + B

	Lk = 0.4 * L_al

	alphak = Lk/(L_al + w0s)


	cap_p = np.polyfit(1./Nfingers, Zcs/pF, 1)
	cap_l = np.polyfit(Nfingers, Zcs/pF, 1)
	p0 = [cap_p[0], cap_p[1]]
	#popt, _ = curve_fit(cap_logfit, Nfingers, Zcs/pF, p0=p0, method='lm')
	#popt2, _ = curve_fit(cap_powfit, Nfingers, Zcs/pF, p0=p0, method='lm')
	#popt3, _ = curve_fit(cap_harmonicfit, Nfingers, Zcs/pF, p0=p0, method='lm')
	x = np.linspace(250, 550, 100)
	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers, Zcs/pF, 'ko')
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
	ax.plot(disttoGNDplane/MHz, w0s/nH, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	plt.savefig(plotdir + 'ind_vs_design_freq.png')

	ind_p = np.polyfit(Nfingers, w0s/nH, 2)

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(Nfingers, w0s/nH, 'ko')
	ax.plot(x, np.polyval(ind_p, x), 'k--',
		label='Quadratic Fit:a={0:1.3e} b={1:1.3f} c={2:1.3f}'.format(*ind_p))
	ax.grid(which='both')
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Number of Fingers')
	ax.set_ylabel(r'Parasitic Inductance [nH]')
	plt.savefig(plotdir + 'ind_vs_nfingers.png')


	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(disttoGNDplane/MHz, alphak, 'ko')
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
	ax.plot(disttoGNDplane/MHz, epsilons*1e5, 'ko')
	ax.grid(which='both')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.set_ylabel(r'Resistance (x10^-5) [Ohms]')
	plt.savefig(plotdir + 'epsilons_vs_design_freq.png')

	plt.close('all')
	L_tot = w0s + L_al
	freqs = 1./np.sqrt(L_tot * Zcs)/2/pi
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

	p = np.polyfit(disttoGNDplane/MHz, freqs/MHz, 1)
	print (p)
	m, b = p
	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(disttoGNDplane/MHz, freqs/MHz)
	ax.plot(disttoGNDplane/MHz, np.polyval(p, disttoGNDplane/MHz), 'k',
			label='m=%1.3f\nb=%1.3f'%(m,b))
	ax.grid(which='both')
	ax.set_ylabel(r'True Frequency [MHz]')
	ax.set_xlabel(r'Design Frequency [MHz]')
	ax.legend(loc='upper left')
	plt.savefig(plotdir + 'design_vs_actual_freq.png')

