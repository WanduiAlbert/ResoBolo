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
		for i in range(11): fp.readline()
		ypos = float(fp.readline().split(",")[0])
		dist2capedge = float(fp.readline().split(" ")[-1])
		distfromedge = gpwidth - dist2capedge
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

	C_est = np.mean(np.gradient(Y11.imag))/np.mean(2*pi*f)
	C_est = 0.024#*pF
	L_est = 0.206#*nH
	R_est = 1e-5
	#print (C_est/pF, L_est/nH, R_est)
	#print (C_est, L_est, R_est)


	p0 = [C_est, L_est, R_est]
	popt, pcov = curve_fit(admittance_model, f, np.abs(Y11), p0=p0, method='lm')
	C_fit, L_fit, R_fit = popt
	y_fit = admittance_model(f, *popt)
	y_est = admittance_model(f, *p0)

	fig, ax =plt.subplots(figsize=(10,10))
	#ax.plot(f/MHz, y_est, 'go', label='Guess')
	ax.plot(f/MHz, np.abs(Y11), 'b', label='Simulation')
	ax.plot(f/MHz, y_fit, 'r-',
				label="Fit: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_fit,
					L_fit, R_fit))
	#ax.plot(f/MHz, np.abs(y_fit), 'k--',
	#			label="Fit Zc = %1.3fOhms v0 = %1.3fMHz eps=%1.3e "%(Zc_fit,\
	#					0.5*w0_fit/MHz, eps_fit))
	#ax.plot(f/MHz, np.abs(y_est), 'r-',
	#			label="Fit Zc = %1.3fOhms v0 = %1.3fMHz eps=%1.3e "%(Zc_est,\
	#					0.5*w0_est/MHz, eps_est))
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
	plt.close('all')

	plt.close('all')
	return C_fit*pF, L_fit*nH, R_fit
	#return Zc_fit, w0_fit, eps_fit

if __name__=="__main__":
	L_al = 10 * nH
	filename = datadir + "simple_cap_withGP_447um_verticalslice_level1.csv"
	#filename = datadir + "simple_cap_withGP_955um_verticalslice_level2.csv"
	#filename = datadir + "simple_cap_withGP_875um_verticalslice_level1.csv"
	#dataloader =  load_data(filename, nports=1, paramtype='Y')
	ypos, xre, xim, yre, yim  = np.loadtxt(filename, delimiter=',', unpack=True)
	Jx = xre + 1j*xim
	Jy = yre + 1j*yim

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(ypos, np.abs(Jy), 'b', label='|Jy|')
	ax.plot(ypos, np.abs(Jx), 'r', label='|Jx|')
	ax.grid(which='both')
	ax.set_xlabel(r'Y position [um]')
	ax.set_ylabel(r'Current Density [A/m]')
	ax.legend(loc='upper right')
	plt.show()


	Jxfft = np.fft.fft(Jx)
	Jyfft = np.fft.fft(Jy)
	freq = np.fft.fftshift(np.fft.fftfreq(ypos.size, np.abs(ypos[1] - ypos[0])))
	#freq = np.fft.fftfreq(ypos.size, np.abs(ypos[1] - ypos[0]))

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(freq, np.abs(Jyfft), 'b', label='|Jy|')
	ax.plot(freq, np.abs(Jxfft), 'r', label='|Jx|')
	ax.grid(which='both')
	ax.set_xlabel(r'Y position [um]')
	ax.set_ylabel(r'Current Density [A/m]')
	ax.legend(loc='upper right')
	plt.show()
