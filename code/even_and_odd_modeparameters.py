#! /usr/bin/env python3
"""
Script to compare the full 2 port parameters of a network with the even and odd
mode parameters extracted from it.

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

def load_data(fn, nports=2, paramtype='Y'):
	dataset = np.loadtxt(fn, skiprows=9, delimiter=',')
	p = paramtype
	if nports==1:
		dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
		p11 = dataset[:, 1] + 1j*dataset[:, 2]
		tup = list(zip(dataset[:, 0], p11))
		return np.array(tup, dtype=dtypes)
	elif nports==2:
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

def compare_even_and_oddmodes(fn2port, fndiff, fncommon):
	full = load_data(fn2port)
	diff = load_data(fndiff, nports=1)
	common = load_data(fncommon, nports=1)

	#fig, ax = plt.subplots(figsize=(10,10))
	#fracdiff = (full['Y11'] - full['Y22'])/full['Y11']
	#ax.plot(full['frequency'], np.abs(fracdiff), label='Re [(Y11 - Y22)/Y11]')
	#ax.grid()
	#ax.set_xlabel(r'Frequency [MHz]')
	#ax.set_ylabel('|(Y11 - Y22)/Y11|')
	#plt.show()

	#fig, ax = plt.subplots(figsize=(10,10))
	#fracdiff = (full['Y21'] - full['Y12'])/full['Y21']
	#ax.plot(full['frequency'], np.abs(fracdiff))
	#ax.grid()
	#ax.set_xlabel(r'Frequency [MHz]')
	#ax.set_ylabel('|(Y21 - Y12)/Y21|')
	#plt.show()

	nY21 = -full['Y21']
	Yeven = full['Y11'] + full['Y21']
	Yodd = full['Y11'] - full['Y21']

	#fig, ax = plt.subplots(figsize=(10,10))
	#ax.plot(full['frequency'], -full['Y21'].imag, label='Im [-Y21]')
	#ax.plot(full['frequency'], full['Y11'].imag, label='Im [Y11]')
	#ax.grid()
	#ax.legend(loc='upper left')
	#ax.set_xlabel(r'Frequency [MHz]')
	#ax.set_ylabel('Y [1/Ohms]')
	#plt.show()


	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(full['frequency'], nY21.imag, label='Im [-Y21]')
	ax.plot(full['frequency'], Yodd.imag/2, label='Im [Y_odd]/2')
	ax.plot(diff['frequency'], diff['Y11'].imag, label='Im [Y_diff]')
	ax.grid()
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Frequency [MHz]')
	ax.set_ylabel('Y [1/Ohms]')
	plt.savefig('waffleinductor_odd_vs_differentialmode.png')
	plt.show()

	fig, ax = plt.subplots(figsize=(10,10))
	#ax.plot(full['frequency'], nY21.imag, label='Im [-Y21]')
	ax.plot(full['frequency'], 2 * Yeven.imag, label='2 Im [Y_even]')
	ax.plot(common['frequency'], common['Y11'].imag, label='Im [Y_common]')
	ax.grid()
	ax.legend(loc='upper left')
	ax.set_xlabel(r'Frequency [MHz]')
	ax.set_ylabel('Y [1/Ohms]')
	plt.savefig('waffleinductor_even_vs_commonmode.png')
	plt.show()

	fig, ax = plt.subplots(figsize=(10,10))
	L = 1./(diff['Y11'].imag*2*pi*common['frequency']*MHz)/nH
	ax.plot(common['frequency'], L)
	ax.grid()
	ax.set_xlabel(r'Frequency [MHz]')
	ax.set_ylabel('L [nH]')
	plt.savefig('waffleinductor_diffmode_inductance.png')
	plt.show()

	fig, ax = plt.subplots(figsize=(10,10))
	C = common['Y11'].imag/(2*pi*common['frequency']*MHz)/pF
	ax.plot(common['frequency'], C)
	ax.grid()
	ax.set_xlabel(r'Frequency [MHz]')
	ax.set_ylabel('C [pF]')
	plt.savefig('waffleinductor_commonmode_capacitance.png')
	plt.show()

if __name__=="__main__":
    fn1 = datadir + 'waffle_inductor_geometric.csv'
    fn2 = datadir + 'waffle_inductor_differential_geometric.csv'
    fn3 = datadir + 'waffle_inductor_common_geometric.csv'

    compare_even_and_oddmodes(fn1, fn2, fn3)
