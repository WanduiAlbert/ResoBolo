#! /usr/bin/env python3
"""
Code for rational function fitting to Y, Z and S parameters
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
plotdir = 'cap_fig/'

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

def hermitian_conjugate(x):
	return np.conjugate(x.T)

def rational(x, p, q):
	return np.polyval(p,x)/np.polyval(q + [1.0],x)

def fit_rationalfn(s, H):
	max_order = 10
	lambdas = []
	orders = {}
	solutions = {}
	numdenom = {}
	for i in range(max_order):
		for j in range(max_order):
			if i + j == 0: continue
			if i+j > max_order: continue
			P, Q = i, j

			B = np.vander(s, P,increasing=True)
			C = -H[:, np.newaxis] * np.vander(s, Q,increasing=True)

			A1 = np.hstack([np.real(B), np.real(C)])
			A2 = np.hstack([np.imag(B), np.imag(C)])
			A = np.vstack([A1, A2])
			D = np.dot(hermitian_conjugate(A), A)
			w, v = np.linalg.eigh(D)
			n = P+Q+2
			l = np.min(np.abs(w))
			q = v[:, 0]
			if n in orders:
				l_old = orders[n]
				if l < l_old:
					orders[n] = l
					solutions[n] = q
					numdenom[n] = (P, Q)
			else:
				orders[n] = l
				solutions[n] = q
				numdenom[n] = (P, Q)

	x = list(orders.keys())
	y = list(orders.values())
	sel_order = x[np.argmin(y)]
	print (sel_order-2, numdenom[sel_order], solutions[sel_order])
	plt.scatter(x, y)
	plt.xlabel('P + Q + 2')
	plt.ylabel('|Lambda min|')
	plt.ylim(bottom=0)
	plt.show()

if __name__=="__main__":
	filenames = glob.glob(datadir + "*_with_boundary.csv")
	filenames.sort(key=lambda x:int(x.split('/')[-1].split('.')[0].split('_')[1][:3]))
	fn = filenames[0]
	savename = os.path.split(fn)[-1].split(".")[0]
	Ydata = load_data(fn)
	f = Ydata['frequency'] #* MHz
	s = 1j*2*pi*f
	fit_rationalfn(s, Ydata['Y11'])







