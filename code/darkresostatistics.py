#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi


"""
Now that we have 3 dark reso wafers, I'm hoping to see if we can build some
reasonable statistics on the scatter in the resonance frequencies to see if
there is something systematic going on
"""
mm = 1e-3
um = 1e-6


firstwafer = np.array([290.6, 303.1, 310.5, 312.9, 328.9, 363.2, 367.5, 382.0,
	386.3, 391.1, 400.4, 458.5])
r675 = np.array([299.8, 308.7, 315.7, 316.6, 339.8, 342.8, 397.4, 400.6, 404.3,
	411.17, 411.45, 413.0, 416.3, 421.8, 423.3, 424.9, 442.0, 445.7, 455.5,
	460.4, 477.18, 477.19, 478.9, 485.2, 493.3, 495.4])
cf190401 = np.array([283.7, 294.9, 299.6, 301.78, 301.83, 320.9, 325.7, 331.6,
	352.8, 359.2, 360.1, 364.7, 368.7, 375.8, 377.8, 379.3, 387.7, 389.7, 435.3,
	436.5, 438.3, 440.3, 440.5, 445.0, 453.4, 459.1, 459.2])
design = np.array([306, 318, 321, 324, 327, 330, 333, 336, 339, 351, 354, 369,
	372, 381, 384, 387, 390, 393, 396, 399, 402, 405, 408, 417, 420, 435, 438,
	450, 453, 456, 459, 462, 465, 468, 471, 483])

print (firstwafer.size)
print (r675.size)
print (cf190401.size)
print (design.size)
xticks = 280 + 10*np.arange(23)

tdesign = 1.018*design - 38.9

plt.figure(figsize=(35,10))
plt.plot(design, np.ones(design.size), marker='o', ms=15,
		linestyle='None', label='Design')
plt.plot(firstwafer, np.ones(firstwafer.size) + 1, marker='d', ms=15,
		linestyle='None', label='First Wafer')
plt.plot(r675, np.ones(r675.size) + 0.5, marker='h', ms=15, linestyle='None',
		label='R675')
plt.plot(cf190401, np.ones(cf190401.size) - 0.5, marker='v', ms=15,
		linestyle='None', label='CF190401')
plt.grid()
plt.xticks(xticks)
plt.xlim((270, 500))
plt.yticks([0, 1, 2])
#plt.axis('tight')
plt.legend(loc='upper left')
plt.xlabel('Frequency [MHz]')
plt.savefig('frequencyscatter_darkreso.png')

#plt.figure(figsize=(35,10))
#plt.plot(design, np.ones(design.size), marker='o', ms=15,
#		linestyle='None', label='Design')
#plt.plot(cf190401, np.ones(cf190401.size) + 0.5, marker='v', ms=15,
#		linestyle='None', label='CF190401')
#plt.grid()
#plt.xticks(xticks)
#plt.xlim((270, 500))
#plt.yticks([0, 1, 2])
##plt.axis('tight')
#plt.legend(loc='upper left')
#plt.xlabel('Frequency [MHz]')
#plt.savefig('frequencyscatter_darkreso_cf190401.png')




#plt.figure(figsize=(35,10))
#plt.plot(design, np.ones(design.size), marker='o', ms=15,
#		linestyle='None', label='Design')
#plt.plot(r675, np.ones(r675.size) - 0.5, marker='h', ms=15, linestyle='None',
#		label='R675')
#plt.grid()
#plt.xticks(xticks)
#plt.xlim((270, 500))
#plt.yticks([0, 1, 2])
##plt.axis('tight')
#plt.legend(loc='upper left')
#plt.xlabel('Frequency [MHz]')
#plt.savefig('frequencyscatter_darkreso_r675.png')


# Want to generate the cumulative distribution functions for the different
# resonators
X = np.r_[270:510:100j]

def generate_cdf(x, y):
	cdf = np.zeros_like(x)
	N = y.size
	for i in range(x.size):
		cdf[i] =  y[y <= x[i]].size
	cdf /= N
	return cdf

def get_KS_statistic(cdf1, cdf2):
	return np.max(np.abs(cdf1 - cdf2))


design_cdf = generate_cdf(X, design)
firstwafer_cdf = generate_cdf(X, firstwafer)
r675_cdf = generate_cdf(X, r675)
cf190401_cdf = generate_cdf(X, cf190401)

firstwafer_KS = get_KS_statistic(design_cdf, firstwafer_cdf)
r675_KS = get_KS_statistic(design_cdf, r675_cdf)
cf190401_KS = get_KS_statistic(design_cdf, cf190401_cdf)

alpha = 0.05
def reject_nullhypothesis(alpha, n, m, Dnm):
	c = np.sqrt(-0.5*np.log(alpha))
	print (c)
	thresh = c * np.sqrt((n+m)/n/m)
	return Dnm > thresh

print (reject_nullhypothesis(alpha, design.size, firstwafer.size, firstwafer_KS))
print (reject_nullhypothesis(alpha, design.size, r675.size, r675_KS))
print (reject_nullhypothesis(alpha, design.size, cf190401.size, cf190401_KS))


plt.figure(figsize=(10,10))
plt.plot(X, design_cdf, label='Design')
plt.plot(X, firstwafer_cdf, label='First Wafer, KS=%1.2f'%(firstwafer_KS))
plt.plot(X, r675_cdf, label='R675, KS=%1.2f'%(r675_KS))
plt.plot(X, cf190401_cdf, label='CF190401, KS=%1.2f'%(cf190401_KS))
plt.legend(loc='upper left')
plt.grid()
plt.xlabel('Resonance Frequency [MHz]')
plt.ylabel('Cumulative Probability')
plt.savefig('darkreso_kolmogorov_comparison.png')
exit()


design = np.array([306, 321, 318, 471, 468, 483, 339, 324, 327, 462, 465, 450,
	336, 333, 330, 459, 456, 453, 351, 390, 393, 396, 399, 438, 354, 387, 384,
	405, 402, 435, 369, 372, 381, 408, 417, 420])


Nrow, Ncol = 6, 6
designarray = design.reshape(Nrow, -1)
designarray = designarray
spacing = 11*mm
print (designarray)

row, col = np.meshgrid(np.arange(Nrow), np.arange(Ncol), indexing='xy')
X0, Y0 = 1.52*mm, 0*mm
X0, Y0 = -spacing/2, -spacing/2
X = X0 - (Ncol - 2*row - 1)/2 * spacing
Y = Y0 - (Nrow - 2*col - 1)/2 * spacing

R = np.sqrt((X + X0)**2 + (Y+Y0)**2)

print (X)
print (Y)

# Now lets model the non-uniformity across the wafer
delta0 = 0e-3
delta1 = -2e-5
#delta = delta0 + delta1*(R/mm)**2
delta = delta0 + delta1*(R/mm)**2

plt.figure(figsize=(10,10))
plt.plot((R/mm).flatten(), delta.flatten()*1e3, 'bo')
plt.xlabel('R [mm]')
plt.ylabel('delta [x10^-3]')
plt.grid()
plt.savefig('delta_vs_R.png')

sigdelta = 2.0e-3



# Lets collect some statistics for the shift in resonance frequency

Ntries = 10000
samplesets = np.zeros((Ntries, design.size))
for i in range(Ntries):
	deltasampled = delta + np.random.randn(Nrow, Ncol)*sigdelta
	fmeasarray = (designarray-0)*(1 + deltasampled)
	samplesets[i] = fmeasarray.flatten()

# First just print the mean and std
sigmas = np.std(samplesets, axis=0)
means = np.mean(samplesets, axis=0)
print (design)
print (means)
print (sigmas)

plt.figure(figsize=(10,10))
plt.plot(design, means, 'bo')
plt.grid()
plt.xlabel('Design Frequency [MHz]')
plt.ylabel('Simulated Actual Frequency [MHz]')
plt.savefig('simulated_freq_shift.png')
plt.close()
cf190401matrix = np.zeros((Nrow*Ncol, cf190401.size))
r675matrix = np.zeros((Nrow*Ncol, r675.size))

# Make sample histograms
for i in range(Nrow*Ncol):
	plt.figure(figsize=(10,10))
	plt.hist(samplesets[:, i], histtype='step', density='True')
	x = np.r_[-5*sigmas[i]:5*sigmas[i]:100j] + means[i]
	pdf = 1./(2*pi)**0.5/sigmas[i]*np.exp(-0.5*(x - means[i])**2/sigmas[i]**2)
	cf190401pdfs = 1./(2*pi)**0.5/sigmas[i]*np.exp(-0.5*(cf190401 - means[i])**2/sigmas[i]**2)
	cf190401pdfs[cf190401pdfs < 0.01] = np.nan
	cf190401matrix[i] = cf190401pdfs
	r675pdfs = 1./(2*pi)**0.5/sigmas[i]*np.exp(-0.5*(r675 - means[i])**2/sigmas[i]**2)
	r675pdfs[r675pdfs < 0.01] = np.nan
	r675matrix[i] = r675pdfs
	plt.plot(x, pdf, 'r')
	plt.plot(cf190401, cf190401pdfs, 'kv', ms=30)
	plt.plot(r675, r675pdfs, 'g^', ms=30)
	plt.grid()
	#plt.ylim((0,1))
	#plt.yscale('log')
	plt.savefig('%dMHz_resohist.png'%design[i])
	plt.close()


# Can now also evaluate the probability
for i in range(cf190401.size):
	#if i > 0: continue
	plt.figure(figsize=(10,10))
	plt.plot(design, cf190401matrix[:, i], 'bo')
	plt.xlabel('Design Frequency [MHz]')
	plt.ylabel('Probabilty ')
	plt.grid()
	plt.xlim((280, 500))
	plt.ylim(bottom=0)
	plt.savefig('cf190401_%1.1fMHzreso_probability.png'%(cf190401[i]))
	plt.close()


