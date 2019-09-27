#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import scipy.optimize as opt
import emcee
import corner

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
#design = np.array([306, 318, 321, 324, 327, 330, 333, 336, 339, 351, 354, 369,
#	372, 381, 384, 387, 390, 393, 396, 399, 402, 405, 408, 417, 420, 435, 438,
#	450, 453, 456, 459, 462, 465, 468, 471, 483])


design = np.array([306, 321, 318, 471, 468, 483, 339, 324, 327, 462, 465, 450,
	336, 333, 330, 459, 456, 453, 351, 390, 393, 396, 399, 438, 354, 387, 384,
	405, 402, 435, 369, 372, 381, 408, 417, 420])


Nrow, Ncol = 6, 6
designarray = design.reshape(Nrow, -1)
designarray = designarray
spacing = 11*mm

row, col = np.meshgrid(np.arange(Nrow), np.arange(Ncol), indexing='xy')
X0, Y0 = 1.52*mm, 0*mm
X0, Y0 = -spacing/2, -spacing/2
X = X0 - (Ncol - 2*row - 1)/2 * spacing
Y = Y0 - (Nrow - 2*col - 1)/2 * spacing

R = np.sqrt((X + X0)**2 + (Y+Y0)**2)

# Now lets model the non-uniformity across the wafer
delta0 = 1e-3
r0 = 0
delta1 = -2e-5
sigma = 1e-3

def lnlike(theta, x, freq_meas, freq_err):
    delta0, r0, delta1, sigma = theta
    freq_design, R = np.split(x, 2)
    delta = delta0*(R - r0)**2 + delta1
    model = freq_design*(1 + delta)
    inv_sigma2 = 1./(sigma**2 + freq_err**2)
    return -0.5*np.sum((freq_meas - model)**2*inv_sigma2)

def lnprior(theta):
    delta0, r0, delta1, sigma = theta
    if -10 < delta0 < 10 and 0 < r0 < 10 and -10 < delta1 < 10 and 0 < sigma < 20:
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

nll = lambda *args: -lnlike(*args)

y = cf190401
yerr = np.ones_like(y.size)
indices = np.arange(Nrow*Ncol)
sampleindices = np.random.choice(indices, y.size, replace=False)

x = np.hstack([designarray.flatten()[sampleindices], R.flatten()[sampleindices]])

p0 = [delta0, r0, delta1,sigma]
result = opt.minimize(nll, p0, args=(x, y, yerr))

print (result["x"])

ndim, nwalkers = 4, 100
pos = [result["x"] + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 800)
chain = sampler.chain
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

print (sampler.chain.shape)

fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True)
ax[0].plot(chain[::2, :, 0], 'k')
ax[1].plot(chain[::2, :, 1], 'k')
ax[2].plot(chain[::2, :, 2], 'k')
ax[3].plot(chain[::2, :, 3], 'k')

ax[3].set_xlabel('Step')
ax[0].set_ylabel('delta0')
ax[1].set_ylabel('r0')
ax[2].set_ylabel('delta1')
ax[3].set_ylabel('sigma')

plt.savefig('darkreso_walkers.png')
exit()

fig = corner.corner(samples, labels=["delta0", "r0", "delta1", "sigma"],
                      truths=result["x"])
fig.savefig("darkreso_triangle.png")

