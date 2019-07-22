#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pdb

class MetropolisHastingsSampler():

    default_steps = 100

    def __init__(self, lnprob, ndim):
        self.lnprob_fn = lnprob
        self.ndim = ndim
        mat = np.random.randn(ndim, ndim)
        self.cov = 0.5 * (mat + mat.T) # covariance matrix for the gaussian
        self.chain = np.zeros((self.ndim, self.default_steps))

    def get_next_sample(self, x_n, *args):
        candidate = np.random.multivariate_normal(x_n, self.cov)
        a1 = self.lnprob_fn(candidate, *args) - self.lnprob_fn(x_n, *args)
        a2 = 0 # Since I'm using a multivariate gaussian as my proposal distribution

        acceptance_ratio = np.min([1, np.exp(a1 + a2)])
        u = np.random.random_sample()
        x_np1 = candidate if u <= acceptance_ratio else x_n

        return x_np1

    def run_mcmc(self, nsteps, x0, *args):
        self.chain = np.zeros((self.ndim, nsteps))
        self.chain[:, 0] = x0
        for t in range(1, nsteps):
            self.chain[:, t] = self.get_next_sample(self.chain[:, t-1], *args)


def ln_rayleigh(x, sigma):
    if x < 0: return -np.inf
    return -0.5*x**2/sigma**2 + np.log(x) - 2*np.log(sigma)

if __name__=="__main__":
    sigma = 1.0
    x0 = np.abs(np.random.randn())
    sampler = MetropolisHastingsSampler(ln_rayleigh, 1)
    sampler.run_mcmc(1000, x0, sigma)

    x = np.linspace(0, 5, 1000)
    y = (x/sigma**2)*np.exp(-0.5*(x/sigma)**2)

    samples, = sampler.chain
    samples = samples[20:]
    plt.figure(figsize=(10,10))
    plt.hist(samples, histtype='step', density=True)
    plt.plot(x, y, 'r')
    plt.ylabel('')
    plt.xlabel('x')
    plt.show()
