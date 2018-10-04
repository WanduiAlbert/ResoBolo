#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.special import wofz

def Lorentzian(x, x0, gamma):
    return gamma/pi * 1/((x-x0)**2 + gamma**2)

def shifted_sum(x, x0, x1, gamma):
    return 0.5 * Lorentzian(x, x0, gamma) + 0.5 * Lorentzian(x, x1, gamma)

def scaled_sum(x, x0, gamma1, A1, gamma2, A2):
    return A1 * Lorentzian(x,x0, gamma1) + A2 * Lorentzian(x, x0, gamma2)

def voigt(x, x0, gamma, sigma, A):
    z = ((x - x0) + 1j*gamma)/sigma/2**0.5
    return A * np.real(wofz(z))/sigma/(2*pi)**0.5

def shifted_chisq(theta, x, y,yerr):
    return 0.5*np.sum((y - shifted_sum(x,*theta))**2)

def scaled_chisq(theta,x,y,yerr):
    return 0.5*np.sum((y - scaled_sum(x,*theta))**2/yerr**2)

def voigt_chisq(theta,x,y,yerr):
    return 0.5*np.sum((y - voigt(x,*theta))**2/yerr**2)

def single_chisq(theta,x,y,yerr):
    x0, gamma, A = theta
    return 0.5*np.sum((y - A*Lorentzian(x,x0,gamma))**2/yerr**2)

f0 = 450
Nchirps = 1500
N = 5000
BW = 2e-3
freq_plot = np.linspace(f0-BW/2, f0+BW/2, N)
freq_range = freq_plot[np.newaxis, :]
f = np.zeros((Nchirps, N)) + freq_range
sigma = 0.01e-3
gamma = 0.1e-3
gaussian = np.random.randn(Nchirps)*sigma
#gaussian = np.zeros(Nchirps)

fr = f0 + gaussian
fr = fr[:, np.newaxis]
x = f - fr
zf = Lorentzian(f, fr, gamma)
zf_sigma = 10.
zf_noise = np.random.randn(Nchirps, N) * zf_sigma
zf += zf_noise
zfc_mean = np.average(zf, axis=0)
zfc_err = np.std(zf_noise, axis=0)
zfc_err = np.ones(N)*zf_sigma

#plt.imshow(zf_noise)
#plt.colorbar()
#plt.show()
#exit()
A1 = 0.7
A2 = A1 * 0.4

#fit = "single"
#fit = "scaled"
fit = "voigt"

if fit == "single":
    p0 = [f0, gamma*2.5, A1]
    result = opt.minimize(single_chisq, p0, args=(freq_plot, zfc_mean, zfc_err), method="Nelder-Mead")
    fr_fit, gamma_fit, A_fit = result["x"]
    z_fit =  A_fit * Lorentzian(freq_plot, fr_fit, gamma_fit)
elif fit == "scaled":
    p0 = [f0, gamma*2.5, A1, gamma/2., A2 ]
    result = opt.minimize(scaled_chisq, p0, args=(freq_plot, zfc_mean, zfc_err), method="Nelder-Mead")
    z_fit =  scaled_sum(freq_plot, *result["x"])
elif fit == "voigt":
    p0 = [f0, gamma*2.5, sigma, A1 ]
    result = opt.minimize(voigt_chisq, p0, args=(freq_plot, zfc_mean, zfc_err), method="Nelder-Mead")
    z_fit =  voigt(freq_plot, *result["x"])

residuals = zfc_mean - z_fit

print (p0)
print (result["x"])

#fig, ax = plt.subplots(figsize=(10,10))
#ax.errorbar(freq_plot, zfc_mean, yerr=zfc_err, fmt='b.', label="Data")
#ax.plot(freq_plot, z_fit, 'r', label="Fit")
#ax.grid(which='both')
#ax.legend(loc='best')
#ax.set_xlabel('Frequency [MHz]')
#ax.set_ylabel('PSD')


fig, ax = plt.subplots(figsize=(10,10))
ax.plot(freq_plot, residuals, 'b')
ax.grid(which='both')
ax.set_xlabel('Frequency [MHz]')
ax.set_ylabel('PSD Residuals')



plt.show()
