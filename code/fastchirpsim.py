#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy import fftpack, optimize

MHz = 1e6
ms = 1e-3

# The transfer function of a single resonator
def S21(f, fr, dQr, dQe):
	x = (f - fr)/fr
	return 1 - dQe/(dQr + 1j * 2 * x)

# The time domain waveform describing the chirp waveform
def zin(t, z0, f0, k, phi0):
	exponent = phi0 + 2*pi*f0*t + pi*k*t**2
	return z0*np.exp(1j*exponent)

def generate_chirp(N, sample_rate, lo):
	f0 =  sample_rate / 2.0
	f1 = -sample_rate / 2.0
	T = N * sample_rate
	k = (f1 - f0) / T

	# Use the full scale ADC output
	z0 = 1<<14
	t = np.arange(N)
	z_chirp = zin(t, z0, lo + f0, k, 0)

	return t, z_chirp

def lorentzian(x, gamma):
	return (gamma/pi)/(x**2 + gamma**2)

def psd_lorentzian(f, f0, gamma, A):
	return A * lorentzian(f - f0, gamma)

def psd_model(f, f0, gamma1, gamma2, A1, A2):
	return A1 * lorentzian(f-f0, gamma1) + A2*lorentzian(f-f0, gamma2)

#def chisq(theta, x, y, yerr):
#	inv_sigma2 = 1/yerr**2
#	return 0.5*np.sum((y - psd_model(x, *theta))**2*inv_sigma2)

def chisq(theta, x, y):
	return 0.5 * np.sum((y - psd_model(x, *theta))**2)


if __name__=="__main__":
	navg = 100
	ntotal = 8192
	nchirp = 1700
	sample_rate = 100e6
	fast_chirp_rate = sample_rate / ntotal
	chirp_rate = sample_rate / navg
	chirp_period = 1e3 / chirp_rate # in ms

	rffreq = 450e6
	t_chirp, z_chirp = generate_chirp(nchirp, sample_rate, rffreq)
	t = np.arange(ntotal)
	z_in = np.zeros_like(t, dtype=np.complex64)
	z_in[:nchirp] = z_chirp
	chirp_mask = (z_in > 0)
	listen_mask = (z_in == 0)

	z_in_fft = fftpack.fftshift(fftpack.fft(z_in))
	freqs = fftpack.fftshift(fftpack.fftfreq(ntotal, 1./sample_rate)) + rffreq
	#plt.plot(t, np.abs(z_in_fft))
	#plt.show()


	Qr = 10000
	Qe = 20000
	Qe *= (1 + 0.2j)
	fr = 440e6
	# Transfer function
	H = S21(freqs, fr, 1./Qr, 1./Qe)
	#plt.plot(t, np.abs(H))
	#plt.grid()
	#plt.show()

	z_out_fft = z_in_fft * H

	z_out = fftpack.ifft(z_out_fft)

	#plt.plot(t[listen_mask], np.real(z_out[listen_mask]))
	#plt.plot(t, np.real(z_out))
	#plt.grid()
	#plt.show()


	x = t[listen_mask]
	y = z_out[listen_mask]

	# Now we have a resonant response from the fast chirp
	resonance = fftpack.fft(y)
	fs = fftpack.fftshift(fftpack.fftfreq(x.size, 1./sample_rate)) + rffreq
	zfc = np.abs(resonance)
	peak = np.max(zfc)
	zfc /= peak
	fs /= MHz
	fw = 1.5
	f0 = fr/MHz
	ok = (f0 - fw < fs) & (fs < f0 + fw)
	#plt.plot(fs/MHz, zfc)
	#plt.grid()
	#plt.show()

	# Now test and see the relation between this resonance and the main one

	gamma1 = f0/(Qr)
	A1 = np.max(zfc)*pi/gamma1
	A2 = A1/4
	gamma2 = 0.5*gamma1
	p0 = [f0, gamma1, gamma2, A1, A2]
	popt, pcov = optimize.curve_fit(psd_model, fs[ok], zfc[ok], p0=p0,
			method = 'lm')
	f0, gamma1, gamma2, A1, A2 = popt
	#result = optimize.minimize(chisq, p0, args=(fs[ok], zfc[ok]))
	#f0, gamma, A = result['x']
	Qr1 = f0/(2*gamma1)
	Qr2 = f0/(2*gamma2)
	print (f0, Qr1, Qr2, A1, A2)
	#chisq_val = chisq(result["x"], fs[ok], zfc[ok])
	fs_fit = np.linspace(fs[ok][0], fs[ok][-1], 1000)

	plt.plot(fs[ok], zfc[ok], 'b', label='Data')
	plt.plot(fs[ok], psd_model(fs[ok], *popt), 'r', label='Fit')
	plt.legend()
	plt.grid()
	plt.show()

















