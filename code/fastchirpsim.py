#! /usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi
from scipy import fftpack, optimize, signal
from scipy.special import erf, fresnel

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

def psd_root_lorentzian(f, f0, gamma, A):
	return (A/pi)/np.sqrt((f-f0)**2 + gamma**2)

def psd_model(f, f0, gamma1, gamma2, A1, A2):
	return A1 * lorentzian(f-f0, gamma1) + A2*lorentzian(f-f0, gamma2)

#def chisq(theta, x, y, yerr):
#	inv_sigma2 = 1/yerr**2
#	return 0.5*np.sum((y - psd_model(x, *theta))**2*inv_sigma2)

def chisq(theta, x, y):
	return 0.5 * np.sum((y - psd_model(x, *theta))**2)

def ringdown(t,  tau, B, t0):
	return B * np.exp(-(t - t0)/tau)

def get_chirpspectrum(nchirp, ntotal, sample_rate, rffreq, *reso_args):
	fr, dQr, dQe = reso_args
	z0 = 1<<14
	t_chirp, z_chirp = generate_chirp(nchirp, sample_rate, rffreq)
	t = np.arange(ntotal)/sample_rate
	z_in = np.zeros_like(t, dtype=np.complex64)
	z_in[0] = 1
	#z_in[:nchirp] = z_chirp
	listen_mask = np.ones(ntotal, dtype=bool)
	listen_mask[:nchirp] = 0
	listen_mask[-20:] = 0
	z_in_fft = fftpack.fftshift(fftpack.fft(z_in))
	freqs = fftpack.fftshift(fftpack.fftfreq(ntotal, 1./sample_rate)) + rffreq

	# Transfer function
	H = S21(freqs, *reso_args)
	z_out_fft = z_in_fft * H
	z_out = fftpack.ifft(fftpack.ifftshift(z_out_fft))
	#plt.plot(np.real(z_out)[2:20], 'k', label='Data')
	#plt.plot(np.imag(z_out)[2:20], 'b', label='Data')
	#plt.grid()
	#plt.show()
	#exit()
	x = t
	y = z_out * listen_mask
	popt, pcov = optimize.curve_fit(ringdown, x[listen_mask], y[listen_mask].real,
			p0=[dQr, np.max(np.abs(y[listen_mask].real)), x[listen_mask][0]],
			method='lm')
	yfit = ringdown(x[listen_mask], *popt)
	tau, A, B0 = popt
	Qr_fit = tau * fr * pi
	print (tau, A, B0)
	#listen_fft = fftpack.fftshift(fftpack.fft(listen_mask))
	#plt.figure()
	#plt.semilogy(freqs, np.abs(listen_fft), 'k', label='listen mask Real')
	#plt.grid()
	#plt.show()
	#exit()
	plt.figure()
	plt.plot(x, y.real, 'k', label='Timestream Real')
	#plt.plot(x[listen_mask], y[listen_mask].real, 'k', label='Timestream Real')
	plt.plot(x[listen_mask], yfit, 'r--',
			label='Ringdown Fit: Qr=%1.0f'%(Qr_fit))
	plt.legend()
	plt.grid()
	plt.show()
	#exit()
	# Now we have a resonant response from the fast chirp
	resonance = fftpack.fftshift(fftpack.fft(y))
	fs = fftpack.fftshift(fftpack.fftfreq(x.size, 1./sample_rate)) + rffreq
	Tc = nchirp/sample_rate
	resonance *= np.exp(2j*pi*fs*Tc)
	zfc = np.abs(resonance)
	peak = np.max(zfc)
	zfc /= peak
	fs /= MHz
	fw = 0.8
	f0 = fr/MHz
	ok = (f0 - fw < fs) & (fs < f0 + fw)

	# Now test and see the relation between this resonance and the main one
	gamma1 = 6381 #f0/(3*Qr)
	A1 = np.max(zfc)*pi/gamma1
	A2 = A1/4
	gamma2 = 675 #0.5*gamma1
	p0 = [f0, gamma1, gamma2, A1, A2]
	popt, pcov = optimize.curve_fit(psd_model, fs[ok], zfc[ok], p0=p0,
			method = 'lm')
	f0, gamma1, gamma2, A1, A2 = popt
	popt2, pcov2 = optimize.curve_fit(psd_lorentzian, fs[ok],
			resonance[ok].real, p0=[f0, f0/2/Qr_fit, A1], method = 'lm')
	fr_realfit, gamma_realfit, A_realfit = popt2
	Qr_realfit = fr_realfit/2/gamma_realfit
	print (fr_realfit, Qr_realfit, A_realfit)
	popt3, pcov3 = optimize.curve_fit(psd_root_lorentzian, fs[ok],
			resonance[ok].real, p0=[f0, f0/2/Qr_fit, A1], method = 'lm')
	fr_sqrtfit, gamma_sqrtfit, A_sqrtfit = popt2
	Qr_sqrtfit = fr_sqrtfit/2/gamma_sqrtfit
	print (fr_sqrtfit, Qr_sqrtfit, A_sqrtfit)
	#result = optimize.minimize(chisq, p0, args=(fs[ok], zfc[ok]))
	#f0, gamma, A = result['x']
	Qr1 = f0/(2*gamma1)
	Qr2 = f0/(2*gamma2)
	print (f0, Qr1, Qr2, A1, A2)
	#chisq_val = chisq(result["x"], fs[ok], zfc[ok])
	fs_fit = np.linspace(fs[ok][0], fs[ok][-1], 1000)

	# my model for the fast chirp resonance is the analytic signal of the
	# resonator.
	f_start = rffreq + sample_rate/2
	f_end = rffreq - sample_rate/2
	k = 1./nchirp
	d_sig = 2*pi*np.abs(f_end - f_start)
	tc = nchirp*sample_rate
	H = S21(fs*MHz, *reso_args)
	H11 = H - 1
	z_trial = -0.5*H11 - 0.5j*np.convolve(H11, 1/fs, mode='same')/pi
	z_trial /= np.max(np.abs(z_trial))
	#z_trial *= -1j

	plt.figure()
	plt.plot(fs[ok], zfc[ok], 'k', label='Sim Abs')
	plt.plot(fs[ok], resonance[ok].real/peak, 'b', label='Sim Real')
	plt.plot(fs[ok], resonance[ok].imag/peak, 'r', label='Sim Imag')
	#plt.plot(fs[ok], np.abs(z_trial)[ok], 'm', label='Trial Abs')
	#plt.plot(fs[ok], z_trial[ok].real, 'c', label='Trial Real')
	#plt.plot(fs[ok], z_trial[ok].imag, 'k', label='Trial Imag')
	#plt.plot(fs[ok], amplitude, 'r', label='Fit')
	#plt.plot(fs[ok], psd_model(fs[ok], *popt), 'r', label='Fit')
	plt.legend()
	plt.grid()
	plt.show()

	plt.figure()
	plt.plot(fs[ok], zfc[ok], 'k', label='Data')
	plt.plot(fs[ok], resonance[ok].real/peak, 'b', label='Data Real')
	plt.plot(fs[ok], resonance[ok].imag/peak, 'r.', label='Data Imag')
	plt.plot(fs[ok], psd_model(fs[ok], *popt), 'm--', label='Fit')
	plt.legend()
	plt.grid()
	plt.show()


if __name__=="__main__":
	navg = 100
	ntotal = 16384 #8192
	nchirp = 3400 #1700
	sample_rate = 100e6
	fast_chirp_rate = sample_rate / ntotal
	chirp_rate = sample_rate / navg
	z0 = 1<<14

	rffreq = 450e6
	Qr = 30000
	Qe = 40000
	#Qe *= (1 + 0.1j)
	fr = 450e6
	get_chirpspectrum(nchirp, ntotal, sample_rate, rffreq, fr, 1./Qr, 1./Qe)
















