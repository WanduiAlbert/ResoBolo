#! /usr/bin/env python3

"""
heatersweepsim.py
-------------------------------------------------------------------------------

Simulating the response of the TKID bolometer to different drive powers on the
heater. Useful for predicting what I should see when doing time constant type
measurements.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate, optimize
from math import pi

def get_periodic_ts(nperiod, ntotal, sample_rate):
	ntile = (ntotal // nperiod) + 1

	ts = np.ones((ntile, nperiod)) * np.arange(nperiod) / sample_rate
	return ts.flatten()[:ntotal]


def squaresweep(Tacq, *args):
	f1, f2, T, phi0, sample_rate = args
	ntotal = int(Tacq*sample_rate)
	t = np.arange(ntotal) / sample_rate
	mask = (t // (T//2))
	print (mask)
	y = np.ones_like(t)
	y[mask % 2 == 1] = -1.0

	return t, y

def boloresponse(f, f3dB):
	return 1./(1 + 1j*(f/f3dB))

def logsweep(Tacq, *args):
	f1, f2, T, phi0, sample_rate = args
	k = np.log(f2/f1)
	nperiod = int(T*sample_rate)
	ntotal = int(Tacq*sample_rate)
	t = get_periodic_ts(nperiod, ntotal, sample_rate)
	phase = 2*pi*f1*T/k*(np.exp(t*k/T) - 1) + phi0
	return np.arange(ntotal)/sample_rate, np.sin(phase)

def getsweepasd(Tacq, func, *args):
	sample_rate = args[-1]
	t, y = func(Tacq, *args)
	nt = len(t)

	asd = np.abs(np.fft.fft(y))[:nt//2]
	fs = np.fft.fftfreq(nt, d=1./sample_rate)[:nt//2]

	return fs, asd

def get_resonator_ts(Tacq, f3dB, func, *args):
	sample_rate = args[-1]
	t, y = func(Tacq, *args)
	ntotal = len(t)
	fs = np.fft.fftfreq(ntotal, d=1./sample_rate)
	ys = np.fft.fft(y) * boloresponse(fs, f3dB)
	noise  = np.random.randn(ntotal) + 1j*np.random.randn(ntotal)
	noise *= 5
	ys += noise
	# Invert the fourier transform
	yout = np.fft.ifft(ys).real
	return t, yout

def thermal_model(f,A,f3db):
	return np.log10(A / np.sqrt(1.0+(f/f3db)**2))

f1 = 1
f2 = 100
T =  10
Tacq = 120
f3dB = 28
phi0 = np.random.randn()
sample_rate = 100e6/100./8192 # get the fastchirp sampling rate
print ("sample rate: %1.3f Hz"%sample_rate)
print ("start freq: %1.3f Hz"%f1)
print ("stop freq: %1.3f Hz"%f2)
print ("sweep time: %1.3f s"%T)
print ("phase offset as fraction of 2 pi: %1.3f"%(phi0/2/pi))
print ("total time: %1.3f s"%Tacq)

#t, y = logsweep(Tacq, f1, f2, T, phi0, sample_rate)
#fs, asd = getsweepasd(T, logsweep, f1, f2, T, phi0, sample_rate)
t, y = get_resonator_ts(Tacq, f3dB, logsweep, f1, f2, T, phi0, sample_rate)
#t, y = squaresweep(Tacq, f1, f2, T, phi0, sample_rate)
#fs, asd = getsweepasd(T, squaresweep, f1, f2, T, phi0, sample_rate)

#plt.figure(1)
#plt.plot(t, y)
#plt.show()

#plt.figure(2)
#plt.plot(fs, asd)
#plt.show()

nperseg = 1024
fs, asd = signal.welch(y, fs=sample_rate, detrend='linear', scaling="density",
		nperseg=nperseg)
t_ref, y_ref = logsweep(Tacq, f1, f2, T, phi0, sample_rate)
fs_ref, asd_ref = signal.welch(y_ref, fs=sample_rate, detrend='linear',
		scaling="density", nperseg=nperseg)
asd /= asd_ref

ly = np.log10(asd)
p0 = (asd[0],10.0)
#thermal_model(f,A,f3db)
popt,pcov = optimize.curve_fit(thermal_model,fs,ly,p0=p0)
asd_fit = 10**thermal_model(fs,*popt)
print("f3db: ",popt[1])
print("A: ",popt[0])
plt.figure(123)
#plt.title('Bolometer Frequency Response at P = 5.00pW')
plt.semilogx(fs, asd,label="Measured")
plt.semilogx(fs, asd_fit,label='Fit f3db=%.1fHz'%popt[1])
plt.grid()
plt.ylim([0, 1.2])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Transfer Function')
plt.legend(loc='upper left')
#plt.savefig('fig/reso_%3.1fMHz_timeconstant.png'%resfreqs[ch])
plt.show()

