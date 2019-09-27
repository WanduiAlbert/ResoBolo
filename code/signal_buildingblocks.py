#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy import misc

x = np.zeros(30)
x[0] = 1


m1 = -1.8
m2 = 0.96

# Low pass
b = np.array([1,2,1])
a = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("Low Pass")
ax.axis('tight')
plt.savefig('standard_lowpass_response.png')
#plt.show()


# Band pass
b = np.array([-1,0,1])
a = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("Band Pass")
ax.axis('tight')
plt.savefig('standard_bandpass_response.png')
#plt.show()

# High pass
b = np.array([1,-2,1])
a = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("High Pass")
ax.axis('tight')
plt.savefig('standard_highpass_response.png')
#plt.show()

# Notch
b = np.array([1,m1/m2**0.5,1])
a = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("Notch")
ax.axis('tight')
plt.savefig('standard_notch_response.png')
#plt.show()

m3 = -1.42

# Low Pass Notch
b = np.array([1,m3,1])
a = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("Low Pass Notch")
ax.axis('tight')
plt.savefig('standard_lowpass_notch_response.png')

# High Pass Notch
b = np.array([1,m3,1])
a = np.array([m2, -m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.set_title("High Pass Notch")
ax.axis('tight')
plt.savefig('standard_highpass_notch_response.png')

# All Pass
a = np.array([1,m1,m2])
b = np.array([m2, m1, 1])
w, h = sig.freqz(b, a)
h /= np.max(np.abs(h))
h = 20*np.log10(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(w, h)
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response [dB]')
ax.axis('tight')
ax.set_ylim((-10,10))
ax.grid()
ax.set_title("All Pass")
plt.savefig('standard_allpass_response.png')

#plt.show()

# Digital Oscillators
w0 = 7./10*pi
b = np.array([0,np.sin(w0),0])
a = np.array([1, -2*np.cos(w0), 1])
x = np.zeros(100)
x[0] = 1
y = sig.lfilter(b, a, x)
w, h = sig.freqz(b, a)
#h /= np.max(np.abs(h))
fig, ax = plt.subplots(figsize=(10,10))
ax.plot(y)
ax.set_xlabel('Iteration')
ax.set_ylabel('Output Signal')
ax.grid()
ax.set_title("Digital Oscillators")
ax.axis('tight')
plt.savefig('standard_digital_oscillators_response.png')

plt.show()












