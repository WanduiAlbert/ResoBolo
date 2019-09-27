#! /usr/bin/env python3

import numpy as np
from math import pi
import scipy.signal as sig
import matplotlib.pyplot as plt

N = 53
Omega_p = 4.
Omega_r = 4.2
Omega_s = 10.0

kp = int(N*(Omega_p/Omega_s))
kr = int(N*(Omega_r/Omega_s))

k = np.arange(1, N//2)
A = np.zeros(N)
A[:(kp+1)] = 1
h = np.zeros(N)

for n in range(N):
    h[n] = A[0] + 2*np.sum((-1)**k * A[k] * np.cos(pi*k*(1+2*n)/N))

h /= N
print (h[:(N//2+1)])

fig, ax = plt.subplots(figsize=(10,10))
ax.stem(h)
ax.set_xlabel('Index')
ax.set_ylabel('h[n]')
ax.grid()
plt.savefig('freqsampled_lpfilter.png')
#plt.show()

H = np.fft.rfft(h, n=512)
f = np.fft.rfftfreq(512, d=2*pi/Omega_s)
Omega = 2*pi*f
HdB = 20*np.log10(np.abs(H))

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(Omega, HdB)
ax.set_xlabel('Frequency [rad/s]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.axis('tight')
plt.savefig('freqsampled_lpfilter_magnituderesponse.png')
#plt.show()

plt.close('all')

# Bandstop filter
N = 81
M = N
Omega_c1 = 2000
Omega_c2 = 4000
Omega_s = 10000

wc1 = 2*pi*(Omega_c1/Omega_s)
wc2 = 2*pi*(Omega_c2/Omega_s)

n = np.arange(1, N//2+1)
h = np.zeros(N)
h0 = 1 - (wc2 - wc1)/pi
haux = (np.sin(wc1*n) - np.sin(wc2*n))/(pi*n)

h[:N//2] = haux[::-1]
h[N//2] = h0
h[(N//2+1):] = haux
print (h[:(N//2+1)])
n = np.arange(N) - N//2
#h = np.stack([np.fliplr(haux), [h0], haux])

# Rectangular Window
wr = np.zeros(N)
mask = np.abs(n) < M//2
wr[mask] = 1

window = wr

#Triangular window
wt = np.zeros(N)
mask = np.abs(n) < M//2
wt[mask] = -2*np.abs(n[mask])/(M+2) + 1
#window = wt

#Bartlett window
wtB = np.zeros(N)
mask = np.abs(n) < M//2
wtB[mask] = -2*np.abs(n[mask])/M + 1
#window = wtB

#Hamming window
alpha = 0.54
wH = np.zeros(M)
mask = np.abs(n) < M//2
wH[mask] = alpha + (1-alpha)*np.cos(2*pi*n[mask]/M)
window = wH

#Hann window
alpha = 0.5
wH = np.zeros(M)
mask = np.abs(n) < M//2
wH[mask] = alpha + (1-alpha)*np.cos(2*pi*n[mask]/M)
#window = wH

#Blackman window
wB = np.zeros(M)
mask = np.abs(n) < M//2
wB[mask] = 0.42 + 0.5*np.cos(2*pi*n[mask]/M) + 0.08*np.cos(4*pi*n[mask]/M)
#window = wB

h *= window
fig, ax = plt.subplots(figsize=(10,10))
ax.stem(h)
ax.set_xlabel('Index')
ax.set_ylabel('h[n]')
ax.grid()
plt.savefig('freqsampled_bandstopfilter.png')
#plt.show()

H = np.fft.rfft(h, n=512)
f = np.fft.rfftfreq(512, d=2*pi/Omega_s)
Omega = 2*pi*f
HdB = 20*np.log10(np.abs(H))

fig, ax = plt.subplots(figsize=(10,10))
ax.plot(Omega, HdB)
ax.set_xlabel('Frequency [rad/s]')
ax.set_ylabel('Magnitude Response [dB]')
ax.grid()
ax.axis('tight')
plt.savefig('freqsampled_bandstopfilter_magnituderesponse.png')
plt.show()
