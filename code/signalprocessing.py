#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy import misc 

x = np.zeros(30)
x[0] = 1

fig, ax = plt.subplots(figsize=(10,10))
fig, ax2 = plt.subplots(figsize=(10,10))
for M in [1, 2, 3, 6, 10, 20]:

    b = np.ones(M)/M
    a = np.array([1])

    z, p, k = sig.tf2zpk(b, a)
    print (z)
    print (p)
    print (k)

    w, h = sig.freqz(b, a)
    index = np.arange(w.size)
    delay = -np.unwrap(np.angle(h), discont=pi)/w
    ax.plot(w/pi, np.abs(h), label="M = %d"%M)
    ax2.plot(w/pi, delay, label="M = %d"%M)
    ax2.hlines((M-1)/2, 0, 1, linestyles='dashed')
ax.set_xlabel('Normalized Frequency [rad/sample]')
ax.set_ylabel('Magnitude Response')
ax.grid()
ax.legend(loc='best')
ax.axis('tight')
plt.savefig('convolution_lowpass_magnituderesponse.png')
ax2.set_xlabel('Normalized Frequency [rad/sample]')
ax2.set_ylabel('Phase Delay [Samples]')
ax2.grid()
ax2.legend(loc='best')
ax2.axis('tight')
plt.savefig('convolution_lowpass_phaseresponse.png')
plt.show()
