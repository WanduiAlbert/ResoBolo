#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from math import pi


def make_switching_noise(t, fs, T, dx, sigma, phi0):
    y = dx*signal.square(2*pi*t/T, duty=0.5)
    roll_dist = int((phi0/2/pi)*fs*T)
    np.roll(y, roll_dist)
    y += sigma*np.random.randn(t.size)
    return y


def main():
    fs = 2e3
    T = 6
    N = 10e6
    t = np.arange(N, dtype=np.float128)/fs
    dx = 10
    sigma = 2
    phi0 = 0.1
    x = make_switching_noise(t, fs, T, dx,sigma, phi0)

    plt.figure(figsize=(10,10))
    plt.plot(t[:100000], x[:100000])
    plt.grid()
    plt.xlabel('Time [s]')
    plt.ylabel('df [Hz]')
    plt.savefig('switching_noise_timestream.png')
    plt.show()

    f, psd = signal.welch(x,fs=fs,window='hann', nperseg=1<<18)
    asd = np.sqrt(psd)
    asd[0] = np.nan

    plt.figure(figsize=(10,10))
    plt.loglog(f, asd)
    plt.grid()
    plt.xlim(left=1e-2)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('ASD [Hz/rtHz]')
    plt.savefig('switching_noise_spectrum.png')
    plt.show()

if __name__=="__main__":
    main()
