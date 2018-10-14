#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi

MHz = 1e6


# The transfer function of a single resonator
def S21(f, fr, dQr, dQe):
    x = (f - fr)/fr
    return 1 - dQe/(dQr + 1j * 2 * x)

# The time domain waveform describing the chirp waveform
def zin(t, z0, f0, k, phi0):
    exponent = phi0 + 2*pi*f0*t + pi*k*t**2
    return z0*np.exp(1j*exponent)

if __name__=="__main__":
    navg = 100
    ntotal = 8192
    nchirp = 1700
    sample_rate = 100e3
    fast_chirp_rate = sample_rate/ntotal
    chirp_rate = fast_chirp_rate / navg
    chirp_period = 1e3 / chirp_rate

    nsamps = int(chirp_period*sample_rate)

    dt = 1./sample_rate
    t = nsamps * np.arange(nsamps)

    rffreq = 450e6
    f_start = rffreq - sample_rate / 2.0
    f_end = rffreq + sample_rate / 2.0

    k = (f_end - f_start) / (sample_rate * nchirp)

    # Use the full scale ADC output
    z0 = 1<<14
    z_chirp = zin(t, z0, f_start, k, 0)

    # plt.plot(t, z_chirp.real)
    # plt.show() 
