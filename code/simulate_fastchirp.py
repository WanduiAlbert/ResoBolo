#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi

np.random.seed(100)

MHz = 1e6
Z0 = 50

calibration = 5.775091506536823e-6  # Vrms per ADC unit

attenuation = 75.2 # dB. Sum of the warm & cold attenuation + VAUnix IL.
usrp_power = -6 #dBm
amplifier_gain = 53 #dB

# Goal is to simulate 120s fast chirp noise timestreams starting from a
# projected sum of two lorentzians that we see in the data.

def lorentzian (x, x0, gamma):
	return (gamma/pi)*1./((x-x0)**2 + gamma**2)

def generate_spectrum(f0, BW, N, Pg):
	"""
	generate_lorentzian(f0, Pg)

	Parameters: f0 - Resonance frequency in MHz
                BW - Bandwidth in MHz
                N  - Number of points
				Pg - Tx gain of the USRP

	returns     f  - Frequency of the output
                L  - Output spectrum in ADC^2/Hz

    Generates a sum of two lorentzians as a model for the fast chirp response.
	The parameters are chosen at random from the unit uniform interval.

	"""
	A1 = np.random.uniform()
	A2 = np.random.uniform()

	gamma1 = np.random.uniform()
	gamma2 = np.random.uniform()

	res = A1 * lorentzian(f, f0, gamma1) + A2 * lorentzian(f, f0, gamma2)

	P_out = Pg + usrp_power - attenuation + amplifier_gain

	# Want the power in Watts instead of in dBm
	P = 1e-3 * 10**(P_out/10.)
    V = P**0.5 * Z0 # in Vrms units
    scale = V/calibration

    freq = np.linspace(f0 - BW/2., f0 + BW/2., N
    return freq, res * scale

