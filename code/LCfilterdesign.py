#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi

""""
Designing a low pass filter to ensure that the input into the heaters is cleaned
of RF interference that leaks into our resonators. The goal is to try and make
the cut off frequency as low as possible at the expense of making the filters
incredibly large.
"""

nH = 1e-9
pF = 1e-12
MHz = 1e6

if __name__=="__main__":
    Z0 = 50 #Ohms. The heater load resistance. Unclear if this is the Z0 to
               # which I should design instead of 50 Ohms.
    wc = 100*MHz
    n = 3 # Order of the Butterworth filter
    k = np.arange(1, n+1)
    g = np.abs(2*np.sin((2*k-1)*pi/2/n))
    L = g[0] * Z0/wc
    C = g[1] * 1/Z0/wc

    print ("L = %1.3fnH, C = %1.3fpF"%(L/nH, C/pF))

    # Want to see the response expected for a filter such as this
    omega = np.linspace(0.1, 2*pi*800, 1000)*MHz
    ZL = 50#0.106 # The resistance of the heaters.
    Z_filt = 1j*omega*L*(1 + 1/(1 - omega**2*L*C))
    H2 = ZL/(Z_filt + ZL) # Transfer function of the filter
    H = ZL/(-1j*omega**3*L**2*C - omega**2*L*C*ZL + 2j*omega*L + ZL)
    H2dB = 10*np.log10(np.abs(H2))
    HdB = 10*np.log10(np.abs(H))

    plt.figure()
    plt.semilogx(omega/2/pi/MHz,  HdB)
    plt.semilogx(omega/2/pi/MHz,  H2dB)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('|H(f)| [dB]')
    plt.grid()
    plt.show()
