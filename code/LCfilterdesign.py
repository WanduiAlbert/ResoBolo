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
    R1 = 50 #Ohms. The heater load resistance. Unclear if this is the Z0 to
               # which I should design instead of 50 Ohms.
    R2 = 0.1 # The resistance of the heaters.
    wc = 500*MHz
    n = 3 # Order of the Butterworth filter
    w1 = wc
    L2 = R2/wc
    C = 1/w1/R1
    L1 = 1/C/w1**2
    print ("L1 = %1.3fnH, C = %1.3fpF, L2 = %1.3fnH"%(L1/nH, C/pF, L2/nH))

    # Want to see the response expected for a filter such as this
    omega = np.linspace(0.1, 2*pi*800, 1000)*MHz
    L2 = 5*nH
    #R2 = 50
    H = R2/(-1j*omega**3*L1*L2*C - omega**2*C*(L2*R1 + L1*R2) + 1j*omega*(L2 +
        R1*R2*C) + R2)
    HdB = 20*np.log10(np.abs(H))

    plt.figure()
    plt.plot(omega/2/pi/MHz,  HdB)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('|H(f)| [dB]')
    plt.grid()
    plt.show()
