#! /usr/bin/env python3

import numpy as np
from scipy import integrate
from math import pi
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0


def cap_integral_rect(x2, y2, y1, x1):
    return 1./np.sqrt((x2 - x1)**2 + (y2 - y1)**2)**3

def cap_integral_polar(theta2, theta1, r2, r1):
    return 1./np.sqrt(r1**2 + r2**2 - 2*r1*r2*np.cos(theta1-theta2))**3


# Model 1: An anular electrode of radius R with a gap distance of size s between
# the connector and the ground plane. I'll normalize distances with R = 1

N = 20
s = np.logspace(-3, 2, N)
R = 1

C = np.zeros(N)
for i in range(N):
    int_lims = [(0, 2*pi), (0, 2*pi), (R+s[i], np.infty), (0, R)]
    val, err = integrate.nquad(cap_integral_polar, int_lims)
    print (s[i], val, err)
    C[i] = val

C /= pi

plt.figure(figsize=(10,10))
plt.plot(s, C)
plt.grid()
plt.savefig('capacitance_of_an_anular_electode.png')
plt.show()

