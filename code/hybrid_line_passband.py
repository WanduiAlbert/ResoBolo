#! /usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.constants import c
from math import pi

MHz = 1e6
GHz = 1e9
um = 1e-6

nu = np.linspace(1, 10e6, 1000)*MHz
er = 4
n = er**0.5
k = 2*pi*nu/n/c
d = 100e-6
eta = 0.94
argument = 0.5*np.cos((1-eta)*k*d)*(1 + np.cos(eta*k*d))
beta = 1/d * np.arccos(argument)


fig, ax = plt.subplots(figsize=(10,10))
ax.plot(nu/GHz, beta)
plt.show()
