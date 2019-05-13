#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk
from math import pi
from scipy.constants import mu_0, epsilon_0

um = 1e-6
nH = 1e-9
pF = 1e-12

w = 1.5 * um
g = 0.5 * um
h = 307 * um
l = 2 * (w + g)
N = 126 * um /l

eps = 2*w/l
ki = np.sin(pi*eps/2)
kip = np.sqrt(1 - ki**2)

ke = 2*eps**0.5/(1+eps)
kep = np.sqrt(1 - ke**2)

Li = mu_0 * h * ellipk(kip)/ellipk(ki)
Le = mu_0 * h * ellipk(kep)/ellipk(ke)

print (Li/nH)
print (Le/nH)

Lt = (N-3)*Li/2 + (Li*Le)/(Li + Le)
print (Lt/nH)
