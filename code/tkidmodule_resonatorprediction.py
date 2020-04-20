#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import pi
from scipy import interpolate
from scipy.constants import h, k, c
from scipy import interpolate, optimize
import scipy.signal as signal
import glob
import sys,os
import pdb
from scipy.constants import epsilon_0, mu_0
import reso_fit
import meshanalysisofresonator as mesh
pi = np.pi

MHz = 1e6
GHz = 1e9
kHz = 1e3
pW = 1e-12
um = 1e-6

nH = 1e-9
pF = 1e-12
MHz = 1e6
Z0 = 50
Y0 = 1./Z0

datadir = '../numerical_sims/'
plotdir = 'TKIDModule_fig/'

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'][1:]

plot_diagnostic = True
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

N = np.array([4, 6, 8, 10, 12, 14, 16])*53
Qcsim = np.array([24101, 29710, 35751, 41224, 46913, 49695, 56001])
frsim = np.array([660.4, 551.6, 480.6, 429.3, 390.1, 358.8, 334.0])*MHz

C = np.array([3.7593, 5.6224, 7.4950, 9.2137, 10.9582, 12.7016, 14.3668])*pF
Cc = np.array([0.09253, 0.1155, 0.1320, 0.1519, 0.1712, 0.1960, 0.2155])*pF
Ca = np.array([2.1492, 2.1581, 2.1907, 2.0866, 2.04726, 2.00413, 1.9076])*pF
Cb = np.array([1.3837, 1.3840, 1.4358, 1.4695, 1.5215, 1.5818, 1.6951])*pF
Cg = Ca - Cc
Ct = Ca + Cb + Cc

Lpar = np.array([2.5685, 2.7041, 2.9943, 3.4146, 3.8231, 4.2223, 4.6096])*nH
Lind = 10.6964*nH
Ltot = Lind + Lpar

w0 = 1./np.sqrt(Ltot*C)
u = 4*Ca*(Cb+Cc)*(Ca + Cb + Cc)
v = Ca*Cb*(Ca+Cb)*Cc**2*Z0**2*w0**2
x = 4*(Ca + Cb + Cc)**2
y = (Ca+Cb)**2*Cc**2*Z0**2*w0**2
#print ("Ratio of v/u ", v/u)
#print ("Ratio of y/x ", y/x)


num = 4*Ca*(Cb + Cc)*(Ca + Cb + Cc)
denom = 4*(Ca + Cb + Cc)**2
Ccorr = num/denom
print ("Additional capacitive loading ", Ccorr/pF)

fr = 1./np.sqrt(Ltot*C)/2/pi
print ("Frequency with no capacitive loading ", fr/MHz)
fr = 1./np.sqrt(Ltot*(C+Ccorr))/2/pi
print ("Frequency with capacitive loading ", fr/MHz)
Qc = 2*C/(Z0*2*pi*fr*Cc**2)
print ("Initial term in Qc ", Qc)
Qc *= (Ct/Ca)**2
print ("Second term in Qc ", Qc)
Qc *= (Ltot/Lind)**2
print ("Third term in Qc ", Qc)

print ("Qc from simulation vs. predicted ", Qcsim/Qc)
print ("fr from simulation vs. predicted ", (frsim - fr)/MHz)

plt.figure()
plt.plot(N, frsim/MHz, 'rs', ms=12, label='From sonnet simulation')
plt.plot(N, fr/MHz, 'bd', ms=12,  label='Best prediction')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Frequency [MHz]')
plt.legend(loc='upper right')
plt.savefig('sonnet_vs_predicted_fr.png')
plt.show()

plt.figure()
plt.plot(N, Qcsim, 'rs', ms=12,  label='From sonnet simulation')
plt.plot(N, Qc, 'bd', ms=12,  label='Best prediction')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Qc')
plt.legend(loc='upper left')
plt.savefig('sonnet_vs_predicted_Qc.png')
plt.show()






