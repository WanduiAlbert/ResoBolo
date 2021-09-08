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
plotdir = 'TKIDModule_fig_2um/'

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'][1:]

plot_diagnostic = True
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

N = np.array([4, 6, 8, 10, 12, 14, 16])*53
N0 = 600
scaling = 1
#scaling = 1.4789
#Qcsim = np.array([24101, 29710, 35751, 41224, 46913, 49695, 56001])*scaling**2
#frsim = np.array([660.4, 551.6, 480.6, 429.3, 390.1, 358.8, 334.0])*MHz*scaling
Qcsim = np.array([17028, 22825, 29589, 36795, 45002, 51000, 61559])
frsim = np.array([990.2, 809.3, 692.2, 608.4, 544.8, 494.4, 454.7])*MHz

C = np.array([3.7593, 5.6224, 7.4950, 9.2137, 10.9582, 12.7016, 14.3668])*pF
Cc = np.array([0.09253, 0.1155, 0.1320, 0.1519, 0.1712, 0.1960, 0.2155])*pF
Ca = np.array([2.1492, 2.1581, 2.1907, 2.0866, 2.04726, 2.00413, 1.9076])*pF
Cb = np.array([1.3837, 1.3840, 1.4358, 1.4695, 1.5215, 1.5818, 1.6951])*pF
Cg = Ca - Cc
Ct = Ca + Cb + Cc

Lpar = np.array([2.5685, 2.7041, 2.9943, 3.4146, 3.8231, 4.2223, 4.6096])*nH
#Lind = 10.6964*nH
Lind = 4.8905*nH
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
print ("fr from simulation vs. predicted ", (frsim/fr))

plt.figure()
plt.plot(N, frsim/MHz, 'rs', ms=12, label='From sonnet simulation')
plt.plot(N, fr/MHz, 'bd', ms=12,  label='Best prediction')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Frequency [MHz]')
plt.legend(loc='upper right')
plt.savefig(plotdir + 'sonnet_vs_predicted_fr.png')
plt.show()

plt.figure()
plt.plot(N, Qcsim, 'rs', ms=12,  label='From sonnet simulation')
plt.plot(N, Qc, 'bd', ms=12,  label='Best prediction')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Qc')
plt.legend(loc='upper left')
plt.savefig(plotdir + 'sonnet_vs_predicted_Qc.png')
plt.show()

plt.figure()
plt.plot(N, Qcsim/Qc, 'rs', ms=12)
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Qcsim/Qc')
plt.savefig(plotdir + 'sonnet_vs_predicted_Qc_ratio.png')
plt.show()

#def freqmodel(x, fa, alpha, beta, epsilon):
#    return fa/np.sqrt(x)/np.sqrt(1 + alpha + beta*x + epsilon/x)

def freqmodel(x, fa, a, b):
    return fa/np.sqrt(x)/np.sqrt(1 + a*x + b/x)


#def freqmodel(x, fa, a1, a2, b1):
#    alpha = a1*b1
#    beta = a1 + a2*b1
#    gamma = a2
#    epsilon = b1
#    return fa/np.sqrt(x)/np.sqrt(1 + alpha + beta*x + gamma*x**2 + epsilon/x)

p0 = [695., 1.5e-02, 7.4e-03]
xfit = np.r_[0.3:1.6:1000j]

popt, pcov = optimize.curve_fit(freqmodel, N/N0, frsim/MHz, p0=p0)
ypred = freqmodel(xfit, *popt)
print (popt)
res  = frsim/MHz - freqmodel(N/N0, *popt)

plt.figure()
plt.plot(N, frsim/MHz, 'ko', ms=12, label='From sonnet simulation')
plt.plot(xfit*N0, ypred, 'r', ms=12,  label='Fit')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Frequency [MHz]')
plt.legend(loc='upper right')
plt.savefig(plotdir + 'frequency_vs_nfingers.png')
plt.show()

plt.figure()
plt.plot(N, res, 'ko', ms=12, label='From sonnet simulation')
plt.grid()
plt.xlabel('Nfingers')
plt.ylabel('Residuals [MHz]')
plt.legend(loc='upper right')
plt.savefig(plotdir + 'residuals_vs_Nfingers.png')
plt.show()
exit()

print (popt)
print (pcov)
