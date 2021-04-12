#! /usr/bin/env python3

import numpy as np
from scipy import optimize
from math import pi
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.integrate import quad
from scipy.constants import h,k,c,e

um = 1e-6
eV = e
mK = 1e-3

N0 = 1.72e10 * um**-3 * eV**-1
Tc = 1.310#*K
alphak = 0.37

t = 0.05 * um
w_trace = 1 * um #width of inductor trace
s_trace = 1 * um #spacing between inductor meanders
N_sq = 16200  # Number of squares
l = N_sq * w_trace # Total length of the inductor (approximately. More exact is +21um + ...)
A = t * w_trace # cross-sectional area
V_sc = l * A


def delta_calc(T, Delta0):
    beta = 1./(k*T)
    Delta = Delta0
    gamma = beta*Delta/2.

    def integrand(x):
        y = np.sqrt(1+x**2)
        return (1./y)*np.exp(-2*gamma*y)/(1 + np.exp(-2*gamma*y))

    for i in range(10):
        print (Delta/(1.763*k))
        gamma = beta*Delta/2.
        integral, _ = quad(integrand, 0, np.infty)
        Delta = Delta0*(1 - integral)

    return Delta

def nqp_thermal(T, Delta):
    eta = (Delta/(k*T))
    nqpth = 2*N0*np.sqrt(2*pi*k*T*Delta)*np.exp(-eta)
    return nqpth


# Give T in Kelvin and E in Joules
def thermal_dist(E, T):
    alpha = (E/k/T)
    return 1./(np.exp(alpha)+1)


# Give T in Kelvin and Delta in Joules
def nqp_thermal_integrated(T, Delta):
    nqp = np.zeros_like(T)
    for i in range(T.size):
        eta = (Delta/(k*T[i]))
        integrand = lambda x: np.exp(-eta*np.sqrt(1+x**2))/(1 + np.exp(-eta*np.sqrt(1+x**2)))
        y, yerr = quad(integrand, 0, np.infty)
        nqp[i] = y

    nqp *= 4*N0*Delta
    return nqp

Delta0 = 1.763*k*Tc
Delta = delta_calc(Tc, Delta0)

exit()


Delta = 1.763*k*Tc
T = np.r_[80:700:1000j]*mK


nqp_approx = nqp_thermal(T, Delta)
nqp_integral = nqp_thermal_integrated(T, Delta)

plt.figure(figsize=(10,10))
plt.plot(T/mK, nqp_approx*um**3, 'b', label='Approx')
plt.plot(T/mK, nqp_integral*um**3, 'r', label='Integrated')
plt.grid()
plt.legend(loc='upper left')
plt.xlabel('T [mK]')
plt.ylabel('Quasiparticle Density [um^-3]')
plt.show()

#def xQMB(self,temps,f0,Tc,alphak, nfloor=0, nstar=nqp_star):
#    T = temps*unit.K
#    Tc = Tc*unit.K
#    f0 = f0*1e6*unit.Hz
#    Delta = 1.76*kB*Tc
#    q = (h*f0/(2*kB*T)).to('')
#    #print (nfloor)
#    nfloor *= unit.um**-3
#    nstar *= unit.um**-3
#    #gamma = nfloor**2/(2*nstar*taumax)
#    nqpth = 2*N0*np.sqrt(2*pi*kB*T*Delta)*np.exp(-(Delta/(kB*T)).m)
#    nqpth = nqpth.to(unit.um**-3)
#    #nqp = nqpth
#    nqp = np.sqrt((nqpth + nstar)**2 + nfloor**2) - nstar
#    #nqp = np.sqrt(nqpth**2 + nfloor**2)
#    nqp = nqp.to(unit.um**-3)
#    lwa = (pi*(4*N0*Delta/nqp)**2).to('').m
#    Teffs = ((2*Delta/kB)/special.lambertw(lwa)).real.to(unit.kelvin)
#    hf2kt = (h * f0 / (2 * kB * Teffs)).to('')
#    S1t = ((2./pi) * np.sqrt(2 * Delta / (pi * kB * Teffs)) * np.sinh(hf2kt) * special.k0(hf2kt)).to('')
#    S2t = (1.0 +  np.sqrt(2 * Delta / (pi * kB * Teffs)) * np.exp(-hf2kt) * special.i0(hf2kt)).to('')
#    QMB = (2*N0*Delta / (alphak*S1t*nqp)).to('').m
#    #S1 = (2./pi)*np.sqrt(2*Delta/(pi*kB*T))*np.sinh(q)*special.k0(q)
#    #S2 = 1. + np.sqrt(2*Delta/(pi*kB*T))*np.exp(-q)*special.i0(q)
#    #print("S2/S1: ",S2/S1)
#    xMB = (-alphak * S2t * nqp / (4.*N0*Delta)).to('').m
#    #xMB = xMB.to('')
#    #xMB = xMB.m
#    ##QMB = 2*N0*Delta/(alphak*S1*nqp)
#    #QMB = QMB.to('')
#    #QMB = QMB.m
#    return xMB,QMB

