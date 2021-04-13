#! /usr/bin/env python3

import numpy as np
from scipy import optimize
from math import pi
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.integrate import quad
from scipy.constants import h,k,c,e
import scipy.optimize as opt

um = 1e-6
eV = e
mK = 1e-3

N0 = 1.72e10 * um**-3 * eV**-1
Tc = 1.310#*K
alphak = 0.37

T_debye = 433#K
Omega_debye = k*T_debye


Delta0 = 1.763*k*Tc
dNV = np.log(2*Omega_debye/Delta0)

t = 0.05 * um
w_trace = 1 * um #width of inductor trace
s_trace = 1 * um #spacing between inductor meanders
N_sq = 16200  # Number of squares
l = N_sq * w_trace # Total length of the inductor (approximately. More exact is +21um + ...)
A = t * w_trace # cross-sectional area
V_sc = l * A

sech = lambda x: np.sqrt(1-np.tanh(x)**2)

def delta_calc(T, Delta0):
    beta = 1./(k*T)
    x_start = beta*Delta0/2
    x_upper = beta*Omega_debye/2


    def integrand(x, x0):
        y = np.sqrt(x**2 + x0**2)
        return np.tanh(y)/y

    def deriv_integrand(x, x0):
        y = np.sqrt(x**2 + x0**2)
        return x0*(sech(y)**2 - np.tanh(y)/y)/y**2

    # Return the scalar function and its derivative
    def gapfactor(xdelta):
        integral, _ = quad(integrand, 0, x_upper, args=(xdelta,))
        dintegral, _ = quad(deriv_integrand, 0, x_upper, args=(xdelta,))
        return integral - dNV, dintegral

    #zerolevel, _ = quad(lambda x: 1.0 if x ==0 else np.tanh(x)/x, 0, x_upper)
    #zerolevel, _ = quad(integrand, 0, x_upper, args=(0,))
    #print (zerolevel/dNV)
    #exit()
    #xd = np.r_[0:x_upper:100j]
    #y = np.zeros_like(xd)
    #for i in range(xd.size):
    #    f, fprime = gapfactor(xd[i])
    #    y[i] = f

    #plt.figure(figsize=(10,10))
    #plt.plot(xd, y)
    #plt.axhline(dNV)
    #plt.grid()
    #plt.show()

    #exit()
    #xsol = x_start
    #for i in range(20):
    #    print (xsol)
    #    f, fprime = gapfactor(xsol)
    #    xsol -= f/fprime
    sol = optimize.root_scalar(gapfactor, x0=x_start, bracket=[0, x_upper],
            fprime=True, method='newton')

    xsol = sol.root

    Delta = 2*np.abs(xsol)/beta
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

T = np.r_[0.01:Tc:1000j]

Delta = np.zeros_like(T)
for i in range(T.size):
    Delta[i] = delta_calc(T[i], Delta0)

plt.figure(figsize=(10,10))
plt.plot(T/Tc, Delta/Delta0)
plt.grid()
plt.xlabel('$T/T_c$')
plt.ylabel('$\Delta/\Delta_0$')
plt.show()


exit()




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

