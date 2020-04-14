#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0
from math import pi
from scipy import special
import mpmath
jtheta = np.vectorize(mpmath.jtheta, 'D')

um = 1e-6
pF = 1e-12
nH = 1e-9
MHz = 1e6
pH = 1e-12

Nlayers = 3
ers = [1, 3.9, 12]
hs = [np.inf, 0.3*um, np.inf]

def interior_capacitance(w, g, h, er):
    eta = w/(w + g)
    lambd = 2*(w + g)

    if h != np.inf:
        # For the finite layer dielectric
        r = h/lambd
        q = np.exp(-4*pi*r)
        #print (r)
        k = (jtheta(2, 0, q)/jtheta(3,0,q))**2
        k = k.real
        t4 = 1./k
        #print (special.ellipk(k)*eta)
        t2,_,_,_ = special.ellipj(special.ellipk(k)*eta, k)
        k1 = t2*np.sqrt((t4**2-1)/(t4**2-t2**2))
        k1p = np.sqrt(1-k1**2)
        Ca = epsilon_0*er*special.ellipk(k1)/special.ellipk(k1p)
    else:
        # For the infinite layer dielectric
        k1inf = np.sin(pi/2*eta)
        k1infp = np.sqrt(1-k1inf**2)
        Ca = epsilon_0*er*special.ellipk(k1inf)/special.ellipk(k1infp)

    #return Ca + Cb
    return Ca


def exterior_capacitance(w, g, h, er):
    eta = w/(w + g)
    lambd = 2*(w + g)

    if h != np.inf:
        r = h/lambd
        # For the finite layer dielectric
        t4 = np.cosh(pi*(eta+1)/(8*r))
        t3 = np.cosh(pi*(1-eta)/(8*r))
        kE = 1./t3*np.sqrt((t4**2-t3**2)/(t4**2-1))
        kEp = np.sqrt(1-kE**2)
        Ca = epsilon_0*er*special.ellipk(kE)/special.ellipk(kEp)
    else:
        # For the infinite layer dielectric
        kEinf = 2*eta**0.5/(1 + eta)
        kEinfp = np.sqrt(1-kEinf**2)
        Ca = epsilon_0*er*special.ellipk(kEinf)/special.ellipk(kEinfp)

    #return Ca + Cb
    return Ca

def capacitance(w, g, h, er, N, L):
    Cint = L*interior_capacitance(w,g,h,er)
    Cext = L*exterior_capacitance(w,g,h,er)
    C = (N-3)*Cint/2 + 2*(Cint*Cext)/(Cint + Cext)
    return C

def interior_inductance(w, g, h, er):
    eta = w/(w + g)
    lambd = 2*(w + g)
    r = h/lambd

    # For the infinite layer dielectric
    k1inf = np.sin(pi/2*eta)
    k1infp = np.sqrt(1-k1inf**2)
    Lb = mu_0*special.ellipk(k1infp)/special.ellipk(k1inf)

    return Lb

def exterior_inductance(w, g, h, er):
    eta = w/(w + g)
    lambd = 2*(w + g)
    r = h/lambd

    # For the infinite layer dielectric
    kEinf = 2*eta**0.5/(1 + eta)
    kEinfp = np.sqrt(1-kEinf**2)
    Lb = mu_0*special.ellipk(kEinfp)/special.ellipk(kEinf)

    return Lb

def inductance(w, g, h, er, N, L):
    Lint = L*interior_inductance(w,g,h,er)
    Lext = L*exterior_inductance(w,g,h,er)
    Lpar = (N+1)*Lint/2 + 2*Lext
    return Lpar

if __name__=="__main__":
    #W = 1000*um
    #L = 500*um
    #h = 500*um
    #N = 310
    W = 1238*um
    L = 140*um
    w = 1.0*um
    g = 1.0*um
    period = 2*(w+g)
    N = int(2*W/period)
    eta = w/(w+g)
    Cint = interior_capacitance(w,g,hs[0], ers[0])
    Cint += interior_capacitance(w,g,hs[1], ers[1] - ers[2])
    Cint += interior_capacitance(w,g,hs[2], ers[2])
    Cext = exterior_capacitance(w,g,hs[0], ers[0])
    Cext += exterior_capacitance(w,g,hs[1], ers[1] - ers[2])
    Cext += exterior_capacitance(w,g,hs[2], ers[2])
    Cint *= L
    Cext *= L
    C = (N-3)*Cint/2 + 2*(Cint*Cext)/(Cint + Cext)
    print ("The capacitance is %1.3fpF"%(C/pF))
    exit()

    Ls = 0.158*pH#/sq
    period = 8*um
    eta = np.r_[0:1:200j][1:-1]
    w = period*eta/2
    g = period/2 - w
    N = 502
    W = N*period/2
    #N = 2*int(W/period)
    er = 11.9
    Nsq = N*L/w

    Cs = capacitance(w, g, h, er, N, L)
    Lpar = inductance(w, g, h, er, N, L)
    L = 10*nH
    fr = 1./(2*pi*np.sqrt(L*Cs))

    fig,ax = plt.subplots(figsize=(10,10))
    ax.plot(eta, Cs/pF)
    ax.axvline(0.5, color='k', ls='--')
    ax.set_xlabel('Fill fraction')
    ax.set_ylabel('C [pF]')
    ax.grid()
    plt.show()

    #fig,ax = plt.subplots(figsize=(10,10))
    #ax.plot(eta, Lpar/nH)
    #ax.axvline(0.5, color='k', ls='--')
    #ax.set_xlabel('Fill fraction')
    #ax.set_ylabel('Lpar [nH]')
    #ax.grid()
    #plt.show()

    fig,ax = plt.subplots(figsize=(10,10))
    ax.plot(eta, fr/MHz)
    ax.axvline(0.5, color='k', ls='--')
    ax.set_xlabel('Fill fraction')
    ax.set_ylabel('Resonance Frequency [MHz]')
    ax.grid()
    plt.show()


    L = 1000*um
    h = 500*um
    period = 8*um
    eta = 0.5
    w = period*eta/2
    g = period/2 - w
    er = 11.9
    N = 502 - np.array([0, 17, 32, 47, 61, 74, 87, 100, 111, 122, 133, 144])
    #N = np.arange(250, 550)

    Cs = capacitance(w, g, h, er, N, L)
    p = np.polyfit(N, Cs/pF, 1)
    Nfit = np.arange(250, 550)
    Cfit = np.polyval(p, Nfit)
    Cfit2 = np.polyval([0.050, 0.788], Nfit)
    print (p)
    #exit()
    L = 10*nH
    fr = 1./(2*pi*np.sqrt(L*Cs))

    fig,ax = plt.subplots(figsize=(10,10))
    ax.plot(N, Cs/pF, 'ko', label='Analytical Calculation')
    ax.plot(Nfit, Cfit, 'r', label='Best Fit: ')
    ax.plot(Nfit, Cfit2, 'b', label='Best Fit to sonnet model')
    ax.set_xlabel('Number of Fingers')
    ax.set_ylabel('C [pF]')
    ax.legend(loc='upper left')
    ax.grid()
    plt.show()

    fig,ax = plt.subplots(figsize=(10,10))
    ax.plot(N, fr/MHz, 'ko')
    ax.set_xlabel('Number of Fingers')
    ax.set_ylabel('Resonance Frequency [MHz]')
    ax.grid()
    plt.show()


    frdes = np.array([266.1, 271.7, 276.9, 282.4, 287.6, 292.7, 298.0, 303.5,
        308.3, 313.3, 318.5, 323.4])

    frseen = np.sqrt(capacitance(2.0, 2.0, h, er, N, L
        )/capacitance(1.5, 2.5, h, er, N, L))*frdes
    print (frseen)





