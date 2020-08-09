#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import scipy.linalg as linalg

N = 10

nH = 1e-9
pF = 1e-12
MHz = 1e6

Nreso = 5
L = 10*nH
f0 = 300*MHz
df = 5*MHz
frs0 = f0 + np.arange(Nreso)*df
Z0 = 50


Cs = 1./((2*pi*frs0)**2*L)
indices = np.argsort(Cs/pF)

Qc = 6000
Ccs = np.sqrt(2*Cs/(2*pi*frs0*Z0*Qc))

#print (Cs[indices]/pF)

Nsteps = 20
C12s = np.linspace(0, 0.3, Nsteps)*pF

frs = np.zeros((Nsteps, Nreso))
Qcs = np.zeros((Nsteps, Nreso))


for istep, C12 in enumerate(C12s):
    cmat = np.diag(Cs + 2*C12)
    cmat[0,0] -= C12
    cmat[-1,-1] -= C12
    cmat += np.diag(-np.array([C12]*(Nreso-1)), k=1)
    cmat += np.diag(-np.array([C12]*(Nreso-1)), k=-1)
    cmat /= pF
    w, v = linalg.eigh(cmat)
    v = v.T
    frs[istep] = 1./(2*pi*np.sqrt(L*w*pF))
    if istep == 0:
        print (frs0/MHz)
        print (frs[istep][indices]/MHz)
    #if istep==3:
    #    print (frs[istep]/MHz)
    #    print (Cs[indices]/pF)
    #    print (Ccs[indices]/pF)
    for ireso in range(Nreso):
        #print ("\n%d zeroth mode: "%i, w[0], v[0])
        #print ("%d first mode: "%i, w[1], v[1])

        num = np.abs(v[ireso])**2*Cs[indices]
        #num = np.diag(np.outer(v[ireso], np.conj(v[ireso])))*Cs[indices]
        outer = np.outer(v[ireso]*Ccs[indices], np.conj(v[ireso])*Ccs[indices])
        Qcs[istep, ireso] = (2*np.sum(num)/(Z0*2*pi*frs[istep][indices][ireso]*np.sum(outer)))


plt.figure(figsize=(10,10))
for ireso in range(Nreso):
    plt.plot(C12s/pF, frs[:, ireso]/MHz, marker='o', ls='None', label='%d'%ireso)
plt.grid()
plt.legend(loc='upper right')
plt.xlabel('C12 [pF]')
plt.ylabel('fr [MHz]')
plt.show()

plt.figure(figsize=(10,10))
for ireso in range(Nreso):
    plt.plot(C12s/pF, Qcs[:, ireso], marker='o', ls='None', label='%d'%ireso)
plt.grid()
plt.legend(loc='upper right')
plt.xlabel('C12 [pF]')
plt.ylabel('Qc')
plt.show()

plt.figure(figsize=(10,10))
for istep in range(Nsteps):
    plt.plot(frs[istep, :]/MHz, Qcs[istep, :], marker='o', ls='None',
            label='C_coup=%1.3fpF'%(C12s[istep]/pF))
plt.grid()
#plt.legend(loc='upper right')
plt.xlabel('fr [MHz]')
plt.ylabel('Qc')
plt.show()

