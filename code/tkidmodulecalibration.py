#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.special import kn, iv
from scipy.constants import h,k,e
K0 = lambda x: kn(0,x)
I0 = lambda x: iv(0,x)

um = 1e-6
N0 = 1.72e10/um**3/e

N = 16
df = 3
f0 = 300
f = f0 + np.arange(N)*df
k = 1000

fs = 2e3
dt = 1./fs

R0 = 0.155
sigr = 0.0019
R = R0 + np.random.randn(N)*sigr

Idc = 5.7e-6
deltaI = 0.05e-6

t = dt*np.arange(k)
ft = 10
I = Idc + deltaI*np.sin(2*pi*ft*t)

P = I**2*R[:, np.newaxis]

T0 = 250e-3
K = 120e-12 + 1e-12*np.random.randn(N)
n = 1.9 + 0.05*np.random.randn(N)
T = (P/K[:, np.newaxis] + T0**(n[:, np.newaxis] + 1))**(1./(n[:, np.newaxis] + 1))

alphak = 0.51 + 0.05*np.random.randn(N)
Tc = 1.388 + 0.005*np.random.randn(N)
Delta = 1.764*k*Tc

nqp = 2*N0*np.sqrt(2*pi*k*T)*np.exp(-Delta[:, np.newaxis]/k/T)
eta = h*f[:, np.newaxis]/(2*k*T)
S2 = 1 + np.sqrt(2*Delta[:, np.newaxis]/(pi*k*T))*np.exp(-eta)*I0(eta)
x = alphak[:, np.newaxis]*S2/(4*N0*Delta[:, np.newaxis])*nqp

Adiag = (x - np.mean(x, axis=1)[:, np.newaxis])/(P - np.mean(P, axis=1)[:, np.newaxis])
Adiag = np.diag(np.mean(Adiag, axis=1))
#Adiag = np.diag(np.mean(x, axis=1)/np.mean(P, axis=1))
# Add in a little bit of crosstalk
crosstalk = np.abs(0.005*np.random.randn(N, N))
np.fill_diagonal(crosstalk, 0)
x = np.dot(np.eye(N) + crosstalk, x)

x0 = np.mean(x, axis=1)[:, np.newaxis]
dx = x - x0
P0 = np.mean(P, axis=1)[:, np.newaxis]
dP = P - P0

# Now I can write my "true" estimate for the gain and hopefully see how well I
# can reconstruct it
A = np.dot(np.eye(N) + crosstalk, Adiag)
Utrue, Strue, Vtrueh = np.linalg.svd(A)
#print (Strue/1e6)

# I only have access to the total power summed along the rows
dPtot = np.sum(dP, axis=0)
# Now I can construct the J vector
J = np.ones(N)/N
covx = np.dot(dx.T, dx)/(N-1)
corrx = np.dot(dx, dx.T)/(N-1)

covA = corrx/(np.sum(dPtot**2)/N)
V, S, Vh = np.linalg.svd(covA)
#print(np.sqrt(S)/1e6)

#Aest = np.dot(U, np.dot(np.diag(np.sqrt(S)), Vh))*(N-1)
#
#print (A)
#print (Aest)
dPest = dPtot/N
Adiag_est = np.mean((dx/dPest), axis=1)

error = (Adiag_est - np.diag(A))/np.diag(A)
print (error*100)


