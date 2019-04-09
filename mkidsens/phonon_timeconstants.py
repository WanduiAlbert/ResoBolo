#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.constants import k

mJ = 1e-3
um = 1e-6
pW = 1e-12
pJ = 1e-12
K = 1
mol = 1
cm = 1e-2
g = 1

# Material properties
rho = 2.7 * g/cm**3
Ar = 26.98 * g/mol
gamma = 1.35 * mJ/mol/K**2
Tc = 1.329 * K
delta = 1.76 * k * Tc
V_sc = 16200*0.05*um**3
Sigma = 0.2e9 # W/m^3/K^5
n = 2.104
K0 = 151.193 * pW/K**(n+1)

T = 0.320 * K
C_is = 0.25*(T/(0.350*K))*pJ/K
G_is = K0*(n+1)*T**n
C_qp = 7.1 * gamma*(rho*V_sc/Ar)*Tc*np.exp(-delta/(k*T))
G_eph = 5*Sigma*V_sc*T**4

#print (G_is/pW)
#print (C_is/pJ)
#print (G_eph/pW)
#print (C_qp/pJ)

G_eff = G_eph*G_is/(G_eph + G_is)

t1 = C_is/G_is
t2 = C_is*G_eff/(G_is*G_eph)
t3 = C_qp/G_eff


f = np.r_[1:1e5:1000j]
w = 2*pi*f

Y = G_eff * ((1 + 1j*w*t1)/(1 + 1j*w*t2) + 1j*w*t3)
Z = 1./Y

plt.semilogx(f, np.real(Z))
plt.xlabel('Frequency [Hz]')
plt.ylabel('| Thermal Z |')
plt.grid()
plt.show()


# Second try
C_t = C_is + C_qp
T1 = C_qp/G_eph
T2 = T1 * C_is/C_t
T3 = C_t/G_is

f1 = 1./(2*pi*T1)
f2 = 1./(2*pi*T2)
f3 = 1./(2*pi*T3)


#t2 = 0
#t3 = 0
print (f1)
print (f2)
print (f3)
T2 = 0
T1 = 0

Y = G_is/(1 + 1./(1j*w*T3)*((1+1j*w*T1)/(1 + 1j*w*T2)))
Z = 1./Y
plt.loglog(f, np.abs(Z))
plt.xlabel('Frequency [Hz]')
plt.ylabel('| Thermal Z |')
plt.grid()
plt.show()


