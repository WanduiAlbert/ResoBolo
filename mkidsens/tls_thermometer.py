#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import digamma, polygamma
from scipy.constants import h,k
from math import pi
from julia.api import Julia
jl = Julia(compiled_modules=False)
jl.using("SpecialFunctions")
from julia import SpecialFunctions


MHz = 1e6
mK = 1e-3
ppm = 1e-6

F = 0.5
delta0 = 2e-3
PcdBm = -95
PdBm = np.r_[-110:-80:10j]
Pc = 1e-3*10**(PcdBm/10)
P = 1e-3*10**(PdBm/10)
f0 = 684.4284*MHz

df = 0.1#Hz
dx = df/f0
dQ = 1e-8

print (h*f0/k/mK)

T = np.r_[250:800:1000j]*mK
labels = ["%.1f dBm"%(_) for _ in PdBm]

delta = F*delta0/np.sqrt(1 + P/Pc)
dd = F*delta0/pi*np.ones_like(delta)
ny = delta.size
y = h*f0/(k*T)
dQi = delta*np.tanh(y[:, np.newaxis]/2)
x = dd*(np.real(digamma(0.5 + y/(2j*pi))) - np.log(y))[:, np.newaxis]
#f = f0*(1+x)
#x = (f - f[-1])/f[-1]

plt.figure(figsize=(10,10))
plt.plot(T/mK, x[:, 0]/ppm)
plt.xlabel('T [mK]')
plt.ylabel('x_tls [ppm]')
plt.grid()
plt.savefig('tls_thermometer_xppm.png')
plt.show()

plt.figure(figsize=(10,10))
for iy in range(ny):
    plt.plot(T/mK, dQi[:, iy], label=labels[iy])
plt.xlabel('T [mK]')
plt.ylabel('dQi_tls')
plt.legend(loc='upper right')
plt.grid()
plt.savefig('tls_thermometer_dQi.png')
plt.show()

dQidT = delta*((1 - np.tanh(y/2)**2)*(-y/T))[:, np.newaxis]
dxdT = np.zeros_like(dQidT)
nx, ny = dxdT.shape

for ix in range(nx):
    dxdT[ix, :] = np.real(SpecialFunctions.trigamma(0.5 + y[ix]/(2j*pi)))
dxdT -= 1./y[:, np.newaxis]
dxdT *= -dd
dxdT *= (y/T)[:, np.newaxis]

dT_fromx = np.abs(dx/dxdT)
dT_fromdQ = np.abs(dQ/dQidT)

plt.figure(figsize=(10,10))
plt.plot(T/mK, dxdT[:, 0]/ppm*mK)
plt.xlabel('T [mK]')
plt.ylabel('dx/dT [ppm/mK]')
plt.grid()
plt.savefig('tls_thermometer_xsensitivity.png')
plt.show()

plt.figure(figsize=(10,10))
for iy in range(ny):
    plt.plot(T/mK, dQidT[:, iy]/ppm*mK, label=labels[iy])
plt.xlabel('T [mK]')
plt.ylabel('dQi/dT [ppm/mK]')
plt.legend(loc='lower right')
plt.grid()
plt.savefig('tls_thermometer_dQisensitivity.png')
plt.show()

plt.figure(figsize=(10,10))
plt.plot(T/mK, dT_fromx[:, 0]/mK)
plt.xlabel('T [mK]')
plt.ylabel('dT [mK]')
plt.grid()
plt.savefig('tls_thermometer_temp_resolution_fromx.png')
plt.show()

plt.figure(figsize=(10,10))
for iy in range(ny):
    plt.plot(T/mK, dT_fromdQ[:, iy]/mK, label=labels[iy])
plt.xlabel('T [mK]')
plt.ylabel('dT [mK]')
plt.legend(loc='lower right')
plt.grid()
plt.savefig('tls_thermometer_temp_resolution_fromx.png')
plt.show()


