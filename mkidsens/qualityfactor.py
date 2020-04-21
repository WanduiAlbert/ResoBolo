import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h,k,c
from scipy.special import kn

pi = np.pi
K0 = lambda x: kn(0, x)

Tc = 1.2
Delta = 1.763*k*Tc

alphak = 0.4

T = np.r_[90:500:1000j]*1e-3
f = 300e6

eta = h*f/(2*k*T)
Qi = pi*np.exp(Delta/(k*T))/(4*alphak*np.sinh(eta)*K0(eta))

plt.figure(figsize=(10,10))
plt.plot(T/1e-3, Qi, 'k', lw=3)
plt.axvline(380, color='r', ls='dashed')
plt.grid()
plt.yscale('log')
plt.xlabel('Temperature [mK]')
plt.ylabel('Qi')
plt.ylim(top=1e7)
plt.xlim(left=100, right=500)
plt.savefig('qualityfactor_vs_temp.pdf')
plt.show()
