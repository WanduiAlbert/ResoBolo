
import pint
import numpy as np
import matplotlib.pyplot as pl
from scipy import special


u = pint.UnitRegistry()
pi = np.pi

f = 300*u.megahertz
T = np.linspace(0.1,0.3)*u.kelvin
h = u.planck_constant
k = u.boltzmann_constant

t1 = np.real(special.digamma(0.5 - h*f/(2j*pi*k*T)))
#t1 = np.log(h*f/(2*pi*k*T))
t2 = -np.log(T.m)

t1 = t1 - t1[0]
t2 = t2 - t2[0]

pl.plot(T,t1,label='digamma')
pl.plot(T,t2,label='log')
pl.grid()
pl.legend()
pl.xlabel('Temperature (K)')
pl.ylabel('Frequency shift')
pl.show()


