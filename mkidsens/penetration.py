
import pint
import numpy as np

u = pint.UnitRegistry()
nm = u.nanometer

d = 400. * nm
lam = 110. * nm
mu0 = 4.*np.pi *1e-7 * u.newtons /u.amps**2

Lk = (mu0 * lam / np.tanh(d/lam)).to('picohenry')

print Lk
