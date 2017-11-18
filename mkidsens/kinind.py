
import numpy as np
import pint
from math import pi

u = pint.UnitRegistry()
k = u.boltzmann_constant
hbar = u.hbar
mu_0 = u.mu_0

Tc = 1.32 * u.K
delta = 1.76 * k * Tc
rhon = 1.15 * u.microohm * u.cm
t = 0.05*u.micrometer
Rs = (rhon/t).to('ohm')
w = 1.0*u.micrometer
Rl = Rs / w
nsquare = 16200.

penetration_depth = np.sqrt(hbar * rhon / (pi * delta * mu_0)).to('nm')

Ls = (hbar * Rs / (pi * delta)).to('pH')
Ll = Ls / w

print "Tc: ",Tc
print "w: ",w
print "t: ",t
print "Sheet resistance of film:",Rs
print "penetration depth of film:",penetration_depth
print "Linear resistance of strip: ",Rl
print "Sheet inductance of film:",Ls
print "Linear inductance of film:",Ll

print "Total resistance: ",nsquare * Rs
print "Total kinetic inductance: ",nsquare * Ls
