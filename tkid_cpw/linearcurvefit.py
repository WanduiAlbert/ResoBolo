import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import astropy.units as u

# Using a best linear fit model to estimate the capacitance to ground for a given length of resonator tank
Y = np.array([0.408505, 0.453653, 0.476200, 0.498798, 0.543803])  # Capacitance in pF
C = Y * u.pF
L = np.array([1300, 1500, 1600, 1700, 1900]) * u.um
A = np.vstack([np.ones(5), [1300, 1500, 1600, 1700, 1900]])  # Lengths are in um
A = A.T

#Covariance matrix is simply the identity

# Performing a linear best  fit Y = A X.

X = np.linalg.solve(A.T.dot(A), A.T.dot(Y)) # X = [b, m]

m = X[1] * u.pF/u.um

C0 = X[0] * u.pF

variances = np.diag(np.linalg.inv(A.T.dot(A)))

# print (variances)

sigm = variances[1]**0.5 * u.pF/u.um
sigb = variances[0]**0.5 * u.pF

print ("The m equals {0} +/- {1} ".format(m, sigm ))

print ("The C0 equals {0} +/- {1} ".format(C0, sigb))

residues = C - (m*L + C0)

# Next part would be to ascertain that this does scale with the number of squares which should lead us to a way of computing the 
# parasitic capacitance for any section of the capacitor tank

.fig, ax = plt.subplots(figsize=(10,10))
ax.plot(L, Y, 'bx', markersize=15)
ax.plot(L, m *L  + C0, 'r')
ax.set_xlabel(r'Length [$\mu $m]')
ax.set_ylabel(r"Capacitance to Ground [pF]")
ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.4f}"))
ax.grid(which='both')
ax.axis('tight')


fig, ax = plt.subplots(figsize=(10, 10))
ax.stem(L, residues)
ax.set_xlabel(r'Length [$\mu $m]')
ax.set_ylabel(r"Residues[pF]")
ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.4e}"))
ax.grid(which='both')
ax.axis('tight')

plt.show()