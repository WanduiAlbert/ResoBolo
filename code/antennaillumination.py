#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.constants import c

mm = 1e-3
um = 1e-6
GHz = 1e9


er = 11.9
nu0 = 150*GHz
lambda0 = c/nu0

a = 302*um
b = 302*um

N = 12

# Solve for theta at FWHM
x = 0.1
for i in range(30):
    x = np.arcsin(2**0.5/N*np.sin(N*x))

print (x)
theta_fwhm = 2**0.5*np.arcsin(lambda0*x/(pi*a*er**0.5))*180/pi
print ("FWHM: ", theta_fwhm)

x = np.r_[-5:5:1000j]*mm
y = np.r_[-5:5:1000j]*mm


theta = np.r_[0:pi:1000j]
phi = np.r_[-pi:pi:1000j]

Theta, Phi = np.meshgrid(theta, phi)

X = Theta*np.cos(Phi)*180/pi
Y = Theta*np.sin(Phi)*180/pi

AF = np.sin(N*pi*a*er**0.5*np.sin(Theta)*np.cos(Phi)/lambda0)
AF /= np.sin(pi*a*er**0.5*np.sin(Theta)*np.cos(Phi)/lambda0)
AF *= np.sin(N*pi*b*er**0.5*np.sin(Theta)*np.sin(Phi)/lambda0)
AF /= np.sin(pi*b*er**0.5*np.sin(Theta)*np.sin(Phi)/lambda0)
AF[np.isnan(AF)] = N**2
AF = np.abs(AF)/N**2
#AF[AF < 0.5] = np.nan

#FF = np.sin(Theta)*np.sin(pi*a*er**0.5*np.cos(Theta)/lambda0)
#FF /= (pi*a*er**0.5*np.cos(Theta)/lambda0)
#FF = np.abs(FF)**2
#AF = FF
AF = 10*np.log10(AF)

plt.figure(figsize=(10,10))
contours = plt.contour(X, Y, AF, levels=10, colors='black')
#plt.clabel(contours, inline=True, fontsize=10)
plt.pcolormesh(X, Y, AF, alpha=0.5)
plt.plot(np.cos(phi)*theta_fwhm, np.sin(phi)*theta_fwhm, 'r--', lw=4)
plt.xlabel('Theta cos(phi) [Deg]')
plt.ylabel('Theta sin(phi) [Deg]')
plt.title('Array Factor %dX%d array'%(N,N))
plt.axis('square')
plt.xlim(left=-30, right=30)
plt.ylim(bottom=-30, top=30)
plt.grid()
plt.colorbar()
plt.savefig('arrayfactor_%dX%d.png'%(N,N))
plt.show()
