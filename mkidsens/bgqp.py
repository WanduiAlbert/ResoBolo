import pint
import numpy as np
from scipy import special
import matplotlib.pyplot as pl
pi = np.pi

'''
Model a resonator that has its frequency shift with temperature modified by an excess source of quasiparticles
A simple fixed background of quasiparticles with S2 won't work because that just shifts the frequency up and down.
Instead we try a source generating QPs at a constant rate and follow Jonas's review paper on how that creates frequency shift and Qi
'''

u = pint.UnitRegistry()

h = u.planck_constant
k = u.boltzmann_constant

Tc = 1.32 * u.kelvin
N0 = 1.7e10 / (u.electron_volt * u.micrometer**3)
taumax = 0.0005 * u.second
nstar = 100 / (u.micrometer**3)
fro = 300*u.MHz

delta0 = (1.76 * k * Tc).to('joule')

alphak = 0.4

To = np.linspace(0.08,0.35,100)*u.kelvin

nth = (2 * N0 * np.sqrt(2. * pi * k * To * delta0) * np.exp(-delta0 / (k * To))).to(1./u.micrometer**3)
delta0 = (1.76 * k * Tc).to('joule')
hf2kt = (h * fro / (2 * k * To)).to('')

y = np.sqrt(2*delta0 / (pi * k * To))
S1 = ((2./pi) * np.sqrt(2 * delta0 / (pi * k * To)) * np.sinh(hf2kt) * special.k0(hf2kt)).to('')
S2 = (1.0 +  np.sqrt(2 * delta0 / (pi * k * To)) * np.exp(-hf2kt) * special.i0(hf2kt)).to('')
print (1./hf2kt)
print (S1)
print (S2)
xth = (alphak * S2 * nth / (4 * N0 * delta0)).to('').m
Qith = (2 * N0 * delta0 / (alphak * S1 * nth)).to('').m
pl.figure(1)
pl.subplot(211)
pl.plot(To,xth*1e6,label='thermal qp')
pl.subplot(212)
pl.plot(To,1./Qith,label='thermal qp')

pl.figure(2)
pl.plot(To,nth,label='thermal qp')
nfloors = np.arange(1, 11, 1)*200
for nfloor in nfloors:
	#nfloor = 400/u.micrometer**3
	#gamma = (nfloor/u.micrometer**3)**2/(2*nstar*taumax)

	nqp = np.sqrt((nth + nstar)**2 + (nfloor/u.micrometer**3)**2)-nstar


	xqp = (alphak * S2 * nqp / (4 * N0 * delta0)).to('').m
	Qiqp = (2 * N0 * delta0 / (alphak * S1 * nqp)).to('').m

	pl.figure(1)
	pl.subplot(211)
	pl.plot(To,xqp*1e6,label="nfloor = %1.1f"%nfloor)
	pl.subplot(212)
	pl.plot(To,1./Qiqp,label="nfloor = %1.1f"%nfloor)

	pl.figure(2)
	pl.plot(To,nqp,label="nfloor = %1.1f"%nfloor)
	pl.legend(loc='upper left')

pl.figure(1)
pl.subplot(211)
pl.legend(loc='upper left')
pl.grid()
pl.ylabel('Frequency shift (ppm)')
pl.subplot(212)
pl.grid()
pl.ylabel('1/Q')
pl.xlabel('T (K)')
pl.gca().set_yscale('log')
pl.savefig('fqex.png')

pl.figure(2)
pl.legend(loc='upper left')
pl.grid()
pl.xlabel('T (K)')
pl.ylabel('Quasiparticle density (1/um^3)')
pl.savefig('nex.png')
pl.show()



