
import numpy as np
from scipy import optimize
from math import pi
import matplotlib.pyplot as pl

def real_of_complex(z):
	''' flatten n-dim complex vector to 2n-dim real vector for fitting '''
	r = np.hstack((z.real,z.imag))
	return r
def complex_of_real(r):
	assert len(r.shape) == 1
	nt = r.size
	assert nt % 2 == 0
	no = nt/2
	z = r[:no] + 1j*r[no:]
	return z

def model(f,f0,A,phi,D,Qr,Qe_re,Qe_im):
	f0 = f0 * 1e6
	f = f * 1e6
	cable_phase = np.exp(-2.j*pi*(1e-6*D*f+phi))
	Qe = Qe_re + 1.j*Qe_im
	x = (f - f0)/f0
	s21 = A*cable_phase*(1. - (Qr/Qe)/(1. + 2.j*Qr*x))
	return s21

def main():
	fn = '/home/lebicep_admin/netanal/20160715_171928.dat'
	data = np.loadtxt(fn,skiprows=5)
	freqs,re,im = data.T
	phase0 = np.mean(np.unwrap(np.arctan2(im,re)))
	z = re + 1j*im
	z *= np.exp(-1j*phase0)
	nt = len(freqs)

	f0 = np.mean(freqs)*1e-6
	A = 0.001
	phi = 0.
	D = 0.
	Qr = 10000
	Qe_re = 10000
	Qe_im = 0
	p0 = (f0,A,phi,D,Qr,Qe_re,Qe_im)

	ydata = np.hstack((z.real,z.imag))
	print data.shape

	popt,pcov = optimize.curve_fit(model,freqs,ydata,p0=p0)
	f0,A,phi,D,Qr,Qe_re,Qe_im = popt

	Qi = 1.0/ (1./Qr - 1./Qe_re)
	print "f0: ",f0
	print "A: ",A
	print "phi: ",phi
	print "D: ",D
	print "Qr: ",Qr
	print "Qe_re: ",Qe_re
	print "Qe_im: ",Qe_im
	print "Qi: ",Qi

	ymodel = model(freqs,*popt)
	ymodel_re = ymodel[:nt]
	ymodel_im = ymodel[nt:]
	freqs *= 1e-6
	pl.plot(freqs,z.real,color='blue',linestyle='.',marker='.',label='re data')
	pl.plot(freqs,z.imag,color='green',linestyle='.',marker='.',label='im data')
	pl.plot(freqs,ymodel_re,color='blue',linestyle='-',label='re fit')
	pl.plot(freqs,ymodel_im,color='green',linestyle='-',label='im fit')
	pl.xlabel('Frequency (MHz)')
	pl.ylabel('S21')
	pl.grid()
	pl.legend(loc='upper left')
	pl.title('158.7MHz resonator\nQi=%d Qr=%d'%(Qi,Qr))
	pl.savefig('fit_158p7MHz.png')
	pl.show()

if __name__=='__main__':
	main()
