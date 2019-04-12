
import numpy as np
from scipy import optimize
from math import pi
import matplotlib.pyplot as pl
import pdb

def real_of_complex(z):
	''' flatten n-dim complex vector to 2n-dim real vector for fitting '''
	r = np.hstack((z.real,z.imag))
	return r
def complex_of_real(r):
	assert len(r.shape) == 1
	nt = r.size
	assert nt % 2 == 0
	no = nt//2
	z = r[:no] + 1j*r[no:]
	return z

def model_linear(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a):
	f0 = f0 * 1e6
	cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
	dQe = dQe_re + 1.j*dQe_im
	x = (f - f0)/f0
	s21 = A*cable_phase*(1. - dQe/(dQr + 2.j*x))
	return real_of_complex(s21)

def model_linear_wslope(f,f0,A,m,phi,D,dQr,dQe_re,dQe_im,a):
	f0 = f0 * 1e6
	cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
	dQe = dQe_re + 1.j*dQe_im
	x = (f - f0)/f0
	s21 = (A + m*(f-f0)*1e-6)*cable_phase*(1. - dQe/(dQr + 2.j*x))
	return real_of_complex(s21)

def get_roots(y0, a):
	coeffs = [4, -4*y0, 1, -(y0+a)]
	roots = np.roots(coeffs)
	return roots[np.isreal(roots)]

def get_y(y0, a, f):
	low_to_high = np.all(np.diff(f) > 0)
	if low_to_high:
		find_roots = lambda x: np.min(get_roots(x, a))
	else:
		find_roots = lambda x: np.max(get_roots(x, a))
	ys = list(map(find_roots, y0))
	return np.array(ys)

def make_model(dQe_re,dQe_im):
	def f(f,f0,A,phi,D,dQr,a):
		return full_model(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a)
	return f

def model(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a):
	f0 = f0 * 1e6
	cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
	dQe = dQe_re + 1.j*dQe_im
	x0 = (f - f0)/f0
	y0 = x0/dQr
	k2 = np.sqrt((y0**3/27. + y0/12. + a/8.)**2 - (y0**2/9. - 1/12.)**3, dtype=np.complex128)
	k1 = np.power(a/8. + y0/12. + k2 + y0**3/27., 1./3)
	eps = (-1. + 3**0.5 * 1j)/2.

	#np.seterr(all='raise')
	#print (np.sum(np.abs(k1) == 0.0))
	y1 = y0/3. + (y0**2/9.-1/12.)/k1 + k1
	y2 = y0/3. + (y0**2/9.-1/12.)/eps/k1 + eps*k1
	y3 = y0/3. + (y0**2/9.-1/12.)/eps**2/k1 + eps**2*k1

	y1[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
	y2[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
	y3[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.

	# Out of the three roots we need to pick the right branch of the bifurcation
	thresh = 1e-4
	low_to_high = np.all(np.diff(f) > 0)
	if low_to_high:
		y = y2.real
		mask = (np.abs(y2.imag) >= thresh)
		y[mask] = y1.real[mask]
	else:
		y = y1.real
		mask = (np.abs(y1.imag) >= thresh)
		y[mask] = y2.real[mask]

	#if np.sum(k1 == 0.0) > 0 :
	#	y[k1 == 0.0] = y0[k1 == 0.0]/3
	#y = get_y(y0, a, f)
	x = y*dQr
	s21 = A*cable_phase*(1. - (dQe)/(dQr + 2.j*x))
	return real_of_complex(s21)

def model_wslope(f,f0,A,m,phi,D,dQr,dQe_re,dQe_im,a):
	f0 = f0 * 1e6
	cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
	dQe = dQe_re + 1.j*dQe_im
	x0 = (f - f0)/f0
	y0 = x0/dQr
	k2 = np.sqrt((y0**3/27. + y0/12. + a/8.)**2 - (y0**2/9. - 1/12.)**3, dtype=np.complex128)
	k1 = np.power(a/8. + y0/12. + k2 + y0**3/27., 1./3)
	eps = (-1. + 3**0.5 * 1j)/2.

	#np.seterr(all='raise')
	#print (np.sum(np.abs(k1) == 0.0))
	y1 = y0/3. + (y0**2/9.-1/12.)/k1 + k1
	y2 = y0/3. + (y0**2/9.-1/12.)/eps/k1 + eps*k1
	y3 = y0/3. + (y0**2/9.-1/12.)/eps**2/k1 + eps**2*k1

	y1[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
	y2[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
	y3[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.

	# Out of the three roots we need to pick the right branch of the bifurcation
	thresh = 1e-4
	low_to_high = np.all(np.diff(f) > 0)
	if low_to_high:
		y = y2.real
		mask = (np.abs(y2.imag) >= thresh)
		y[mask] = y1.real[mask]
	else:
		y = y1.real
		mask = (np.abs(y1.imag) >= thresh)
		y[mask] = y2.real[mask]

	#if np.sum(k1 == 0.0) > 0 :
	#	y[k1 == 0.0] = y0[k1 == 0.0]/3
	#y = get_y(y0, a, f)
	x = y*dQr
	s21 = (A + m*(f-f0)*1e-6)*cable_phase*(1. - (dQe)/(dQr + 2.j*x))
	return real_of_complex(s21)

def chisq(theta,x,y,yerr):
	inv_sigma2 = 1/yerr**2
	return 0.5*np.sum((y - model(x, *theta))**2*inv_sigma2)

def model_Qi(f,f0,A,phi,D,Qr,Qi,Qe_im):
	f0 = f0 * 1e6
	cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
	Qe_re = 1./(1./Qr - 1./Qi)
	Qe = Qe_re + 1.j*Qe_im
	x = (f - f0)/f0
	s21 = A*cable_phase*(1. - (Qr/Qe)/(1. + 2.j*Qr*x))
	return real_of_complex(s21)


def do_fit(freq,re,im,plot=False,get_cov=False,verbose=False):
	nt = len(freq)
	#pdb.set_trace()
	mag = np.sqrt(re*re+im*im)
	phase = np.unwrap(np.arctan2(im,re))

	m,b = np.polyfit(freq,phase,1)
	D = m*1e6/(2*pi)
	D = 0.
	f0 = freq[np.argmin(mag)]*1e-6
	A = np.max(mag)
	m = 0
	phi = np.mean(phase)/(2*pi)
	#phi = b
	dQe_im = 0
	#p0 = (f0,A,phi,D,dQr,dQe_re,dQe_im)
	a = 0.1
	pBound = ([-np.inf]*9, [np.inf]*9)
	pBound[0][5] = 0
	pBound[0][6] = 0
	#pBound[1][5] = 1e8

	ydata = np.hstack((re,im))

	pdQr = 1./36000
	pdQe_re = 1./20000
	p0 = (f0,A,m,phi,D,pdQr,pdQe_re,dQe_im,a)
	while True:
		if plot:
			popt = p0
			break
		else:
			popt,pcov = optimize.curve_fit(model_wslope,freq,ydata,method='lm',
					p0=p0)#, bounds=pBound)
			f0,A,m,phi,D,dQr,dQe_re,dQe_im,a = popt
			if dQr < 0 or dQe_re < 0 or dQr < dQe_re:
				pdQr *= 2
				pdQe_re *= 1.5
				p0 = (f0,A,m,phi,D,pdQr,pdQe_re,dQe_im,a)
			else:
				break

	dQe = dQe_re + 1j*dQe_im
	Qe = 1/dQe
	Qr = 1/dQr
	Qi = 1.0/ (dQr - dQe_re)
	if verbose:
		print("f0: ",f0)
		print("A: ",A)
		print("m: ",m)
		print("phi: ",phi)
		print("D: ",D)
		print("Qr: ",Qr)
		print("Qe_re: ",Qe.real)
		print("Qe_im: ",Qe.imag)
		print("Qi: ",Qi)
		print("a: ", a)

	ymodel = model_wslope(freq,*popt)
	ymodel_re = ymodel[:nt]
	ymodel_im = ymodel[nt:]
	freq = freq * 1e-6

	if plot:
		'''
		pl.plot(freq,re,color='blue',linestyle='.',marker='.',label='re data')
		pl.plot(freq,im,color='green',linestyle='.',marker='.',label='im data')
		pl.plot(freq,ymodel_re,color='blue',linestyle='-',label='re fit')
		pl.plot(freq,ymodel_im,color='green',linestyle='-',label='im fit')
		'''
		pl.plot(freq,re*re+im*im,color='blue',linestyle='.',marker='.',label='abs data')
		pl.ylim(ymin=0)
		ax2 = pl.gca().twinx()
		ax2.plot(freq,np.arctan2(im,re),color='green',linestyle='.',marker='.',label='phase data')
		pl.plot(freq,ymodel_re*ymodel_re+ymodel_im*ymodel_im,color='blue',linestyle='-',label='abs fit')
		ax2.plot(freq,np.arctan2(ymodel_im,ymodel_re),color='green',linestyle='-',label='phase fit')
		pl.xlabel('Frequency (MHz)')
		pl.ylabel('S21')
		pl.grid()
		pl.legend(loc='upper left')
		pl.title('Qi=%d Qr=%d Qere=%d'%(Qi,Qr,Qe_re))
		pl.show()
		#exit()

	if get_cov:
		return f0, A, m, phi, D, Qi,Qr, Qe.real, Qe.imag, a, ymodel_re,ymodel_im, pcov
	return f0,Qi,Qr, Qe.real, Qe.imag, a,ymodel_re,ymodel_im

def main():
	import sys
	fn=sys.argv[1]
	data = np.loadtxt(fn,skiprows=5).T
	freq = data[0]
	res = data[1::2]
	ims = data[2::2]
	nm = res.shape[0]

	print(freq.shape)
	print(res.shape)
	print(ims.shape)

	f0s = []
	Qis = []
	Qrs = []
	powers = np.linspace(-10,16,nm) - 90
	for i in range(nm):
			f0,Qi,Qr,Qe_re = do_fit(freq,res[i],ims[i])
			f0s.append(f0)
			Qis.append(Qi)
			Qrs.append(Qr)
	f0s = np.array(f0s)
	Qis = np.array(Qis)
	Qrs = np.array(Qrs)
	#pl.plot(powers,f0s*1e6)
	pl.plot(powers,Qrs,label='Qr')
	pl.plot(powers,Qis,label='Qi')
	pl.xlabel('Power (dBm)')
	pl.ylabel('Q')
	pl.grid()
	pl.title('Q vs. power at %d MHz'%(f0*1e6))
	pl.legend()
	pl.ylim(0,100e3)
	pl.savefig('qpower_%dMHz.png'%(f0*1e6))
	pl.show()

if __name__=='__main__':
	main()
