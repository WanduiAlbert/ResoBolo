
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import reso_fit
from scipy import optimize

def do_fit(fs,s21):
	f0 = fs[np.argmin(np.abs(s21))]
	A = 0.5
	phi = 0.
	D = 0.
	Qr = 30000
	Qe_re = 60000
	Qe_im = 0
	p0 = (f0,A,phi,D,Qr,Qe_re,Qe_im)

	def resid_func(p):
		resid = reso_fit.real_of_complex(s21 - reso_fit.model(fs,*p))
		return resid

	'''
	pl.plot(fs,np.abs(s21))
	ym = reso_fit.model(fs,*p0)
	pl.plot(fs,np.abs(ym))
	pl.show()
	exit()
	'''

	p,ier = optimize.leastsq(resid_func,p0)
	f0,A,phi,D,Qr,Qe_re,Qe_im = p
	Qi = 1.0 / (1.0/Qr - 1.0/Qe_re)

	print Qr,Qi,Qe_re,Qe_im

def main():

	mainC = 40e-12
	mainL = 6e-9
	mainQ = 20000
	mainR = mainQ * np.sqrt(mainL/mainC)
	mainf0 = 1.0 / (2.0 * np.pi * np.sqrt(mainC*mainL))

	fs = np.linspace(mainf0 * 0.996,mainf0 * 1.004,10000)
	ws = 2.0 * np.pi * fs

	I = 1.0j
	Zmc = 1.0 / (I * ws * mainC)
	Zml = I * ws * mainL
	Zm = 1.0 / (1.0/Zmc + 1.0/Zml + 1.0/mainR)
	
	coupC = 0.1e-12
	Zcc = 1.0 / (I * ws * coupC)

	parasiticC = 0.3e-12
	Zp = 1.0 / (I * ws * parasiticC)

	Ztank1 = Zcc + Zm + Zp
	Ztank2 = Zcc + 1.0 / (1.0/Zp + 1.0 / (Zm + Zp))

	Z0 = 50.0

	Ztank1 = 1.0 / (1.0/Ztank1 + 1.0/Z0)
	t1 = Ztank1 / (Ztank1 + Z0)

	Ztank2 = 1.0 / (1.0/Ztank2 + 1.0/Z0)
	t2 = Ztank2 / (Ztank2 + Z0)

	fs *= 1e-6

	do_fit(fs,t1)
	do_fit(fs,t2)

	plt.plot(fs,np.abs(t1),'b')
	plt.plot(fs,np.abs(t2),'r')
	plt.grid()
	plt.show()

if __name__=='__main__':
	main()
