import numpy as np
import matplotlib.pyplot as pl
import reso_fit
from scipy.interpolate import RectBivariateSpline
pi = np.pi
from scipy import optimize
from scipy.constants import c, epsilon_0

# transmission line
# Nb/Nb, 0.3um ILD, 20um MS, 80um CPW 20um wide slot

Z0 = 50
xline_eps = 13.2 - 0.0081j
xline_n = np.sqrt(xline_eps)
xline_Z0 = 45.
xline_ltot = 0.502

def abcd_shuntimp(Z):
	zero = np.zeros_like(Z)
	one = np.ones_like(Z)
	abcd = np.array([[one,zero],[1./Z,one]])
	return abcd

def abcd_xline(fs,Z0,n,l):
	gamma = 1.0j*2*pi*fs*n/c
	cosh = np.cosh(gamma*l)
	sinh = np.sinh(gamma*l)
	abcd = np.array([[cosh,Z0*sinh],[(1./Z0)*sinh,cosh]])
	return abcd

def s_of_abcd(abcd,Z0):
	[[a,b],[c,d]] = abcd
	b = b/Z0
	c = c*Z0
	denom = a + b + c + d
	s11 = (a + b - c - d) / denom
	s12 = 2*(a*d - b*c) / denom
	s21 = 2. / denom
	s22 = (-a + b - c + d) / denom
	s = np.array([[s11,s12],[s21,s22]])
	return s

def periodic_model(t, A, B, nu, phi):
	return A + B*np.cos(2*pi*nu*t + phi)

def fabryperot_model(t, A, B, C, D, nu, phi):
	#return A + B/(1 + C - 2*D*np.cos(2*pi*nu*t + phi))#*1/(1 + D*np.sin(2*pi*nu*t + phi)**2)
	return A + B/( C+np.cos(2*pi*nu*t + phi)**2 + D*np.sin(2*pi*nu*t + phi))

def qc_prediction(fs,C, Cc, Zc, n, l):
	w0 = 2*pi*fs
	gamma = 1.0j*2*pi*fs*n/c
	prefactor = 2*(C)/(Cc**2*Z0*w0)#*(1 + (Cc*Z0*w0/2)**2)
	cosh_tot = np.real(np.cosh(gamma*xline_ltot))
	sinh_tot = np.imag(np.sinh(gamma*xline_ltot))
	cosh_in = np.real(np.cosh(2*gamma*l))
	sinh_in = np.imag(np.sinh(2*gamma*l))
	cosh_out = np.real(np.cosh(2*gamma*(xline_ltot-l)))
	sinh_out = np.imag(np.sinh(2*gamma*(xline_ltot-l)))
	lambd = Zc/Z0
	num = -(4*lambd**2*cosh_tot**2 + (1+lambd**2)**2*sinh_tot**2)
	denom = lambd**2*(-2*(1+lambd**2) + (-1+lambd**2)*cosh_in +\
			(1-lambd**2)*cosh_out)
	return prefactor * num/denom

ncf = 11
nfrac = 101
fracs = np.linspace(0,1,nfrac)
cfs = np.linspace(290e6,500e6,ncf)
#Qes = np.zeros(nfrac,dtype=np.complex)

cmap = pl.get_cmap('jet')
all_Qes = []
all_f0s = []
for k in range(ncf):
	print (k)
	f0 = cfs[k]
	L0 = 10e-9
	Qi = 100000
	#Cc = 0.1e-12
	Qc = 40000.
	C0 = 1/((2*pi*f0)**2*L0)
	Cc = np.sqrt(C0/(pi*f0*Qc*50.))
	C0 -= Cc
	#Qc = C0/(pi*f0*Cc**2*50.)
	R0 = Qi * np.sqrt(L0/C0)
	#c = 3e8
	Qr = 1/(1/Qi + 1/Qc)

	df = 10*f0/Qr
	fs = np.linspace(f0-df,f0+df,1001)
	ws = 2*pi*fs
	YC0 = 1.0j * ws * C0
	YL0 = 1.0/(1.0j*ws*L0)
	YR0 = 1.0/R0
	ZCc = 1.0/(1.0j*ws*Cc)


	Yreso = YC0 + YL0 + YR0
	Zreso = 1.0 / Yreso
	Zresocc = Zreso + ZCc
	Qes = [0]*nfrac
	f0s = [0]*nfrac
	for j in range(nfrac):
		frac = fracs[j]
		abcd_xline1 = abcd_xline(fs,xline_Z0,xline_n,xline_ltot*frac)
		abcd_shunt = abcd_shuntimp(Zresocc)
		abcd_xline2 = abcd_xline(fs,xline_Z0,xline_n,xline_ltot*(1-frac))
		assert np.all(np.isfinite(abcd_xline1))

		abcd = np.zeros_like(abcd_shunt)
		for i in range(abcd_xline1.shape[-1]):
			tmp = np.dot(abcd_shunt[:,:,i],abcd_xline2[:,:,i])
			abcd[:,:,i] = np.dot(abcd_xline1[:,:,i],tmp)

		s = s_of_abcd(abcd,50.)

		def dbmag(x):
			return 10*np.log10(np.abs(x))

		s21 = s[1,0]
		baseline = 0.5*(s21[0] + s21[-1])
		s21 /= baseline
		f0f, Qif, Qrf, Qe_re, Qe_im, a, yf_re, yf_im = reso_fit.do_fit(fs,
				s[1,0].real, s[1,0].imag)
		Qe = Qe_re + 1j*Qe_im
		dQe = 1./Qe
		Qc = 1./np.real(dQe)
		yfit = yf_re + 1j*yf_im
		Qes[j] = Qc
		f0s[j] = f0f*1e6 - f0

	# We now have Qe as a function of the position along the feedline nfrac
	Qes = np.array(Qes)
	xfine = np.linspace(0,1,1<<16)
	Qcs_pred = qc_prediction(f0,C0, Cc, xline_Z0, xline_n, xline_ltot*xfine)
	Qesinterp = np.interp(xfine, fracs, Qes, period=1)
	pl.figure()
	pl.title("Freq %1.3fMHz"%(f0/1e6))
	pl.plot(xfine, Qcs_pred,'r-')
	pl.plot(xfine, Qesinterp,'k--')
	pl.plot(fracs, Qes, 'bh', ls='None')
	pl.grid()
	#pl.ylim((0, 50000))
	pl.show()
	continue
	#exit()
	#pl.figure()
	#pl.plot(xfine, np.abs(Qesinterp),'r-')
	#pl.plot(fracs, np.abs(Qes), 'bh', ls='None')
	#pl.grid()
	#pl.show()
	#exit()
	#Qefft = np.fft.fft(Qesinterp - np.mean(Qesinterp))
	#fracfft = np.fft.fftfreq(Qesinterp.size, xfine[1] - xfine[0])
	#Qefft = np.fft.fft(Qesinterp)
	#fracfft = np.fft.fftshift(np.fft.fftfreq(Qesinterp.size, xfine[1] - xfine[0]))
	#pl.figure()
	#pl.plot(fracfft, np.abs(Qefft))
	#pl.grid()
	#pl.yscale('log')
	#pl.xlim(left=0, right=50)
	#pl.show()
	#exit()
	#A = np.max(np.abs(Qes))
	#B = A - np.min(np.abs(Qes))
	#C = B/70
	#p0 = [A, -B, C, 0.15*C, 3.8, -pi/9]
	#Qerefine = fabryperot_model(xfine, *p0)
	#pl.figure()
	#pl.plot(xfine, Qerefine,'r-')
	#pl.plot(fracs, np.abs(Qes), 'bh', ls='None')
	#try:
	#	popt, pcov = optimize.curve_fit(fabryperot_model,fracs, np.abs(Qes), p0=p0)
	#	#print (popt)
	#	Qeabsfine = fabryperot_model(xfine, *popt)
	#	pl.plot(xfine, Qeabsfine,'g--')
	#except RuntimeError:
	#	pass
	#	#continue
	#pl.grid()
	#pl.show()
	#exit()
	#p0 = [np.min(Qes.imag), np.mean(Qes.imag), 1, 0]
	#popt, pcov = optimize.curve_fit(periodic_model,fracs, Qes.imag, p0=p0)
	#print (popt)



	color = cmap(k/(ncf-1.))
	pl.subplot(211)
	pl.plot(fracs,np.abs(Qes),label='%d MHz'%(f0*1e-6),color=color)
	pl.subplot(212)
	pl.plot(fracs,np.angle(Qes)*180./pi,label='%d MHz'%(f0*1e-6),color=color)
	all_Qes.append(Qes)
	all_f0s.append(f0s)

pl.subplot(211)
pl.title('Qe variation with resonator position l=502mm eps=13 Z=30')
pl.grid()
pl.legend()
pl.ylabel('|Qe|')
pl.subplot(212)
lgd = pl.legend(ncol=3, loc='upper left', bbox_to_anchor=(1.0,1.0))
pl.grid()
pl.xlabel('Resonator position along line')
pl.ylabel('angle(Qe) (deg)')
pl.savefig('mismatch02.png', bbox_extra_artists=(lgd,), bbox_inches='tight')


# Also want to make a scatter plot to overlay the design frequencies of the
# resonators vs their position on the feedline
#design_freqs = np.array([306, 318, 321, 324, 327, 330, 333, 336, 339, 351, 354,
#	369, 372, 381, 384, 387, 390, 393, 396, 399, 402, 405, 408, 417, 420, 435,
#	438, 450, 453, 456, 459, 462, 465, 468, 471, 483])
design_freqs = np.array([420, 417, 408, 381, 372, 369, 354, 387, 384, 405, 402,
	435, 438, 399, 396, 393, 390, 351, 336, 333, 330, 459, 456, 453, 450, 465,
	462, 327, 324, 339, 306, 321, 318, 471, 468, 483])
num_indices = design_freqs.size
pos_indices = np.arange(num_indices)
pos_indices[6:] += 1
pos_indices[2*6:] += 1
pos_indices[3*6:] += 1
pos_indices[4*6:] += 1
pos_indices[5*6:] += 1



mm = 1e-3
termination = 25.5*mm
reso_segment = 11*mm
pos = termination + (pos_indices + 0.5)*reso_segment
total_len = 2*termination + (num_indices + 5)*reso_segment
frac_position = pos/total_len

print (frac_position)


X, Y = np.meshgrid(fracs, cfs*1e-6)
all_Qes = np.array(all_Qes)
all_phics = np.angle(all_Qes)
all_Qcs = np.abs(all_Qes)/np.cos(all_phics)
all_f0s = np.array(all_f0s)
Z = np.abs(all_Qes)
Z2 = all_f0s
Z3 = all_phics*180/pi

interpolator = RectBivariateSpline(fracs, cfs*1e-6, all_Qcs)
predicted_Qcs = interpolator.ev(frac_position, design_freqs)
interpolator2 = RectBivariateSpline(fracs, cfs*1e-6, all_phics)
predicted_phics = interpolator.ev(frac_position, design_freqs)

pl.figure(figsize=(10,10))
pl.pcolor(X, Y, Z, cmap='jet')
pl.xlabel('Position along the line')
pl.ylabel('Resonance Frequency')
pl.colorbar()
#ax2 = pl.gca().twinx().twiny()
pl.scatter(frac_position, design_freqs, c='k')
pl.title("|Qe|")
pl.savefig("absQe_vs_linepos_vs_freq.png")

pl.figure(figsize=(10,10))
pl.pcolor(X, Y, all_Qcs, cmap='jet')
pl.xlabel('Position along the line')
pl.ylabel('Resonance Frequency')
pl.colorbar()
#ax2 = pl.gca().twinx().twiny()
pl.scatter(frac_position, design_freqs, c='k')
pl.title("Qc")
pl.savefig("Qc_vs_linepos_vs_freq.png")

pl.figure(figsize=(10,10))
pl.pcolor(X, Y, Z2, cmap='jet')
pl.xlabel('Position along the line')
pl.ylabel('Resonance Frequency')
pl.colorbar()
pl.scatter(frac_position, design_freqs, c='k')
pl.title("Predicted frequency shift")
pl.savefig("measfr_vs_linepos_vs_freq.png")

pl.figure(figsize=(10,10))
pl.pcolor(X, Y, Z3, cmap='jet')
pl.xlabel('Position along the line')
pl.ylabel('Resonance Frequency')
pl.colorbar()
pl.scatter(frac_position, design_freqs, c='k')
pl.title("angle Qe [deg]")
pl.savefig("angleQe_vs_linepos_vs_freq.png")


pl.figure(figsize=(10,10))
pl.scatter(design_freqs, predicted_Qcs, c='k')
pl.xlabel('Resonance Frequency [MHz]')
pl.ylabel('Qc')
pl.savefig('predQc_vs_reso_freq.png')

pl.figure(figsize=(10,10))
pl.scatter(frac_position, predicted_Qcs, c='k')
pl.xlabel('Resonance Frequency [MHz]')
pl.ylabel('Qc')
pl.savefig('predQc_vs_linepos.png')

pl.show()
data = np.vstack([design_freqs, frac_position, predicted_Qcs, predicted_phics]).T
np.savetxt('predicted_mismatchedQc.txt', data)
# I'm curious to see the positions of the resonators along the line
#freqs = np.linspace(290e6,450e6,1000)
#abcd_xline_tot = abcd_xline(freqs, xline_Z0, xline_n, xline_ltot)
#s21_xline_tot = s_of_abcd(abcd_xline_tot, 50)
#s_tot = s21_xline_tot[1,0]
#pl.figure(figsize=(10,10))
#pl.plot(freqs, np.abs(s_tot), 'b')
#abcd_xline_tot = abcd_xline(cfs, xline_Z0, xline_n, xline_ltot)
#s21_xline_tot = s_of_abcd(abcd_xline_tot, 50)
#s_tot = s21_xline_tot[1,0]
#pl.scatter(cfs, np.abs(s_tot), c='k')
#pl.show()
