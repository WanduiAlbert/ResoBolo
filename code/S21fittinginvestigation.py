import numpy as np
import matplotlib.pyplot as pl
import reso_fit
import emcee
import corner
from scipy.interpolate import RectBivariateSpline
pi = np.pi
from scipy import optimize
from scipy.constants import c, epsilon_0, k

cmap = pl.get_cmap('jet')

MHz = 1e6
all_Qes = []
all_f0s = []
noise_temp = 5
IFBW = 1e3
PgdBm = -110#dBm
Pg = 1e-3*10**(PgdBm/10)
noisestdev = np.sqrt(k*noise_temp/(2*Pg))*IFBW
ncf = 101
ntries = 101
cfs = np.linspace(290e6,500e6,ncf)
f0k = cfs[0]
#Qc = 20e3
#Qcs = Qc*np.ones(ncf)
Qcs = np.linspace(20e3,100e3,ncf)
phi_csim = 0.0
phics = phi_csim*np.ones(ncf)
#phics = np.linspace(-pi,pi,ncf)

#delay = 45.0e-9
delay = 0.0
fitted_Qcs = np.zeros((ncf, ntries))
fitted_phics = np.zeros((ncf, ntries))
fitted_as = np.zeros((ncf, ntries))
#Qcsim = Qc
def dbmag(x):
	return 10*np.log10(np.abs(x))

for k in range(ncf):
	print (k)
	#f0k = cfs[k]
	Qcsim = Qcs[k]
	phi_csim = phics[k]
	Qi = 100000
	Qr = 1/(1/Qi + 1/Qcsim)
	dQe_im = np.tan(phi_csim)/Qcsim
	p0 = [f0k/MHz,1,0,0,delay,1/Qr,1./Qcsim,dQe_im,0.0]
	df = 0.5e6#1*f0k/Qr
	fs = np.linspace(f0k-df,f0k+df,1001)
	phase_corr = 1e6*fs*delay*2.*pi
	s21sim = reso_fit.complex_of_real(
			reso_fit.model_linear_wslope(
				fs,f0k/MHz,1,0,0,delay,1/Qr,1./Qcsim,dQe_im,0.0))

	allpopts = []
	for atry in range(ntries):
		#noise = np.random.randn(fs.size)*noisestdev + 1j*np.random.randn(fs.size)*noisestdev
		#s21sim += noise


		#s21simdB = dbmag(s21sim)
		#pl.plot(fs/MHz, s21simdB)
		#pl.show()
		#continue

		x = s21sim.real + np.random.randn(fs.size)*noisestdev
		y = s21sim.imag + np.random.randn(fs.size)*noisestdev
		mag = (x*x + y*y)**0.5
		magdB = 10*np.log10(mag)

		# Lets try and do the fitting now
		f0,A,m,phi,D,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov = reso_fit.do_fit(fs,x,\
				y, get_cov=True)
		Qe = Qe_re + 1j * Qe_im
		dQe = 1/Qe
		popt = np.array([f0,A,m,phi,D,1./Qr,dQe.real,dQe.imag,a])
		allpopts.append(popt)
		#popt[0] /= 1e6
		#pcov[0, :] /= 1e6
		#pcov[:, 0] /= 1e6
		Qi = abs(Qi)
		Qr = abs(Qr)
		Qc = 1./np.real(dQe)
		phi_c = np.arctan2(Qe_im, Qe_re)
		fitted_Qcs[k, atry] = Qc
		fitted_phics[k, atry] = phi_c
		fitted_as[k, atry] = a
		print ("Qc from sim, ", Qcsim)
		print ("Qc from fit, ", Qc)
		print ("phic from fit, ", phi_c)
		print ("a from fit, ", a)
		mmag = (mre*mre + mim*mim)**0.5
		mmagdB = 10*np.log10(mmag)

		#pl.figure()
		#pl.plot(fs/1e6, magdB, label='Data')
		#pl.plot(fs/1e6, mmagdB, label='Fit')
		#pl.grid()
		#pl.xlabel('Frequency [MHz]')
		#pl.ylabel('|S21|')
		#pl.savefig("cornerplots/%3.2fMHz_fitvsdata_Qc%d.png"%(f0k/MHz, Qc))
		#pl.close()

		#samples = np.random.multivariate_normal(popt, pcov, 2000000)
		#samples_useful = []
		#for i in range(samples.shape[0]):
		#	if samples[i, -1] < 0: continue
		#	samples_useful.append(samples[i])
		#samples = np.array(samples_useful)
		#print (samples.shape)
		#ndim, nwalkers = len(popt), 100
		#pos = [popt + 3*np.sqrt(np.diag(pcov))*np.random.randn(ndim) for i in range(nwalkers)]

		#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(temps,
		#	freqs, sig_freqs*100))
		#sampler.run_mcmc(pos, 1800)
		#samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
		#sampler.reset()

		print (p0)
		print (popt)
		print (np.sqrt(np.diag(pcov)))

		ndim = len(popt)

	allpopts = np.array(allpopts)
	header = "f_r [MHz] A m phi D delta_r delta_ere delta_eim a"
	np.savetxt("cornerplots/%3.2fMHz_popt_Qc%d.txt"%(f0k/MHz, Qcsim), allpopts,
			delimiter='    ', header=header)
	print (allpopts.shape)
	popt_avg = np.mean(allpopts, axis=0)
	pcov_avg = np.cov(allpopts.T)
	samples = np.random.multivariate_normal(popt_avg, pcov_avg, 2000000)
	samples_useful = []
	for i in range(samples.shape[0]):
		if samples[i, -1] < 0: continue
		samples_useful.append(samples[i])
	samples = np.array(samples_useful)
	print (samples.shape)
	figlabels = [r"$f_r$", "$A$", "$m$", "$\phi$", "$D$", "$\delta_r$",
			"$\delta_{e,re}$", "$\delta_{e,im}$", "$a$"]
	fig = corner.corner(samples, labels=figlabels, truths=popt_avg,
			quantiles=(0.16, 0.84), levels=(1 - np.exp(-0.5),), show_titles=True, fontsize=12,
			title_fmt=".4e", plot_contours=False, title_kwargs={"fontsize":20})

	# Extract the axes
	axes = np.array(fig.axes).reshape((ndim, ndim))

	# Loop over the diagonal
	for i in range(ndim):
		ax = axes[i, i]
		ax.axvline(p0[i], color="r")

	# Loop over the histograms
	for yi in range(ndim):
		for xi in range(yi):
			ax = axes[yi, xi]
			ax.axvline(p0[xi], color="r")
			ax.axhline(p0[yi], color="r")
			ax.plot(p0[xi], p0[yi], "sr")

	fig.savefig("cornerplots/%3.2fMHz_triangle_Qc%d.png"%(f0k/MHz, Qcsim))
	#pl.show()
	pl.close()
fitted_Qcavg = np.mean(fitted_Qcs, axis=1)
fitted_Qcstd = np.std(fitted_Qcs, axis=1)
Qcbias = Qcs - fitted_Qcavg
pl.figure(1)
pl.errorbar(Qcs, Qcbias/Qcs, yerr=fitted_Qcstd/Qcs, fmt='ko')
pl.grid()
pl.xlabel(r'Predicted Qc')
pl.ylabel(r'(Predicted Qc - Fitted Qc)/Predicted_Qc')
pl.savefig('cornerplots/fitted_vs_predicted_Qcs.png')
pl.show()

fitted_phicavg = np.mean(fitted_phics, axis=1)
fitted_phicstd = np.std(fitted_phics, axis=1)
phi_bias = np.angle(np.exp(1j*(phics - fitted_phicavg)))
pl.figure(2)
pl.errorbar(Qcs, phi_bias, yerr=fitted_phicstd, fmt='ko')
pl.grid()
pl.xlabel(r'Predicted Qc')
pl.ylabel(r'Fitted phic - predicted phic')
pl.savefig('cornerplots/fitted_vs_predicted_phics.png')
pl.show()

fitted_aavg = np.mean(fitted_as, axis=1)
fitted_astd = np.std(fitted_as, axis=1)
pl.figure(2)
pl.errorbar(Qcs, fitted_aavg, yerr=fitted_astd, fmt='ko')
pl.grid()
pl.xlabel(r'Predicted Qc')
pl.ylabel(r'a from fit')
pl.savefig('cornerplots/fitted_vs_predicted_as.png')
pl.show()
