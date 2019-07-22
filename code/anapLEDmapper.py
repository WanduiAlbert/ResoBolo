#!/usr/bin/env python import sys
import numpy as np
import matplotlib.pyplot as pl
from scipy import optimize, special
from scipy.constants import k, h
import os, sys
import reso_fit
import kid
import tempfitting
from analyze_tempsweep import FitTempSweep
import pdb
pi = np.pi

plot_diagnostic = True
diag_dir = 'fig/diagnostics'
frequencies = [305.8, 318.3, 337.4, 351.5, 414.6, 428.9, 434.1, 451.4, 473.4, 671.7]
warm_attenuation = 10
cold_attenuation = 60
total_attenuation = warm_attenuation + cold_attenuation	 # dB
# Extract Q vs. power dependence for resonators
# Data taken with HP8753C and autonetanal.py
if not os.path.exists('fig'):
	os.mkdir('fig')

if not os.path.exists('fig/diagnostics'):
	os.mkdir('fig/diagnostics')

if not os.path.exists('diagnostics'):
	os.mkdir('diagnostics')


def load_datachunk(fn):
        array_index = 0
	f = open(fn)
	words = f.readline().split()
	netanals = []
	temps = []
	powers = []
	assert words[0] == 'CONST'
	lines = f.readlines()
	for line in lines:
        words = line.split()
        if words[0] = 'CONST': continue
		words = line.split()
		if words[0] == 'OFF':
            yield array_index, np.average(temps), np.average(powers), np.array(netanals)
            netanals = []
            temps = []
            powers = []
            array_index += 1
            continue
		assert words[0] == 'netanal'
		datafn = 'raw/' + words[1] + '.dat'
		temp = float(words[2])
		power = float(words[5])
		data = np.loadtxt(datafn, skiprows=5)
		netanals.append(data)
		temps.append(temp)
		powers.append(power)
	return

def main():
	fn = sys.argv[1]
	datasets = load_datachunk(fn)
	temps = []
	powers = []
    frequencies = []
	sigma_frequencies = []
    #netanals = []
    internal_qs = []
    sigma_dqis = []


	for array_i, t, p, n in datasets:
		if t < 0: t = 0.08
        nres, npt, _ = netanals.shape
        f0s = np.zeros(nres)
        Qis = np.zeros(nres)
        Qrs = np.zeros(nres)
        sigma_fr = np.zeros(nres)
        sigma_dqi = np.zeros(nres)
        fs = netanals[:,:,0]
        re = netanals[:,:,1]
        im = netanals[:,:,2]
        mag = re*re+im*im
        delay = 42.0e-9
        phase_corr = 1e6*fs*delay*2.*pi
        phi = np.arctan2(im,re)+phase_corr
        phi = np.unwrap(phi,axis=-1)
        z = np.sqrt(mag)*np.exp(1.j*phi)
        re = z.real
        im = z.imag

        for ires in range(nres):
            f = fs[ires,:]
            x = re[ires,:]
            y = im[ires,:]
            imag = mag[itemp,ires,:]
            dbmag = 10*np.log10(imag)
            try:
                f0,_,_,_,_,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov =\
                    reso_fit.do_fit(f*1e6,x,y,get_cov=True)
                sigmas = np.diag(pcov)**0.5
            except RuntimeError:
                print ("Error obtaining the fit")
                continue
            Qe = Qe_re + 1j*Qe_im
            phi_c = np.arctan2(Qe_im, Qe_re)
            dQe = 1./Qe
            Qc = 1./np.real(dQe)
            print(f0,Qi,Qr,Qe_re)
            f0s[ires] = f0
            Qis[ires] = Qi
            Qrs[ires] = Qr
            sigma_fr[ires] = sigmas[0]
            sigma_dqi[ires] = np.sqrt(sigmas[5]**2 + sigmas[6]**2)
            mmag = mre*mre + mim*mim
            mdbmag = 10*np.log10(mmag)
            freq = "%3.1f" %(frequencies[ires]) + 'MHz'

            if plot_diagnostic:
                pl.figure(923 + ires + int(temps[itemp]*1e3))
                pl.plot(f,dbmag, 'b', label='data')
                pl.xlim(right=max(f)+0.25)
                pl.ylim(top=max(dbmag)+0.5)
                pl.plot(f,mdbmag,'r', label='fit')
                pl.grid()
                pl.xlabel('Frequency (MHz)')
                pl.ylabel('S21')
                pl.title('Resonator '+ freq + ' Array Index %d'%(array_i))
                pl.legend(loc='upper right')
                path = os.path.join(diag_dir,'reso_'+freq + '_array_index_%d.png'%(array_i))
                #pl.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
                pl.savefig(path)
                pl.close()

        temps.append(t)
        powers.append(p)
        frequencies.append(f0s)
        sigma_frequencies.append(sigma_fr)
        #netanals.append()
        internal_qs.append(Qis)
        sigma_dqis.append(sigma_dqi)


        temps = np.array(temps)
        powers = np.array(powers)
        frequencies = np.array(frequencies)
        sigma_frequencies = np.array(sigma_frequencies)
        internal_qs = np.array(internal_qs)
        sigma_dqis = np.array(sigma_dqis)

    df = (frequencies - frequencies[0])/frequencies[0]
    dQ = (1./internal_qs - 1./internal_qs[0])

    array_N, _ = frequencies.shape

    for i in range(array_N):
        fig, ax = plt.subplots()
        ax.hist(df[i], nbins=nres, histtype='bar')
        ax.set_xlabel('Resonator Index')
        ax.set_ylabel('dfr/fr')
        ax.set_title('Array Index %d'%(i))
        ax.grid()
        pl.savefig('fig/df_shift_array_index_%d'%(i))


if __name__=='__main__':
	main()
