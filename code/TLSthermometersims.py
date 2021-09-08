#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import pi
from scipy import interpolate
from scipy.constants import h, k, c
from scipy import interpolate, optimize
import scipy.signal as signal
import glob
import sys,os
import pdb
#sys.path.append("/home/wanduialbert/bicep/code/instruments/copper")
import reso_fit
import meshanalysisofresonator as mesh
pi = np.pi

MHz = 1e6
GHz = 1e9
kHz = 1e3
pW = 1e-12
um = 1e-6

nH = 1e-9
pF = 1e-12
MHz = 1e6
Z0 = 50
Y0 = 1./Z0

datadir = '../numerical_sims/'
plotdir = 'TLSthermometer_fig/'

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'][1:]

plot_diagnostic = True
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
doS21analysis = False

def load_data(fn, nports=1, paramtype='Y'):
    p = paramtype
    if nports==1:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
    elif nports==2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    elif nports==3:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "13", np.complex128), (p + "21", np.complex128),\
            (p + "22", np.complex128), (p + "23", np.complex128), (p +"31", np.complex128),\
            (p + "32", np.complex128), (p  +"33", np.complex128)])
    elif nports==4:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "13", np.complex128), (p + "14", np.complex128),\
            (p + "21", np.complex128), (p + "22", np.complex128), (p + "23", np.complex128),\
            (p + "24", np.complex128), (p +"31", np.complex128), (p + "32", np.complex128),\
            (p + "33", np.complex128), (p  +"34", np.complex128), (p +"41", np.complex128),\
            (p + "42", np.complex128), (p + "43", np.complex128), (p  +"44", np.complex128)])
    with open(fn) as fp:
        for i in range(8): fp.readline()
        ildthickness = float(fp.readline().split(" ")[-1])
        ilddielectricconst = float(fp.readline().split(" ")[-1])
        #print (distfromedge)
        dataset = []
        freqs = []
        while True:
            line = fp.readline()
            #print (line)
            #print (line.startswith(" "))
            #print (line)
            if line is None: break
            #if line is "":
            #    pdb.set_trace()
            if line.startswith(" ") or line is "":
                #tup = list(zip(dataset[:, 0], p11, p12, p21, p22))
                if dataset:
                    dataset = np.asarray(dataset)
                    if nports==1:
                        tup = list(zip(np.array(freqs), dataset))
                    elif nports==2:
                        p11 = dataset[:, 0]
                        p12 = dataset[:, 1]
                        p21 = dataset[:, 2]
                        p22 = dataset[:, 3]
                        tup = list(zip(np.array(freqs), p11, p12, p21, p22))
                    elif nports==3:
                        p11 = dataset[:, 0]
                        p12 = dataset[:, 1]
                        p13 = dataset[:, 2]
                        p21 = dataset[:, 3]
                        p22 = dataset[:, 4]
                        p23 = dataset[:, 5]
                        p31 = dataset[:, 6]
                        p32 = dataset[:, 7]
                        p33 = dataset[:, 8]
                        tup = list(zip(np.array(freqs), p11, p12, p13, p21, p22, p23, p31, p32, p33))
                    elif nports==4:
                        p11 = dataset[:, 0]
                        p12 = dataset[:, 1]
                        p13 = dataset[:, 2]
                        p14 = dataset[:, 3]
                        p21 = dataset[:, 4]
                        p22 = dataset[:, 5]
                        p23 = dataset[:, 6]
                        p24 = dataset[:, 7]
                        p31 = dataset[:, 8]
                        p32 = dataset[:, 9]
                        p33 = dataset[:, 10]
                        p34 = dataset[:, 11]
                        p41 = dataset[:, 12]
                        p42 = dataset[:, 13]
                        p43 = dataset[:, 14]
                        p44 = dataset[:, 15]
                        tup = list(zip(np.array(freqs), p11, p12, p13, p14,
                            p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43,
                            p44))
                    #tup = list(zip(np.array(freqs), dataset))
                    yield ildthickness, ilddielectricconst, np.array(tup, dtype=dtypes)
                if line is "": break
                dataset = []
                for i in range(7): fp.readline()
                ildthickness = float(fp.readline().split(" ")[-1])
                ilddielectricconst = float(fp.readline().split(" ")[-1])
                continue
            if nports==1:
                freq, re, im = list(map(lambda x: float(x), line.split(",")))
                freqs.append([freq])
                dataset.append([re + 1j*im])
            elif nports==2:
                freq,\
                re11,im11,\
                re12, im12,\
                re21, im21,\
                re22, im22 = list(map(
                    lambda x: float(x), line.split(",")))
                freqs.append([freq])
                dataset.append([
                        re11 + 1j*im11,\
                        re12 + 1j*im12,\
                        re21 + 1j*im21,\
                        re22 + 1j*im22])
            elif nports==3:
                freq,\
                re11,im11, re12, im12, re13, im13,\
                re21, im21, re22, im22, re23, im23,\
                re31, im31, re32, im32, re33, im33 = list(map(
                    lambda x: float(x), line.split(",")))
                freqs.append([freq])
                dataset.append([
                        re11 + 1j*im11, re12 + 1j*im12, re13 + 1j*im13,\
                        re21 + 1j*im21, re22 + 1j*im22, re23 + 1j*im23,\
                        re31 + 1j*im31, re32 + 1j*im32, re33 + 1j*im33])
            elif nports==4:
                freq,\
                re11,im11, re12, im12, re13, im13, re14, im14,\
                re21, im21, re22, im22, re23, im23, re24, im24,\
                re31, im31, re32, im32, re33, im33, re34, im34,\
                re41, im41, re42, im42, re43, im43, re44, im44 = list(map(
                    lambda x: float(x), line.split(",")))
                freqs.append([freq])
                dataset.append([
                        re11 + 1j*im11, re12 + 1j*im12, re13 + 1j*im13, re14 + 1j*im14,\
                        re21 + 1j*im21, re22 + 1j*im22, re23 + 1j*im23, re24 + 1j*im24,\
                        re31 + 1j*im31, re32 + 1j*im32, re33 + 1j*im33, re34 + 1j*im34,\
                        re41 + 1j*im41, re42 + 1j*im42, re43 + 1j*im43, re44 + 1j*im44])
    return




def load_data_single(fn, nports=2, paramtype='S'):
    dataset = np.loadtxt(fn, skiprows=10, delimiter=',')
    p = paramtype
    if nports==1:
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        tup = list(zip(dataset[:, 0], p11))
        return np.array(tup, dtype=dtypes)
    elif nports==2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        p12 = dataset[:, 3] + 1j*dataset[:, 4]
        p21 = dataset[:, 5] + 1j*dataset[:, 6]
        p22 = dataset[:, 7] + 1j*dataset[:, 8]
        tup = list(zip(dataset[:, 0], p11, p12, p21, p22))
        return np.array(tup, dtype=dtypes)
    elif nports==3:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "13", np.complex128), (p + "21", np.complex128),\
            (p + "22", np.complex128), (p + "23", np.complex128), (p +"31", np.complex128),\
            (p + "32", np.complex128), (p  +"33", np.complex128)])
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        p12 = dataset[:, 3] + 1j*dataset[:, 4]
        p13 = dataset[:, 5] + 1j*dataset[:, 6]
        p21 = dataset[:, 7] + 1j*dataset[:, 8]
        p22 = dataset[:, 9] + 1j*dataset[:,10]
        p23 = dataset[:,11] + 1j*dataset[:,12]
        p31 = dataset[:,13] + 1j*dataset[:,14]
        p32 = dataset[:,15] + 1j*dataset[:,16]
        p33 = dataset[:,17] + 1j*dataset[:,18]
        tup = lipt(zip(dataset[:, 0], p11, p12, p13, p21, p22, p23, p31, p32, p33))
        return np.array(tup, dtype=dtypes)
    elif nports==4:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "13", np.complex128), (p + "14", np.complex128),\
            (p + "21", np.complex128), (p + "22", np.complex128), (p + "23", np.complex128),\
            (p + "24", np.complex128), (p +"31", np.complex128), (p + "32", np.complex128),\
            (p + "33", np.complex128), (p  +"34", np.complex128), (p +"41", np.complex128),\
            (p + "42", np.complex128), (p + "43", np.complex128), (p  +"44", np.complex128)])
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        p12 = dataset[:, 3] + 1j*dataset[:, 4]
        p13 = dataset[:, 5] + 1j*dataset[:, 6]
        p14 = dataset[:, 7] + 1j*dataset[:, 8]
        p21 = dataset[:, 9] + 1j*dataset[:, 10]
        p22 = dataset[:, 11] + 1j*dataset[:, 12]
        p23 = dataset[:, 13] + 1j*dataset[:, 14]
        p24 = dataset[:, 15] + 1j*dataset[:, 16]
        p31 = dataset[:, 17] + 1j*dataset[:, 18]
        p32 = dataset[:, 19] + 1j*dataset[:, 20]
        p33 = dataset[:, 21] + 1j*dataset[:, 22]
        p34 = dataset[:, 23] + 1j*dataset[:, 24]
        p41 = dataset[:, 25] + 1j*dataset[:, 26]
        p42 = dataset[:, 27] + 1j*dataset[:, 28]
        p43 = dataset[:, 29] + 1j*dataset[:, 30]
        p44 = dataset[:, 31] + 1j*dataset[:, 32]
        tup = list(zip(dataset[:, 0], p11, p12, p13, p14,
            p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43,
            p44))
        return np.array(tup, dtype=dtypes)

if __name__=="__main__":
    filenames = glob.glob(datadir + "TKID_Module_Thermometer_*.csv")
    filenames = [filenames[index] for index in [4, 5, 7, 8, 9, 10, 12, 13]]
    #for fn in filenames:
    #    print (fn)
    #exit()
    frs = np.zeros(len(filenames))
    Qrs = np.zeros(len(filenames))
    Qis = np.zeros(len(filenames))
    Qcs = np.zeros(len(filenames))
    phics = np.zeros(len(filenames))
    labels = ['400umGPpad', 'morecoupling', 'wlosstangent', 'symmetric', 'samesidecoupling',\
            '12um_boundary', 'lateralcoupcaps', 'addedline']
    fr_guess = np.array([927.1, 876.7, 684.4, 878.5, 747.1, 1057.3, 748.0, 684.4])

    plt.figure(1, figsize=(10,10))
    plt.figure(2, figsize=(10,10))
    for isim,fn in enumerate(filenames):
        print (isim, fn)
        section = labels[isim]
        Sparams = load_data_single(fn, nports=2, paramtype='S')
        f = Sparams['frequency']
        re = Sparams['S21'].real
        im = Sparams['S21'].imag
        funcre = interpolate.interp1d(f, re, kind='cubic')
        funcim = interpolate.interp1d(f, im, kind='cubic')

        finterp = np.r_[f[0]:f[-1]:5000j]
        reinterp = funcre(finterp)
        iminterp = funcim(finterp)
        maginterp = np.sqrt(reinterp*reinterp + iminterp*iminterp)
        dBmaginterp = 20*np.log10(maginterp)

        mag = np.sqrt(re*re + im*im)
        z = re + 1j*im
        dbmag = 20*np.log10(mag)

        mask = np.ones_like(f, dtype=bool)
        mask[f > fr_guess[isim]+1.5] = False
        mask[f < fr_guess[isim]-1.5] = False
        midpt = f.size //2  + 1
        ffine = np.r_[f[mask][0]:f[mask][-1]:10000j]
        try:
            f0,A,m,phi,D,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov =\
            reso_fit.do_fit(f[mask]*1e6,re[mask],im[mask], get_cov=True)
            sigmas = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print ("Error obtaining the fit")
            #f0 = frequencies[isim]
            continue
        Qe = Qe_re + 1j * Qe_im
        dQe = 1/Qe
        phi_c = np.arctan2(Qe_im, Qe_re)
        Qc = 1./np.real(dQe)
        popt = [f0,A,m,phi,D,1./Qr,dQe.real,dQe.imag,a]
        print("Best fit results: ", section, f0, Qi, Qc, phi_c)
        frs[isim] = f0
        Qcs[isim] = Qc
        phics[isim] = phi_c
        Qrs[isim] = Qr
        Qis[isim] = Qi
        mmag = mre*mre + mim*mim
        mdbmag = 10*np.log10(mmag)

        if plot_diagnostic:
            plt.figure(isim)
            plt.plot(f[mask],dbmag[mask], 'bs', ms=12,\
                    label='Qr=%d Qi=%d Qc=%d phi_c=%1.3f'%(Qr,Qi,Qc, phi_c))
            plt.plot(ffine,mdbmag,'r', label='fit')
            plt.xlim(left=f0-0.5, right=f0+0.5)
            plt.grid()
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('|S21|')
            lgd = plt.legend(loc="upper left", bbox_to_anchor=(1.0,1.0))
            path = os.path.join(plotdir,'reso_%s'%section + ".png")
            plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
            #plt.show()
            plt.close()

            plt.figure(300 + isim )
            scale = np.max(mmag)**0.5
            plt.plot(re[mask]/scale,im[mask]/scale, 'bs', ms=12, label='data')
            plt.plot(mre/scale,mim/scale,'r', label='fit')
            plt.grid()
            plt.axis('square')
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.legend(loc='upper right')
            path = os.path.join(plotdir,'reso_IQ_%s'%section + ".png")
            plt.savefig(path)
            #plt.show()
            plt.close()


