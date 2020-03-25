#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import pi
from scipy import interpolate
from scipy.constants import h, k, c
from scipy import interpolate
import scipy.signal as signal
import glob
import sys,os
import pdb
#sys.path.append("/home/wanduialbert/bicep/code/instruments/copper")
import reso_fit
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
plotdir = 'extractQc_fig/'

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'][1:]

plot_diagnostic = True
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

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
    dataset = np.loadtxt(fn, skiprows=8, delimiter=',')
    p = paramtype
    if nports==2:
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


# Connect matched 50 Ohm loads to ports 3 and 4 of the capacitor
def reduce_Yparameters(capY):
    delta = (Y0 + capY['Y33'])*(Y0 + capY['Y44']) - capY['Y43']*capY['Y34']
    deltap = capY['Y33']*capY['Y44'] - capY['Y43']*capY['Y34']


    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    p11 = capY['Y11']
    #p11 += (capY['Y34']*capY['Y41'] - capY['Y31']*(Y0 + capY['Y44']))/delta
    p11 += (capY['Y34']*capY['Y41'] - capY['Y31']*capY['Y44'])/deltap
    p12 = capY['Y12']
    #p12 += (capY['Y34']*capY['Y42'] - capY['Y32']*(Y0 + capY['Y44']))/delta
    p12 += (capY['Y34']*capY['Y42'] - capY['Y32']*capY['Y44'])/deltap
    p21 = capY['Y21']
    #p21 += (capY['Y43']*capY['Y31'] - capY['Y41']*(Y0 + capY['Y33']))/delta
    p21 += (capY['Y43']*capY['Y31'] - capY['Y41']*capY['Y33'])/deltap
    p22 = capY['Y22']
    #p22 += (capY['Y43']*capY['Y32'] - capY['Y42']*(Y0 + capY['Y33']))/delta
    p22 += (capY['Y43']*capY['Y32'] - capY['Y42']*capY['Y33'])/deltap
    tup = list(zip(capY['frequency'], p11, p12, p21, p22))
    return np.array(tup, dtype=dtypes)



def combine_Sparameters(capS, indS):
    alpha = 1 - indS['S11']*capS['S11'] - indS['S12']*capS['S21']
    beta = -indS['S11']*capS['S12'] - indS['S12']*capS['S22']
    gamma = -indS['S21']*capS['S11'] - indS['S22']*capS['S21']
    delta = 1 - indS['S21']*capS['S12'] - indS['S22']*capS['S22']

    A = indS['S11']*capS['S13'] + indS['S12']*capS['S23']
    B = indS['S11']*capS['S14'] + indS['S12']*capS['S24']
    C = indS['S21']*capS['S13'] + indS['S22']*capS['S23']
    D = indS['S21']*capS['S14'] + indS['S22']*capS['S24']

    p = 'S'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    p11 = capS['S33']
    p11 += capS['S31']*(delta*A-beta*C)/(delta*alpha-beta*gamma)
    p11 += capS['S32']*(gamma*A-alpha*C)/(beta*gamma - delta*alpha)
    p12 = capS['S34']
    p12 += capS['S31']*(delta*B-beta*D)/(delta*alpha-beta*gamma)
    p12 += capS['S32']*(gamma*B-alpha*D)/(beta*gamma - delta*alpha)
    p21 = capS['S43']
    p21 += capS['S41']*(delta*A-beta*C)/(delta*alpha-beta*gamma)
    p21 += capS['S42']*(gamma*A-alpha*C)/(beta*gamma - delta*alpha)
    p22 = capS['S44']
    p22 += capS['S41']*(delta*B-beta*D)/(delta*alpha-beta*gamma)
    p22 += capS['S42']*(gamma*B-alpha*D)/(beta*gamma - delta*alpha)
    tup = list(zip(capS['frequency'], p11, p12, p21, p22))
    return np.array(tup, dtype=dtypes)


def get_Yin_parallel(capY, indY):

    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
    p11 = (capY['Y11'] + indY['Y11']) - (capY['Y21'] + indY['Y21'])
    tup = list(zip(capY['frequency'], p11))
    return np.array(tup, dtype=dtypes)

def get_Yin_series(capY, indY):
    capZ = convert_Y_to_Z(capY, nports=2)
    indZ = convert_Y_to_Z(indY, nports=2)
    p = 'Z'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    p11 = capZ['Z11'] + indZ['Z11']
    p12 = capZ['Z12'] + indZ['Z12']
    p21 = capZ['Z21'] + indZ['Z21']
    p22 = capZ['Z22'] + indZ['Z22']
    tup = list(zip(capY['frequency'], p11, p12, p21, p22))
    seriesZ =  np.array(tup, dtype=dtypes)
    seriesY = convert_Z_to_Y(seriesZ, nports=2)

    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
    p11 = seriesY['Y11'] - seriesY['Y21']
    tup = list(zip(capY['frequency'], p11))
    return np.array(tup, dtype=dtypes)


def resonator_1port(capS, indS):

    numerator = 1 - indS['S22']*capS['S22'] - indS['S12']*capS['S21']
    denominator = indS['S11'] - indS['S11']*capS['S22']*indS['S22'] + indS['S12']*capS['S22']*indS['S22']

    p = 'S'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
    p11 = numerator/denominator
    tup = list(zip(capS['frequency'], p11))
    return np.array(tup, dtype=dtypes)


def convert_Z_to_Y(Zparams, nports=2):
    p = 'Y'
    if nports == 1:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
        p11 = 1./Zparams['Z11']
        tup = list(zip(Zparams['frequency'], p11))
    elif nports == 2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
        delta = Zparams['Z11']*Zparams['Z22'] - Zparams['Z21']*Zparams['Z12']

        p11 = Zparams['Z22']/delta
        p12 = -Zparams['Z12']/delta
        p21 = -Zparams['Z21']/delta
        p22 = Zparams['Z11']/delta
        tup = list(zip(Zparams['frequency'], p11, p12, p21, p22))
    return np.array(tup, dtype=dtypes)


def convert_Y_to_Z(Yparams, nports=2):
    p = 'Z'
    if nports == 1:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
        p11 = 1./Yparams['Y11']
        tup = list(zip(Yparams['frequency'], p11))
    elif nports==2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
        delta = Yparams['Y11']*Yparams['Y22'] - Yparams['Y21']*Yparams['Y12']

        p11 = Yparams['Y22']/delta
        p12 = -Yparams['Y12']/delta
        p21 = -Yparams['Y21']/delta
        p22 = Yparams['Y11']/delta
        tup = list(zip(Yparams['frequency'], p11, p12, p21, p22))
    return np.array(tup, dtype=dtypes)

#fns = glob.glob(datadir + "*_Sparams.csv")
#frequencies = list(map(lambda x: int(os.path.split(x)[-1].split(".")[0][:3]), fns))
#for ires,fn in enumerate(fns):
#    Sparams = load_data(fn)
#    f = Sparams['frequency']
#    re = Sparams['S21'].real
#    im = Sparams['S21'].imag
#    funcre = interpolate.interp1d(f, re)
#    funcim = interpolate.interp1d(f, im)
#
#    finterp = np.r_[f[0]:f[-1]:5000j]
#    reinterp = funcre(finterp)
#    iminterp = funcim(finterp)
#    maginterp = np.sqrt(reinterp*reinterp + iminterp*iminterp)
#    dBmaginterp = 10*np.log10(maginterp)
#
#    #f = finterp
#    #re = reinterp
#    #im = iminterp
#    mag = np.sqrt(re*re + im*im)
#    delay = 0.0e-9
#    phase_corr = 1e6*f*delay*2.*pi
#    phi = np.arctan2(im,re)+phase_corr
#    phi = np.unwrap(phi,axis=-1)
#    z = np.sqrt(mag)*np.exp(1.j*phi)
#    re = z.real
#    im = z.imag
#    dbmag = 10*np.log10(mag)
#
#    #mask = f < 331.2
#    mask = np.ones_like(f, dtype=bool)
#    midpt = f.size //2  + 1
#    #mask[89:90] = False
#    try:
#        f0,A,m,phi,D,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov =\
#        reso_fit.do_fit(f[mask]*1e6,re[mask],im[mask], get_cov=True)
#        sigmas = np.sqrt(np.diag(pcov))
#    except RuntimeError:
#        print ("Error obtaining the fit")
#        #f0 = frequencies[ires]
#        continue
#    Qe = Qe_re + 1j * Qe_im
#    dQe = 1/Qe
#    phi_c = np.arctan2(Qe_im, Qe_re)
#    Qc = 1./np.real(dQe)
#    popt = [f0,A,m,phi,D,1./Qr,dQe.real,dQe.imag,a]
#    print(fn.split(".")[0], f0, Qi, Qc)
#    res_re = mre - re[mask]
#    res_im = mim - im[mask]
#    res_mag = res_re*res_re + res_im*res_im
#    resdbmag = 10*np.log10(res_mag)
#    mmag = mre*mre + mim*mim
#    mdbmag = 10*np.log10(mmag)
#    freq = frequencies[ires]
#    nmask = mask
#
#    if plot_diagnostic:
#        plt.figure(ires )
#        plt.plot(f[mask],mdbmag,'r', label='fit')
#        plt.plot(f[mask],dbmag[mask], 'b', label='Qr=%d Qi=%d Qc=%d phi_c=%1.3f'%(Qr,Qi,Qc, phi_c))
#        #plt.plot(finterp, dBmaginterp, 'k')
#        #plt.ylim(top=max(dbmag)+0.5, bottom=min(dbmag)-12.5)
#        plt.grid()
#        plt.xlabel('Frequency (MHz)')
#        plt.ylabel('|S21|')
#        lgd = plt.legend(loc="upper left", bbox_to_anchor=(1.0,1.0))
#        path = os.path.join(plotdir,'reso_%d'%freq + "MHz.png")
#        #plt.savefig(path)
#        plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
#        plt.close()
#
#        plt.figure(300 + ires )
#        scale = np.max(mmag)**0.5
#        plt.plot(re/scale,im/scale, 'b', label='data')
#        plt.plot(mre/scale,mim/scale,'r', label='fit')
#        plt.grid()
#        plt.axis('square')
#        plt.xlabel('I')
#        plt.ylabel('Q')
#        plt.legend(loc='upper right')
#        path = os.path.join(plotdir,'reso_IQ_%d'%freq + "MHz.png")
#        plt.savefig(path)
#        plt.close()

if __name__=="__main__":
    capfilename = datadir + "Cap_300MHz_with_boundary_and_coupcaps_variableILDthickness_Yparams.csv"
    indfilename = datadir + "inductor_with_GPboundary_variableILDthickness_Yparams.csv"
    capdataloader =  load_data(capfilename, nports=4, paramtype='Y')
    inddataloader =  load_data(indfilename, nports=2, paramtype='Y')
    ildthicknesses = []
    ilddielectricconsts = []
    capdataset = []
    inddataset = []
    print ("Loading data for the capacitor")
    for t, epsilon_r, capdata in capdataloader:
        ildthicknesses.append(t)
        ilddielectricconsts.append(epsilon_r)
        capdataset.append(capdata)

    print ("\n\n\n")
    print ("Loading data for the inductor")
    for t, epsilon_r, inddata in inddataloader:
        inddataset.append(inddata)

    ildthicknesses = np.asarray(ildthicknesses)
    ilddielectricconsts = np.asarray(ilddielectricconsts)
    #capdataset = np.asarray(capdataset)
    #inddataset = np.asarray(inddataset)

    print (ildthicknesses)
    print (ilddielectricconsts)
    #print (capdataset.shape)
    #print (inddataset.shape)


    N = ildthicknesses.size
    frs = np.zeros(N)
    Qs = np.zeros(N)

    #plt.figure(1, figsize=(10,10))
    plt.figure(2, figsize=(10,10))
    for i in range(N):
        #y13 = np.mean(capdataset[i]['Y13'])
        #y14 = np.mean(capdataset[i]['Y14'])
        #y23 = np.mean(capdataset[i]['Y23'])
        #y24 = np.mean(capdataset[i]['Y24'])
        #print (y13.real, y13.imag)
        #print (y14.real, y14.imag)
        #print (y23.real/y13.real, y23.imag)
        #print (y24.real/y14.real, y24.imag)
        #print ("\n\n")
        capY = reduce_Yparameters(capdataset[i])
        finalY = get_Yin_parallel(capY, inddataset[i])
        #finalY = get_Yin_series(capY, inddataset[i])
        Yin = finalY['Y11']
        #finalS = combine_Sparameters(capdataset[i], inddataset[i])
        #Gammain = resonator_1port(capdataset[i], inddataset[i])
        #Yin = 1./Yin['Y11']
        f = capY['frequency']
        f_fine = np.r_[200:400:1000j]
        mask = f > 200
        mask[f > 400] = False
        tck = interpolate.splrep(f[mask], Yin[mask].imag, s=0)
        retck = interpolate.splrep(f[mask], Yin[mask].real, s=0)
        fr, = interpolate.sproot(tck)
        slope = interpolate.splev(fr, tck, der=1)
        dR,  = interpolate.splev([fr], retck)
        R = 1./dR
        Q = slope*R*fr/2
        print (fr, Q)
        frs[i] = fr
        Qs[i] = Q
        #plt.figure(1)
        #p = plt.plot(finalS['frequency'], np.abs(finalS['S21']), ls='None',
        #        ms=12, marker='o')
        #plt.axvline(fr, color=p[0].get_color(), ls='--')
        plt.figure(2)
        p = plt.plot(capY['frequency'], np.imag(Yin), ls='solid', ms=12, marker='o')
        plt.axvline(fr, color=p[0].get_color(), ls='--')
    #plt.figure(1)
    #plt.grid()
    #plt.xlim(left=200, right=400)
    #plt.xlabel('Frequency [MHz]')
    #plt.ylabel('|S21|')

    plt.figure(2)
    plt.grid()
    plt.ylim(top=0.2, bottom=-0.2)
    plt.xlim(left=200, right=400)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Im Yin [1/Ohm]')
    plt.savefig(plotdir + 'Input_admittance_vs_Frequency.png')
    plt.close()

    Nthick = 6
    Ndielectric = 11
    x = np.arange(Ndielectric)*0.1 + 3.5
    y = np.arange(1, Nthick+1)*0.1
    X, Y = np.meshgrid(x,y)
    Z = frs.reshape((-1, Ndielectric))
    Z2 = Qs.reshape((-1, Ndielectric))

    fr0 = Z[2, 4]
    df = (Z - fr0)
    xppm = df/fr0*1e6

    plt.figure(1)
    for i in range(Nthick):
        plt.plot(x, Z[i], ls='None', marker='o', label='%1.1f [um]'%y[i])
    plt.grid()
    plt.xlabel('Dielectric Constant')
    plt.ylabel('Frequency [MHz]')
    lgd = plt.legend(loc='upper left', title='Thickness',\
            bbox_to_anchor=(1.0, 1.0))
    plt.savefig(plotdir + 'frs_vs_dielectric_const.png',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.figure(3)
    for i in range(Nthick):
        plt.plot(x, Z2[i], ls='None', marker='o', label='%1.1f [um]'%y[i])
    plt.grid()
    plt.xlabel('Dielectric Constant')
    plt.ylabel('Q')
    lgd = plt.legend(loc='upper left', title='Thickness',\
            bbox_to_anchor=(1.0, 1.0))
    plt.savefig(plotdir + 'Qs_vs_dielectric_const.png',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.figure(4)
    for i in range(Nthick):
        plt.plot(x, xppm[i], ls='None', marker='o', label='%1.1f [um]'%y[i])
    plt.grid()
    plt.xlabel('Dielectric Constant')
    plt.ylabel('x (ppm)')
    lgd = plt.legend(loc='upper left', title='Thickness',\
            bbox_to_anchor=(1.0, 1.0))
    plt.savefig(plotdir + 'xppm_vs_dielectric_const.png',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.figure(5)
    for i in range(Nthick):
        plt.plot(x, df[i], ls='None', marker='o', label='%1.1f [um]'%y[i])
    plt.grid()
    plt.xlabel('Dielectric Constant')
    plt.ylabel('Frequency Shift [MHz]')
    lgd = plt.legend(loc='upper left', title='Thickness',\
            bbox_to_anchor=(1.0, 1.0))
    plt.savefig(plotdir + 'df_vs_dielectric_const.png',
            bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()

    exit()

    plt.figure(3)
    plt.plot(ildthicknesses, frs)
    plt.grid()
    plt.xlabel('ILD thickness')
    plt.ylabel('Frequency [MHz]')
    plt.show()


    plt.figure(4)
    plt.plot(ildthicknesses, Qs)
    plt.grid()
    plt.xlabel('ILD thickness')
    plt.ylabel('Qs')
    plt.show()
