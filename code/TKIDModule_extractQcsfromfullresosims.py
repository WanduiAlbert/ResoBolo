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
plotdir = 'TKIDModule_fig/'

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


# Connect matched 50 Ohm loads to ports 3 and 4 of the capacitor
def reduce_Yparameters(capY):
    delta = (Y0 + capY['Y33'])*(Y0 + capY['Y44']) - capY['Y43']*capY['Y34']
    #delta = (Y0 - capY['Y33'])*(Y0 - capY['Y44']) - capY['Y43']*capY['Y34']
    deltap = capY['Y33']*capY['Y44'] - capY['Y43']*capY['Y34']


    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    p11 = capY['Y11']
    p11 += (capY['Y34']*capY['Y41'] - capY['Y31']*(Y0 + capY['Y44']))/delta
    #p11 += (capY['Y34']*capY['Y41'] + capY['Y31']*(Y0 - capY['Y44']))/delta
    #p11 += (capY['Y34']*capY['Y41'] - capY['Y31']*capY['Y44'])/deltap
    p12 = capY['Y12']
    p12 += (capY['Y34']*capY['Y42'] - capY['Y32']*(Y0 + capY['Y44']))/delta
    #p12 += (capY['Y34']*capY['Y42'] + capY['Y32']*(Y0 - capY['Y44']))/delta
    #p12 += (capY['Y34']*capY['Y42'] - capY['Y32']*capY['Y44'])/deltap
    p21 = capY['Y21']
    p21 += (capY['Y43']*capY['Y31'] - capY['Y41']*(Y0 + capY['Y33']))/delta
    #p21 += (capY['Y43']*capY['Y31'] + capY['Y41']*(Y0 - capY['Y33']))/delta
    #p21 += (capY['Y43']*capY['Y31'] - capY['Y41']*capY['Y33'])/deltap
    p22 = capY['Y22']
    p22 += (capY['Y43']*capY['Y32'] - capY['Y42']*(Y0 + capY['Y33']))/delta
    #p22 += (capY['Y43']*capY['Y32'] + capY['Y42']*(Y0 - capY['Y33']))/delta
    #p22 += (capY['Y43']*capY['Y32'] - capY['Y42']*capY['Y33'])/deltap
    tup = list(zip(capY['frequency'], p11, p12, p21, p22))
    return np.array(tup, dtype=dtypes)

def get_resonator_Sparams(capY, indY):
    y11 = capY['Y11'] + indY['Y11']
    y12 = capY['Y12'] + indY['Y12']
    y21 = capY['Y21'] + indY['Y21']
    y22 = capY['Y22'] + indY['Y22']
    #delta = y11*y22 - y12*y21
    delta = y11- y12

    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])

    #p11 = capY['Y33'] - (y22*capY['Y13'] - y12*capY['Y23'])/delta
    #p12 = capY['Y34'] - (y22*capY['Y14'] - y12*capY['Y24'])/delta
    #p21 = capY['Y43'] - (-y21*capY['Y13'] + y11*capY['Y23'])/delta
    #p22 = capY['Y44'] - (-y21*capY['Y14'] + y11*capY['Y24'])/delta
    p11 = capY['Y33'] + capY['Y13']*(capY['Y32'] - capY['Y31'])/delta
    p12 = capY['Y34'] + capY['Y14']*(capY['Y32'] - capY['Y31'])/delta
    p21 = capY['Y43'] + capY['Y13']*(capY['Y42'] - capY['Y41'])/delta
    p22 = capY['Y44'] + capY['Y14']*(capY['Y42'] - capY['Y41'])/delta
    tup = list(zip(capY['frequency'], p11, p12, p21, p22))

    combinedY = np.array(tup, dtype=dtypes)

    return convert_Y_to_S(combinedY, nports=2)


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
    delta = (Y0 + capY['Y33'])*(Y0 + capY['Y44']) - capY['Y43']*capY['Y34']

    p11 = (capY['Y11'] + indY['Y11']) - (capY['Y12'] + indY['Y12'])
    #p11.real = (indY['Y11'] - indY['Y12']).real
    #p11.real = capY['Y12'].real
    p11 += capY['Y13']*(capY['Y34']*(capY['Y41'] - capY['Y42']) - (Y0 +
        capY['Y44'])*(capY['Y31'] - capY['Y32']))/delta
    p11 += capY['Y14']*(capY['Y43']*(capY['Y31'] - capY['Y32']) - (Y0 +
        capY['Y33'])*(capY['Y41'] - capY['Y42']))/delta
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

def convert_Y_to_S(Yparams, nports=2):
    p = 'S'
    if nports == 1:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128)])
        p11 = (Y0 - Yparams['Y11'])/(Y0 + Yparams['Y11'])
        tup = list(zip(Yparams['frequency'], p11))
    elif nports==2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
        delta = (Yparams['Y11'] + Y0)*(Yparams['Y22'] + Y0) - Yparams['Y21']*Yparams['Y12']

        p11 = ((Y0 - Yparams['Y11'])*(Y0 + Yparams['Y22']) +
                Yparams['Y12']*Yparams['Y21'])/delta
        p12 = (-2*Yparams['Y12']*Y0)/delta
        p21 = (-2*Yparams['Y21']*Y0)/delta
        p22 = ((Y0 + Yparams['Y11'])*(Y0 - Yparams['Y22']) +
                Yparams['Y12']*Yparams['Y21'])/delta
        tup = list(zip(Yparams['frequency'], p11, p12, p21, p22))
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

def admittance_model(x, C, L, R):
    w = 2*pi*x
    num = (1 - w**2*L*C)**2 + (w*R*C)**2
    denom = R**2*(1-w**2*L*C)**2
    #return np.log(np.sqrt(num/denom) )
    return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
    #return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

def admittance_model_parallel(x, C, L, R):
    w = 2*pi*x
    L = 0
    return np.log(np.abs(1j*w*C + 1./R) )
    #return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))
def tline_model(x, Zc, w0, eps):
    w = 2*pi*x
    return np.abs(-1/Zc/(eps*np.cos(w/w0) + 1j*np.sin(w/w0)))

def chisq(theta, x, y):
    return np.sum((y - admittance_model(x, *theta))**2)

def get_S21(Y):
    Y0 = 1/Z0
    DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
    return -2 * Y['Y21']*Y0/DY

def get_cap_params(Ydata, savename):
    f = Ydata['frequency'] * MHz
    nY21 = -Ydata['Y21']
    Yp1 = Ydata['Y11'] + Ydata['Y21']
    Yp2 = Ydata['Y22'] + Ydata['Y21']
    Yodd = Ydata['Y11'] - Ydata['Y21']

    wpeak = 2*pi*f[np.argmax(nY21.imag)]
    C_est = nY21.imag[0]/(2*pi*f[0])
    L_est = 1/(wpeak**2*C_est)
    fpeak = wpeak/2/pi
    #L_est = 1e-20
    R_est = 1e-8
    #C1_est = 2.0*pF
    nY21 = nY21[f < fpeak]
    f = f[f < fpeak]

    p0 = [C_est, L_est, R_est]
    #pdb.set_trace()
    popt, pcov = optimize.curve_fit(admittance_model, f, np.log(np.abs(nY21)), p0=p0, method='lm')
    #popt, pcov = curve_fit(tline_model, f, np.log(np.abs(nY21)), p0=p0, method='lm')
    #result = minimize(chisq, p0, args=(f, np.abs(nY21)),
    #        method='Nelder-Mead')
    C_fit, L_fit, R_fit = popt
    #C_fit, L_fit, R_fit, C1_fit = result['x']
    #print (p0)
    #print (popt)
    y_est = np.exp(admittance_model(f, C_est, L_est, R_est))

    y_fit = np.exp(admittance_model(f, C_fit, L_fit, R_fit))
    fig, ax =plt.subplots(num = 345, figsize=(10,10))
    ax.plot(f/MHz, np.abs(nY21), 'b', label='Simulation')
    #ax.plot(f/MHz, y_est, 'g', label='Guess')
    ax.plot(f/MHz, y_fit, 'k--',
                label="Fit C = %1.3fpF L = %1.3fnH R=%1.3e Ohms "%(C_fit/pF,
                    L_fit/nH, R_fit))
    #ax.plot(f/MHz, y_est, 'r-',
        #        label="Guess: C = %1.3fpF L = %1.3fnH R=%1.3e Ohms"%(C_est/pF,
        #            L_est/nH, R_est))
    ax.legend()
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|-Y21| [1/Ohms]')
    ax.grid()
    ax.axis('tight')
    plt.savefig(plotdir + savename + "Y21.png")
    plt.close()

    fig, ax = plt.subplots(num=346, figsize=(10,10))
    ax.scatter(f/MHz, np.abs(nY21) - y_fit,
                label="Fit C = %1.3fpF L = %1.3fnH R=%1.3e Ohms "%(C_fit/pF,
                    L_fit/nH, R_fit))
    ax.legend()
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Fit Residuals [1/Ohms]')
    ax.grid()
    ax.axis('tight')
    plt.savefig(plotdir + savename + "Y21_residuals.png")
    #plt.show()
    plt.close()

    return C_fit, L_fit, R_fit




if __name__=="__main__":
    if not doS21analysis:
        capfilenames = glob.glob(datadir +
                "TKID_Module_Capacitor_with_feedline_*outof16.csv")
        capfilenames.sort(key=lambda x: int(x.split("/")[-1].split('.')[0].split('_')[-1][:-7]))
        indfilename = datadir + "TKID_Module_Inductor.csv"
        indY = load_data_single(indfilename, nports=2, paramtype='Y')
        frs = np.zeros(len(capfilenames))
        Qs = np.zeros(len(capfilenames))
        Cs = np.zeros(len(capfilenames))
        Ls = np.zeros(len(capfilenames))
        Rs = np.zeros(len(capfilenames))
        #numsections = np.array([2, 3, 4, 6, 8, 10, 12, 14, 16])
        numsections = np.array([4, 6, 8, 10, 12, 14, 16])
        #numsections = np.array([4, 6])
        Nfingers = numsections*53
        #labels = ['2', '3', '4', '6', '8', '10', '12', '14', '16']
        labels = ['4', '6', '8', '10', '12', '14', '16']
        #labels = ['4', '6']

        plt.figure(1, figsize=(10,10))
        plt.figure(2, figsize=(10,10))
        for i, capfn in enumerate(capfilenames):
            savename = capfn.split("/")[-1].split('.')[0].split('_')[-2]
            capY = load_data_single(capfn, nports=4, paramtype='Y')
            Ceff = -capY['Y21'].imag/(2*pi*capY['frequency']*MHz)/pF
            resoS = get_resonator_Sparams(capY, indY)
            finalY = get_Yin_parallel(capY, indY)
            capY = reduce_Yparameters(capY)
            C,L,R = get_cap_params(capY, savename)
            Cs[i] = C
            Ls[i] = L
            Rs[i] = R
            #print (C/pF, L/nH, R)

            Yin = finalY['Y11']
            f = capY['frequency']
            mask = f > 200
            #mask[f > 400] = False
            tck = interpolate.splrep(f[mask], Yin[mask].imag, s=0)
            retck = interpolate.splrep(f[mask], Yin[mask].real, s=0)
            fr = interpolate.sproot(tck)
            fr = fr[0]
            f_fine = np.r_[-10:10:20000j] + fr
            #yl = mesh.get_Vreso_mesh(1j*2*pi*f_fine*MHz)
            #s21 = mesh.get_S21_mesh(1j*2*pi*f_fine*MHz)
            slope = interpolate.splev(fr, tck, der=1)
            dR,  = interpolate.splev([fr], retck)
            R = 1./dR
            #print (slope, R)
            Q = slope*fr/(2*dR)
            imyinterp = np.abs(dR*2*Q)*(f_fine - fr)/fr
            #print (R)
            #Q = slope*R*fr/2
            #print (labels[i], fr, Q)
            frs[i] = fr
            Qs[i] = Q
            plt.figure(1)
            p = plt.plot(resoS['frequency'], np.abs(resoS['S21']), ls='solid',
                    ms=12, marker='o', label=labels[i])
            #plt.plot(f_fine, np.abs(s21), 'k')
            plt.axvline(fr, color=p[0].get_color(), ls='--')
            plt.figure(2)
            p = plt.plot(capY['frequency'], np.imag(Yin), ls='solid',
                    ms=12, marker='o', label=labels[i])
            #plt.plot(f_fine, imyinterp, 'k')
            #plt.plot(f_fine, np.imag(yl), 'k')
            plt.axvline(fr, color=p[0].get_color(), ls='--')
            plt.figure(3)
            p = plt.plot(capY['frequency'], np.real(Yin), ls='solid',
                    ms=12, marker='o', label=labels[i])
            #plt.plot(f_fine, np.real(yl), 'k')
            plt.axvline(fr, color=p[0].get_color(), ls='--')
            plt.figure(4)
            p = plt.plot(capY['frequency'], Ceff, ls='solid', label=labels[i])
            #plt.plot(f_fine, np.real(yl), 'k')

        print (frs)
        print (Qs)
        print ("Capacitance ", Cs/pF)
        print ("Parasitic inductance ", Ls/nH)
        print (Rs)
        exit()

        plt.figure(1)
        plt.grid()
        plt.legend(loc='upper left')
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('|S21|')
        plt.savefig(plotdir + 'S21_vs_Frequency.png')

        plt.figure(2)
        plt.grid()
        plt.legend(loc='upper left')
        #plt.ylim(top=0.2, bottom=-0.2)
        #plt.xlim(left=200, right=400)
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Im Yin [1/Ohm]')
        plt.savefig(plotdir + 'Input_admittance_vs_Frequency.png')

        plt.figure(3)
        plt.grid()
        plt.legend(loc='upper left')
        #plt.ylim(top=0.2, bottom=-0.2)
        #plt.xlim(left=200, right=400)
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Re Yin [1/Ohm]')
        plt.savefig(plotdir + 'Input_realadmittance_vs_Frequency.png')

        plt.figure(4)
        plt.grid()
        plt.legend(loc='upper left')
        plt.ylim(top=25, bottom=0)
        plt.xlim(left=100, right=600)
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Ceffective [pF]')
        plt.savefig(plotdir + 'ceff_vs_Frequency.png')
        #plt.show()
        plt.close('all')
        #exit()

        #f0 = 600
        #def fr_vs_nfingers(N, N0, a, b):
        #    return np.log(f0) + a*np.log(N) + b*np.log(N0-N)

        N0 = 600
        def fr_vs_nfingers2(N, fa, alpha, beta, gamma, delta, epsilon):
            x = N/N0
            fmain = fa/np.sqrt(x)
            return fmain/np.sqrt(1 + alpha + beta*x + gamma*x**2 + delta*x**3 +
                    epsilon/x)

        def fr_vs_nfingers(N, fa, alpha, beta, epsilon):
            x = N/N0
            fmain = fa/np.sqrt(x)
            return fmain/np.sqrt(1 + alpha + beta*x + epsilon/x)

        def fr_wrapper(theta, N):
            fa, alpha, beta, epsilon = theta
            return fr_vs_nfingers(N, fa, alpha, beta, epsilon)

        mask = np.ones_like(Nfingers, dtype=bool)
        #mask[:2] = False
        p0 = [403, 0.007, 0.2538, 0.1144]
        #p02 = [403, 0.007, 0.2538, 0.0065, 0.0098, 0.1144]
        bounds = ([0,0,0,0], [np.inf]*4)
        #bounds2 = ([0,0,0,0,0,0], [np.inf]*6)
        popt, pcov = optimize.curve_fit(fr_vs_nfingers, Nfingers[mask],
                frs[mask], absolute_sigma=False, p0=p0, bounds=bounds, method='trf')
        sigma = np.sqrt(np.diag(pcov))
        epsilon = sigma
        #popt2, pcov2 = optimize.curve_fit(fr_vs_nfingers2, Nfingers[mask], frs[mask],
        #        p0=p02, bounds=bounds2)
        #sigma2 = np.sqrt(np.diag(pcov2)
        print ("\n\n")
        print ("The best fit parameters are:", popt)
        print ("The errors on the best fit parameters are:", sigma)
        #print ("The best fit parameters are:", popt2)
        #print ("The errors on the best fit parameters are:", sigma2)
        print ("\n\n")
        Nfine = np.arange(50, 1000, 1)
        frfine = fr_vs_nfingers(Nfine, *popt)
        grad = np.array(list(map(lambda x: optimize.approx_fprime(popt,
            fr_wrapper, epsilon, x), Nfine)))
        df = np.sqrt(np.diag(np.dot(grad, np.dot(pcov, grad.T))))/MHz
        #frfine2 = fr_vs_nfingers2(Nfine, *popt2)
        L = 10*nH
        w0s = 2*pi*frs*MHz
        Cs = 1./(w0s**2*L)
        Ca = 2.1355*pF
        Cb = 1.3889*pF
        Cc = 0.0751*pF
        Cg = Ca - Cc
        Ct = Ca + Cb + Cc
        Qcs = 2*Cs/(w0s*Z0*Cc**2)*(Ct/Ca)**2
        Qct = 20000
        c = np.sqrt(w0s*Z0*Qct/(2*Cs))
        a = Cb + Cg
        b = Cg
        Ccs = (2 - b*c)/(2*c) + np.sqrt(4 + 4*(a-b)*c + b**2*c**2)/(2*c)
        #print (Nfingers)
        #print (Ccs/pF)
        #exit()
        residuals = (frs - fr_vs_nfingers(Nfingers, *popt))
        chisq = np.sum(residuals**2/(frs.size - popt.size))
        print ("Chi squared of the best fit is ", chisq)
        print ("Mean of the residuals: ", np.mean(residuals))
        print ("Std. dev of the residuals: ", np.std(residuals))
        #residuals2 = (frs - fr_vs_nfingers2(Nfingers, *popt2))

        plt.figure()
        plt.plot(Nfingers, frs, 'ko', ls='None', ms=12)
        plt.plot(Nfine, frfine, 'r-')
        plt.fill_between(Nfine, frfine + df, frfine - df, color='gray', alpha=0.2)
        #plt.plot(Nfine, frfine2, 'b-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('fr [MHz]')
        plt.savefig('frequency_vs_Nfingers.png')
        #plt.close()
        #plt.show()

        plt.figure()
        plt.plot(Nfingers, residuals, 'ko', ls='None', ms=12)
        #plt.plot(Nfingers, residuals2, 'rs', ls='None', ms=12)
        plt.fill_between(Nfine, df, -df, color='gray', alpha=0.2)
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('Residuals [MHz]')
        plt.savefig('frequencyresiduals_vs_Nfingers.png')
        #plt.close()
        plt.show()

        plt.figure()
        plt.semilogy(Nfingers, np.abs(Qs), 'ko', ls='None', ms=12)
        #plt.semilogy(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('Qc')
        plt.savefig('extractedQc_vs_Nfingers.png')
        plt.close()
        #plt.show()

        plt.figure()
        plt.semilogy(frs, np.abs(Qs), 'ko', ls='None', ms=12)
        #plt.semilogy(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('Qc')
        plt.savefig('extractedQc_vs_frequency.png')
        plt.close()
        #plt.show()

        plt.figure()
        plt.plot(frs, Cs/pF, 'ko', ls='None', ms=12)
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('C [pF]')
        plt.savefig('Cs_vs_frequency.png')

        p = np.polyfit(Nfingers, Cs/pF, 2)
        print (p)
        Cfine = np.polyval(p, Nfine)
        plt.figure()
        plt.plot(Nfingers, Cs/pF, 'ko', ls='None', ms=12)
        plt.plot(Nfine, Cfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('C [pF]')
        plt.savefig('Cs_vs_Nfingers.png')

        plt.figure()
        plt.plot(frs, Ls/nH, 'ko', ls='None', ms=12)
        #plt.semilogy(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('Lpar [nH]')
        plt.savefig('Lpar_vs_frequency.png')

        p = np.polyfit(Nfingers, Ls/nH, 2)
        print (p)
        Lfine = np.polyval(p, Nfine)
        plt.figure()
        plt.plot(Nfingers, Ls/nH, 'ko', ls='None', ms=12)
        plt.plot(Nfine, Lfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('Lpar [nH]')
        plt.savefig('Lpar_vs_Nfingers.png')

        plt.figure()
        plt.plot(frs, Rs, 'ko', ls='None', ms=12)
        #plt.semilogy(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('Rs [Ohms]')
        plt.savefig('loss_vs_frequency.png')

        plt.figure()
        plt.plot(Nfingers, Rs, 'ko', ls='None', ms=12)
        #plt.semilogy(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('Rs [Ohms]')
        plt.savefig('loss_vs_Nfingers.png')
        #plt.show()

    else:
        fns = glob.glob(datadir + "TKID_Module_FullResonator_*outof16.csv")
        fns.sort(key=lambda x: int(x.split("/")[-1].split('.')[0].split('_')[-1][:-7]))
        nsections = list(map(lambda x: int(x.split("/")[-1].split('.')[0].split('_')[-1][:-7]), fns))
        print (nsections)
        frs = np.zeros(len(nsections))
        Qcs = np.zeros(len(nsections))
        phics = np.zeros(len(nsections))
        Nfingers = np.array(nsections)*53
        for ires,fn in enumerate(fns):
            section = nsections[ires]
            Sparams = load_data_single(fn, nports=2, paramtype='S')
            f = Sparams['frequency']
            re = Sparams['S21'].real
            im = Sparams['S21'].imag
            funcre = interpolate.interp1d(f, re)
            funcim = interpolate.interp1d(f, im)

            finterp = np.r_[f[0]:f[-1]:5000j]
            reinterp = funcre(finterp)
            iminterp = funcim(finterp)
            maginterp = np.sqrt(reinterp*reinterp + iminterp*iminterp)
            dBmaginterp = 10*np.log10(maginterp)

            #f = finterp
            #re = reinterp
            #im = iminterp
            mag = np.sqrt(re*re + im*im)
            #delay = 0.0e-9
            #phase_corr = 1e6*f*delay*2.*pi
            #phi = np.arctan2(im,re)+phase_corr
            #phi = np.unwrap(phi,axis=-1)
            #z = np.sqrt(mag)*np.exp(1.j*phi)
            #re = z.real
            #im = z.imag
            z = re + 1j*im
            dbmag = 10*np.log10(mag)

            #mask = f < 331.2
            mask = np.ones_like(f, dtype=bool)
            midpt = f.size //2  + 1
            #mask[89:90] = False
            try:
                f0,A,m,phi,D,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov =\
                reso_fit.do_fit(f[mask]*1e6,re[mask],im[mask], get_cov=True)
                sigmas = np.sqrt(np.diag(pcov))
            except RuntimeError:
                print ("Error obtaining the fit")
                #f0 = frequencies[ires]
                continue
            Qe = Qe_re + 1j * Qe_im
            dQe = 1/Qe
            phi_c = np.arctan2(Qe_im, Qe_re)
            Qc = 1./np.real(dQe)
            popt = [f0,A,m,phi,D,1./Qr,dQe.real,dQe.imag,a]
            print("Best fit results: ", section, f0, Qi, Qc, phi_c)
            frs[ires] = f0
            Qcs[ires] = Qc
            phics[ires] = phi_c
            res_re = mre - re[mask]
            res_im = mim - im[mask]
            res_mag = res_re*res_re + res_im*res_im
            resdbmag = 10*np.log10(res_mag)
            mmag = mre*mre + mim*mim
            mdbmag = 10*np.log10(mmag)
            nmask = mask

            if plot_diagnostic:
                plt.figure(ires)
                plt.plot(f[mask],mdbmag,'r', label='fit')
                plt.plot(f[mask],dbmag[mask], 'b', label='Qr=%d Qi=%d Qc=%d phi_c=%1.3f'%(Qr,Qi,Qc, phi_c))
                #plt.plot(finterp, dBmaginterp, 'k')
                #plt.ylim(top=max(dbmag)+0.5, bottom=min(dbmag)-12.5)
                plt.xlim(left=f0-1.5, right=f0+1.5)
                plt.grid()
                plt.xlabel('Frequency (MHz)')
                plt.ylabel('|S21|')
                lgd = plt.legend(loc="upper left", bbox_to_anchor=(1.0,1.0))
                path = os.path.join(plotdir,'reso_%d'%section + ".png")
                #plt.savefig(path)
                plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close()

                plt.figure(300 + ires )
                scale = np.max(mmag)**0.5
                plt.plot(re/scale,im/scale, 'b', label='data')
                plt.plot(mre/scale,mim/scale,'r', label='fit')
                plt.grid()
                plt.axis('square')
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.legend(loc='upper right')
                path = os.path.join(plotdir,'reso_IQ_%d'%section + ".png")
                plt.savefig(path)
                #plt.show()
                plt.close()

        f0 = 500
        N0 = 600
        p = np.polyfit(Nfingers/N0, Qcs, 2)
        print (p)
        Nfine = np.arange(50, 1000, 1)
        Qcfine = np.polyval(p, Nfine/N0)
        Qcresiduals = (Qcs - np.polyval(p, Nfingers/N0))/Qcs

        plt.figure(123)
        plt.plot(Nfingers, Qcs, 'ko', ls='None', ms=12)
        plt.plot(Nfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('Qc ')
        plt.savefig('Qc_vs_Nfingers.png')

        p = np.polyfit(Nfingers/N0, phics, 2)
        print (p)
        phicfine = np.polyval(p, Nfine/N0)
        phicresiduals = (phics - np.polyval(p, Nfingers/N0))/phics
        plt.figure(124)
        plt.plot(Nfingers, phics, 'ko', ls='None', ms=12)
        plt.plot(Nfine, phicfine, 'r-')
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('phi_c [radians]')
        plt.savefig('phic_vs_Nfingers.png')

        plt.figure(125)
        plt.plot(Nfingers, Qcresiduals, 'ko', ls='None', ms=12)
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('fractional Qc residuals')
        plt.savefig('Qcresiduals_vs_Nfingers.png')

        plt.figure(126)
        plt.plot(Nfingers, phicresiduals, 'ko', ls='None', ms=12)
        plt.grid()
        plt.xlabel('Nfingers')
        plt.ylabel('fractional phic residuals')
        plt.savefig('phicresiduals_vs_Nfingers.png')
        #plt.show()
        plt.close('all')

        p = np.polyfit(frs/f0, Qcs, 2)
        print (p)
        frfine = np.r_[200:800:2000j]
        Qcfine = np.polyval(p, frfine/f0)
        Qcresiduals = (Qcs - np.polyval(p, frs/f0))/Qcs

        plt.figure(123)
        plt.plot(frs, Qcs, 'ko', ls='None', ms=12)
        plt.plot(frfine, Qcfine, 'r-')
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('Qc ')
        plt.savefig('Qc_vs_fr.png')

        p = np.polyfit(frs/f0, phics, 2)
        print (p)
        phicfine = np.polyval(p, frfine/f0)
        phicresiduals = (phics - np.polyval(p, frs/f0))/phics
        plt.figure(124)
        plt.plot(frs, phics, 'ko', ls='None', ms=12)
        plt.plot(frfine, phicfine, 'r-')
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('phi_c [radians]')
        plt.savefig('phic_vs_fr.png')

        plt.figure(125)
        plt.plot(frs, Qcresiduals, 'ko', ls='None', ms=12)
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('fractional Qc residuals')
        plt.savefig('Qcresiduals_vs_fr.png')

        plt.figure(126)
        plt.plot(frs, phicresiduals, 'ko', ls='None', ms=12)
        plt.grid()
        plt.xlabel('fr [MHz]')
        plt.ylabel('fractional phic residuals')
        plt.savefig('phicresiduals_vs_fr.png')

        #plt.show()

