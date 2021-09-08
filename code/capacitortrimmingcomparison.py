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
plotdir = 'capacitor_trims/'

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
    a = (capY['Y34']*capY['Y41'] - capY['Y31']*(Y0 + capY['Y44']))/delta
    b = (capY['Y34']*capY['Y42'] - capY['Y32']*(Y0 + capY['Y44']))/delta
    c = (capY['Y43']*capY['Y31'] - capY['Y41']*(Y0 + capY['Y33']))/delta
    d = (capY['Y43']*capY['Y32'] - capY['Y42']*(Y0 + capY['Y33']))/delta

    p11 = capY['Y11'] + capY['Y13'] * a + capY['Y14'] * c
    p12 = capY['Y12'] + capY['Y13'] * b + capY['Y14'] * d
    p21 = capY['Y21'] + capY['Y23'] * a + capY['Y24'] * c
    p22 = capY['Y22'] + capY['Y23'] * b + capY['Y24'] * d


    p = 'Y'
    dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
        (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
    p11 = capY['Y11']
    p11 += capY['Y13']*a + capY['Y14']*c
    p12 = capY['Y12']
    p12 += capY['Y13']*b + capY['Y14']*d
    p21 = capY['Y21']
    p21 += capY['Y23']*a + capY['Y24']*c
    p22 = capY['Y22']
    p22 += capY['Y23']*b + capY['Y24']*d

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
    return np.log(np.sqrt(num/denom) )
    #return np.log(w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2) )
    return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

def admittance_model_cap(x, C, R):
    w = 2*pi*x
    return np.log(np.abs(1j*w*C + 1./R) )



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

def get_capacitance_to_gnd(capY, savename):
    f = capY['frequency'] * MHz
    deltap = capY['Y33']*capY['Y44'] - capY['Y43']*capY['Y34']

    a = (capY['Y34']*capY['Y41'] - capY['Y31']*capY['Y44'])/deltap
    b = (capY['Y34']*capY['Y42'] - capY['Y32']*capY['Y44'])/deltap
    c = (capY['Y43']*capY['Y31'] - capY['Y41']*capY['Y33'])/deltap
    d = (capY['Y43']*capY['Y32'] - capY['Y42']*capY['Y33'])/deltap

    p11 = capY['Y11'] + capY['Y13'] * a + capY['Y14'] * c
    p12 = capY['Y12'] + capY['Y13'] * b + capY['Y14'] * d
    p21 = capY['Y21'] + capY['Y23'] * a + capY['Y24'] * c
    p22 = capY['Y22'] + capY['Y23'] * b + capY['Y24'] * d

    Yp1 = p11 + p21
    Yp2 = p22 + p21

    p0_p1 = [np.mean(np.diff(Yp1.imag)/2/pi/np.diff(f)), 1e15]
    p0_p2 = [np.mean(np.diff(Yp2.imag)/2/pi/np.diff(f)), 1e15]

    popt_p1, pcov_p1 = optimize.curve_fit(admittance_model_cap, f,\
            np.log(np.abs(Yp1)), p0=p0_p1, method='lm')
    popt_p2, pcov_p2 = optimize.curve_fit(admittance_model_cap, f,\
            np.log(np.abs(Yp2)), p0=p0_p2, method='lm')

    C_p1, R_p1 = popt_p1
    C_p2, R_p2 = popt_p2

    yp1_guess = np.exp(admittance_model_cap(f, *popt_p1))
    yp2_guess = np.exp(admittance_model_cap(f, *popt_p2))

    yp1_fit = np.exp(admittance_model_cap(f, C_p1, R_p1))
    yp2_fit = np.exp(admittance_model_cap(f, C_p2, R_p2))


    fig, ax =plt.subplots(num = 343, figsize=(10,10))
    ax.plot(f/MHz, np.abs(Yp1), 'b', label='Simulation')
    #ax.plot(f/MHz, yp1_guess, 'k--')#, label="Fit C = %1.3fpF R=%1.3e Ohms "%(C_p1, R_p1))
    ax.plot(f/MHz, yp1_fit, 'k--',
                label="Fit C = %1.3fpF R=%1.3e Ohms "%(C_p1/pF, R_p1))
    ax.legend()
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Yp1| [1/Ohms]')
    ax.grid()
    ax.axis('tight')
    plt.savefig(plotdir + savename + "Yp1.png")
    plt.close()

    fig, ax =plt.subplots(num = 344, figsize=(10,10))
    ax.plot(f/MHz, np.abs(Yp2), 'b', label='Simulation')
    #ax.plot(f/MHz, y_est, 'g', label='Guess')
    ax.plot(f/MHz, yp2_fit, 'k--',
                label="Fit C = %1.3fpF R=%1.3e Ohms "%(C_p1/pF, R_p1))
    ax.legend()
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('|Yp2| [1/Ohms]')
    ax.grid()
    ax.axis('tight')
    plt.savefig(plotdir + savename + "Yp2.png")
    plt.close()

    return C_p1, C_p2


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
    print ("fpeak", fpeak/MHz)
    #L_est = 1e-20
    R_est = 1e8
    #C1_est = 2.0*pF
    p0 = [C_est, L_est, R_est]
    nY21 = nY21[f < fpeak]
    f = f[f < fpeak]


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
    #plt.show()
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
    capfilenames = glob.glob(datadir +
            "full_capacitor_*_trimming.csv")
    capfilenames = [capfilenames[i] for i in [0, 4, 3, 1, 2]]
    #capfilenames.sort(key=lambda x: int(x.split("/")[-1].split('.')[0].split('_')[-1][:-7]))
    #indfilename = datadir + "TKID_Module_2um_Inductor.csv"
    indfilename = datadir + "dark_res_2um_inductor.csv"
    indY = load_data_single(indfilename, nports=2, paramtype='Y')
    frs = np.zeros(len(capfilenames))
    Cs = np.zeros(len(capfilenames))
    Cp1s = np.zeros(len(capfilenames))
    Cp2s = np.zeros(len(capfilenames))
    Ls = np.zeros(len(capfilenames))
    Rs = np.zeros(len(capfilenames))
    labels = ['notrim', 'simple', 'nostub', '4umstub', '16umstub']

    plt.figure(1, figsize=(10,10))
    plt.figure(2, figsize=(10,10))
    for i, capfn in enumerate(capfilenames):
        savename = capfn.split("/")[-1].split('.')[0].split('_')[-2]
        capY = load_data_single(capfn, nports=2, paramtype='Y')
        Ceff = -capY['Y21'].imag/(2*pi*capY['frequency']*MHz)/pF
        C,L,R = get_cap_params(capY, savename)
        Cs[i] = C
        Ls[i] = L
        Rs[i] = R
        #print (C/pF, L/nH, R)

    print (Cs/pF)
    print (Ls/nH)

    Ltot = Ls + 3.2952*nH
    fr = 1./(2*pi*np.sqrt(Ltot*Cs))/MHz

    deltaC = (Cs[1:] - Cs[1])/2/Cs[1]
    deltaL = (Ls[1:] - Ls[1])/2/Ls[1]

    print (deltaC)
    print (deltaL)
