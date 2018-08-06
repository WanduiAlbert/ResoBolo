#! /usr/bin/env python3

from math import pi
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('/home/wanduialbert/Desktop/research_stuff/resobolo/masks/code')
from capacitors import IDC
import glob
from scipy.constants import epsilon_0

datadir = '/home/wanduialbert/Desktop/research_stuff/resobolo/numerical_sims/'

MHz = 1e6
pF = 1e-12
mm = 1e-3
um = 1e-6

er_si = 11.8
er_G10 = 5
d_si = 0.5*mm
d_G10 = 0.0651*25.4*mm
d = d_si + d_G10
er_eff = er_si * er_G10 * d / (er_si * d_G10 + er_G10 * d_si)

def get_default_colors():
    prop_cycle = plt.rcParams['axes.prop_cycle']
    return prop_cycle.by_key()['color']
colors = get_default_colors()

def load_admittancedata(fn, i):
    delim = ' ' if i ==0 else ','
    dataset = np.loadtxt(fn, skiprows=8, delimiter=delim)
    dtypes = np.dtype([("frequency", np.float64), ("Y11", np.complex128),\
        ("Y12", np.complex128), ("Y21", np.complex128), ("Y22", np.complex128)])
    y11 = dataset[:, 1] + 1j*dataset[:, 2]
    y12 = dataset[:, 3] + 1j*dataset[:, 4]
    y21 = dataset[:, 5] + 1j*dataset[:, 6]
    y22 = dataset[:, 7] + 1j*dataset[:, 8]
    tup = list(zip(dataset[:, 0], y11, y12, y21, y22))
    return np.array(tup, dtype=dtypes)


def get_capacitances(Y):
    w = 2*pi*Y["frequency"]*MHz
    C12 = (-Y['Y12']/w).imag/pF
    C1g = ((Y['Y11'] + Y['Y12'])/w).imag/pF
    C2g = ((Y['Y22'] + Y['Y12'])/w).imag/pF

    return C12, C1g, C2g

def model_capacitance(f, npairs):
    l = 1000
    w = g = 2
    nfingers = npairs*2
    fgap = 2
    cap = IDC(1.0)
    cap.set_dimensions(w,g,l,fgap,nfingers)
    print (cap.width, cap.height)
    C12 = np.ones_like(f)*cap.capacitance()/pF
    C2gnd_under = parallelplatecapacitance(cap, f)/pF
    C2gnd_over = pointchargecapacitance(cap, f)/pF
    return C12, C2gnd_under, C2gnd_over

def set_axproperties(ax, xlabel, ylabel, title, xlim, uselegend=True):
    if uselegend: ax.legend(loc='best')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(which="both")
    ax.set_xlim(100,1000)

def get_capacitancearea(cap):
    w = cap.trace_width
    g = cap.gap_width
    A_contact = 20*cap.height*um**2
    A_fingers  = (cap.nfinger*cap.finger_length*w)*um**2
    return A_contact + A_fingers

def parallelplatecapacitance(cap, f):
    area = get_capacitancearea(cap)
    print (area**0.5/d)
    return  np.ones_like(f)*epsilon_0*er_eff*area/d
    R = (area/pi)**0.5
    y = epsilon_0*er_eff*(np.log(16*pi*R/d)-1)
    return x

def pointchargecapacitance(cap, f):
    return np.ones_like(f)*8*pi*epsilon_0*er_eff*d

if __name__=="__main__":
    datafiles = glob.glob(datadir + '*_fingerpairs.csv')
    datafiles.sort(key=lambda x: int(x.split('/')[-1].split('_')[0]))
    N = len(datafiles)
    npairs = np.zeros(N, dtype=int)
    Ys = []
    for i in range(N):
        print (datafiles[i].split('/')[-1])
        Ys.append(load_admittancedata(datafiles[i], i))
        npairs[i] = int(datafiles[i].split('/')[-1].split('_')[0])

    # Extract the different capacitances as a function of the number of fingers
    fig, ax =   plt.subplots(num=1, figsize=(10,10))
    fig2, ax2 = plt.subplots(num=2, figsize=(10,10))
    for i in range(N):
        C12, C1g, C2g = get_capacitances(Ys[i])
        f = Ys[i]["frequency"]
        model_C12, model_C2gnd_u, model_C2gnd_o = model_capacitance(f, npairs[i])
        #print (model_C12)
        ax.plot(f, C12, label="Npairs={0:d}".format(npairs[i]))
        ax.plot(f, model_C12, color=colors[i], ls="--",\
            label="model")
        ax2.plot(f, C1g, color=colors[i], label="Npairs={0:d}".format(npairs[i]))
        ax2.plot(f, model_C2gnd_u, color=colors[i], ls="--",\
            label="underestimate")
    ax2.plot(f, model_C2gnd_o, color=colors[N], ls="-.",\
        label="overestimate")
    set_axproperties(ax, "Frequency [MHz]", "Capacitance [pF]",\
        "Capacitance btn ports 1 and 2", [100,1000])
    set_axproperties(ax2, "Frequency [MHz]", "Capacitance to GND [pF]",\
        "Capacitance to GND", [100,1000])
    ax2.set_yscale('log')
    plt.figure(1)
    plt.savefig("Cap_12_simulated_vs_modelled.png")
    #plt.close()
    plt.figure(2)
    plt.savefig("Cap_to_GND.png")
    plt.show()
    
