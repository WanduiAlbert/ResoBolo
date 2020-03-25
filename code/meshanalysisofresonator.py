#! /usr/bin/env/python3

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, epsilon_0, mu_0
from math import pi
import reso_fit

MHz = 1e6
nH = 1e-9
pF = 1e-12
um = 1e-6

Z0 = 50
Y0 = 1./Z0

plotdir = 'extractQc_fig/'

# Circuit parameters. Obtained from pi model calculations of the circuit
# parameters between 310 and 340 MHz for the 4 port capacitor sim and the 2 port
# inductor sim.
# Relevant files:
# ../numerical_sims/Cap_300MHz_with_boundary_and_coupcaps_variableILDthickness_Yparams.csv
# ../numerical_sims/inductor_with_GPboundary_variableILDthickness_Yparams.csv

l = 3000*um # length of the feedline

# Resonator parameters
C = 32.1534*pF      # Capacitance between port 1 and 2
L = 9.8529*nH       # Inductance between port 1 and 2
Cp = 0.02636*pF     # Parasitic capacitance between 1 and 2 for the inductor
Lp = 25.7523*nH     # Parasitic inductance between 1 and 2 for the capacitor
Rc = 1e-6           # Loss in the capacitor
Ri = 1e-6           # Loss in the inductor

#Coupling to feedline
C3 = 0.4050*pF      # Capacitance to port 4 from port 1. Nominally Cc
C5 = 0.1612*pF      # Capacitance to port 3 from port 2. Nominally absent

#Coupling to GND
C1 = 0.2118*pF      # Capacitance from port 1 to GND. Nominaly absent.
C4 = 0.6693*pF      # Capacitance from port 2 to GND. Nominally Cc

#Feedline parameters
Lf = 1.3100*nH      # Through inductance of the feedline
C6 = 0.4811*pF      # Capacitance to GND from port 3.
C7 = 0.4759*pF      # Capacitance to GND from port 4.

def get_resonator_impedance(s):
    ycap = s*C + 1./(Rc + s*Lp)
    yind = s*Cp + 1./(Ri + s*L)
    return 1./(ycap + yind)

def get_meshmatrix(s):
    N = s.size
    meshmat = np.zeros((N,6,6), dtype=np.complex128)
    zreso = get_resonator_impedance(s)
    meshmat[:, 0,0] = Z0 + 1./(s*C6)
    meshmat[:, 1,1] = 1./(s*C6) + 1./(s*C5) + 1./(s*C1)
    meshmat[:, 2,2] = s*Lf + 1./(s*C5) + 1./(s*C3) + zreso
    meshmat[:, 3,3] = 1./(s*C1) + zreso + 1./(s*C4)
    meshmat[:, 4,4] = 1./(s*C3) + 1./(s*C4) + 1./(s*C7)
    meshmat[:, 5,5] = Z0 + 1./(s*C7)

    meshmat[:, 0,1] = -1./(s*C6)
    meshmat[:, 1,2] = -1./(s*C5)
    meshmat[:, 1,3] = -1./(s*C1)
    meshmat[:, 2,3] = -zreso
    meshmat[:, 2,4] = -1./(s*C3)
    meshmat[:, 3,4] = -1./(s*C4)
    meshmat[:, 4,5] = -1./(s*C7)

    meshmat[:, 1,0] = meshmat[:, 0,1]
    meshmat[:, 2,1] = meshmat[:, 1,2]
    meshmat[:, 3,1] = meshmat[:, 1,3]
    meshmat[:, 3,2] = meshmat[:, 2,3]
    meshmat[:, 4,2] = meshmat[:, 2,4]
    meshmat[:, 4,3] = meshmat[:, 3,4]
    meshmat[:, 5,4] = meshmat[:, 4,5]

    return meshmat


if __name__=="__main__":
    f = np.r_[300:350:2000j]*MHz
    #f = np.r_[330:332:8000j]*MHz
    omega = 2*pi*f
    s = 1j*omega

    zreso = get_resonator_impedance(s)
    meshmatrix = get_meshmatrix(s)
    voltages = np.array([1,0,0,0,0,0])
    invmeshmatrix = np.linalg.inv(meshmatrix)
    meshcurr = np.dot(invmeshmatrix, voltages)
    Vreso = (meshcurr[:, 3] - meshcurr[:, 4])*zreso
    S21 = 2*meshcurr[:, -1]*Z0

    plt.figure()
    plt.plot(f/MHz, np.abs(S21))
    plt.grid()
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('|S21|')
    plt.show()

    re = S21.real
    im = S21.imag
    mag = np.sqrt(re*re + im*im)
    delay = 0.0e-9
    phase_corr = 1e6*f*delay*2.*pi
    phi = np.arctan2(im,re)+phase_corr
    phi = np.unwrap(phi,axis=-1)
    z = np.sqrt(mag)*np.exp(1.j*phi)
    re = z.real
    im = z.imag
    dbmag = 10*np.log10(mag)

    mask = np.ones_like(f, dtype=bool)
    try:
        f0,A,m,phi,D,Qi,Qr,Qe_re,Qe_im,a,mre,mim,pcov =\
        reso_fit.do_fit(f[mask],re[mask],im[mask], get_cov=True)
        sigmas = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print ("Error obtaining the fit")
        exit()
    Qe = Qe_re + 1j * Qe_im
    dQe = 1/Qe
    phi_c = np.arctan2(Qe_im, Qe_re)
    Qc = 1./np.real(dQe)
    popt = [f0,A,m,phi,D,1./Qr,dQe.real,dQe.imag,a]
    print(f0, Qi, Qc)
    res_re = mre - re[mask]
    res_im = mim - im[mask]
    res_mag = res_re*res_re + res_im*res_im
    resdbmag = 10*np.log10(res_mag)
    mmag = mre*mre + mim*mim
    mdbmag = 10*np.log10(mmag)
    freq = 334.0
    nmask = mask
    plot_diagnostic=True
    ires = 0

    if plot_diagnostic:
        plt.figure(ires )
        plt.plot(f[mask]/MHz,mdbmag,'r', label='fit')
        plt.plot(f[mask]/MHz,dbmag[mask], 'b', label='Qr=%d Qi=%d Qc=%d phi_c=%1.3f'%(Qr,Qi,Qc, phi_c))
        #plt.plot(finterp, dBmaginterp, 'k')
        #plt.ylim(top=max(dbmag)+0.5, bottom=min(dbmag)-12.5)
        plt.grid()
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('|S21|')
        lgd = plt.legend(loc="upper left", bbox_to_anchor=(1.0,1.0))
        path = os.path.join(plotdir,'reso_%d'%freq + "MHz.png")
        #plt.savefig(path)
        plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()
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
        path = os.path.join(plotdir,'reso_IQ_%d'%freq + "MHz.png")
        plt.savefig(path)
        plt.show()
        plt.close()





