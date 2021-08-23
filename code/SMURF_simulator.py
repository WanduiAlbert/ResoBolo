#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import reso_fit

MHz = 1e6
pi = np.pi

f0 = 536.8721433585404
A = 0.02926699022089219
m = 0.00026318394797172613
phi = 0.15724347592765547
D = 0.000707779518986162
Qi = 49327.01319462302
Qr = 14027.67694585141
Qe_re = 19475.00575110928
Qe_im = -1573.6654131915298
a = 0.06717203150925918

Qe = Qe_re + 1j*Qe_im
dQe = 1/Qe
Qc = 1./np.real(dQe)
dQc = 1/Qc
dQr = 1./Qr

f = (f0 + np.r_[-1:1:5000j])*MHz

phics = np.r_[-pi:pi:10j]

dfr = 100e3/MHz
f_t = 2

t = np.r_[0:10:1000j]

dfr_t = dfr*np.sin(2*pi*f_t*t)

plt.figure(1, figsize=(10,10))
#plt.figure(2, figsize=(10,10))
for phic in phics:
    phic = 0

    dQe = dQc*(1 + 1j*np.tan(phic))
    s21 = reso_fit.complex_of_real(reso_fit.model_linear_wslope(
        f,f0,A,m,phi,D,dQr,dQe.real,dQe.imag,a))

    Ap = A*np.exp(2.j*pi*(1e-6*D*(f-f0*MHz)+phi))

    s21 /= Ap

    mag = np.abs(s21)
    min_idx = np.argmin(mag)

    df = 10e3
    plus_idx = np.argmin(np.abs(f-f[min_idx]-df))
    minus_idx = np.argmin(np.abs(f-f[min_idx]+df))
    s21_plus = s21[plus_idx]
    s21_minus = s21[minus_idx]
    eta = 2*df/(s21_plus - s21_minus)
    eta_mag = np.abs(eta)
    eta_norm = eta/eta_mag
    eta_phase = np.angle(eta)

    s21_prime = eta_norm*s21

    plt.figure(1)
    plt.plot(s21.real, s21.imag, 'b', label='phi_c = {:1.2f}'.format(phic))
    plt.plot(s21_prime.real, s21_prime.imag, 'r',
            label='phi_eta = {0:1.2f}'.format(eta_phase))
    plt.legend(loc='upper right')
    plt.grid()
    plt.axis('square')
    plt.xlabel('I')
    plt.ylabel('Q')
    plt.show()

    #plt.figure(2)
    #plt.plot(f/MHz, mag)
    ##plt.plot(f/MHz, s21.real, 'r', label='I')
    ##plt.plot(f/MHz, s21.imag, 'b', label='Q')
    #plt.axvline(f[plus_idx]/MHz, color='k')
    #plt.axvline(f[minus_idx]/MHz, color='k')
    #plt.axvline(f[min_idx]/MHz, color='k')
    #plt.axvline(f0, color='r')
    #plt.grid()
    ##plt.legend(loc='upper right')
    #plt.xlabel('f [MHz]')
    #plt.ylabel('I,Q')
    #plt.show()

    f0_smurf = f[min_idx]
    Q0 = np.imag(s21[min_idx]*eta)
    Ap0 = A*np.exp(2.j*pi*(1e-6*D*(f0_smurf-f0*MHz)+phi))
    s21_output = reso_fit.complex_of_real(reso_fit.model_linear_wslope(
        f0_smurf,f0+dfr_t,A,m,phi,D,dQr,dQe.real,dQe.imag,a))
    s21_output /= Ap0
    Q_output = np.imag(s21_output*eta)/MHz
    Q_output -= np.mean(Q_output)

    plt.figure(3)
    plt.plot(t, dfr_t, 'b', label='fr(t)')
    plt.plot(t, Q_output, 'r', label='SMURF output')
    plt.grid()
    plt.title('phi_c = {0:1.2f}, phi_eta = {1:1.2f}'.format(phic, eta_phase))
    plt.legend(loc='upper right')
    plt.ylabel('f [MHz]')
    plt.xlabel('t [s]')
    plt.show()


    exit()



