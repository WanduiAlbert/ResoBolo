import numpy as np
import matplotlib.pyplot as pl
import reso_fit
from scipy.interpolate import RectBivariateSpline
pi = np.pi
from scipy import optimize
from scipy.constants import c, epsilon_0

# transmission line
# Nb/Nb, 0.3um ILD, 20um MS, 80um CPW 20um wide slot

pF = 1e-12


Z0 = 50
xline_eps = 13.4 - 0.0081j
xline_n = np.sqrt(xline_eps)
xline_Z0 = 38.
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
    prefactor = 2*(C+Cc)/(Cc**2*Z0*w0)#*(1 + (Cc*Z0*w0/2)**2)
    lambd = Zc/Z0

    theta_in = gamma*l
    theta_out = gamma*(xline_ltot - l)
    theta_tot = theta_in + theta_out
    cos_tot = np.cos(2*np.imag(theta_tot))
    cosh_tot = np.cosh(2*np.real(theta_tot))
    sinh_tot = np.sinh(2*np.real(theta_tot))
    num = -(1-lambd**2)**2*cos_tot + (1 + 6*lambd**2 + lambd**4)*cosh_tot +\
            4*lambd*(1+lambd**2)*sinh_tot
    denom = 2*lambd*(1-lambd**2)*(
            np.cos(2*np.imag(theta_in))*np.cosh(2*np.real(theta_out)) +
            np.cos(2*np.imag(theta_out))*np.cosh(2*np.real(theta_in)))
    denom += 4*lambd*(1 + lambd**2)*np.cosh(2*np.real(theta_tot))
    denom += (1 - lambd**4)*(
            np.cos(2*np.imag(theta_in))*np.sinh(2*np.real(theta_out)) +
            np.cos(2*np.imag(theta_out))*np.sinh(2*np.real(theta_in)))
    denom += (1 + 6*lambd**2 + lambd**4)*np.sinh(2*np.real(theta_tot))
    denom *= lambd
    return prefactor * np.real(num/denom)

def qc_prediction_monotonic(fs,C, Cc, Zc, n):
    w0 = 2*pi*fs
    gamma = 1.0j*2*pi*fs*n/c
    prefactor = 2*(C+Cc)/(Cc**2*Z0*w0)#*(1 + (Cc*Z0*w0/2)**2)
    lambd = Zc/Z0

    l = (fs/400e6 - 1)*xline_ltot

    theta_in = gamma*l
    theta_out = gamma*(xline_ltot - l)
    theta_tot = theta_in + theta_out
    cos_tot = np.cos(2*np.imag(theta_tot))
    cosh_tot = np.cosh(2*np.real(theta_tot))
    sinh_tot = np.sinh(2*np.real(theta_tot))
    num = -(1-lambd**2)**2*cos_tot + (1 + 6*lambd**2 + lambd**4)*cosh_tot +\
            4*lambd*(1+lambd**2)*sinh_tot
    denom = 2*lambd*(1-lambd**2)*(
            np.cos(2*np.imag(theta_in))*np.cosh(2*np.real(theta_out)) +
            np.cos(2*np.imag(theta_out))*np.cosh(2*np.real(theta_in)))
    denom += 4*lambd*(1 + lambd**2)*np.cosh(2*np.real(theta_tot))
    denom += (1 - lambd**4)*(
            np.cos(2*np.imag(theta_in))*np.sinh(2*np.real(theta_out)) +
            np.cos(2*np.imag(theta_out))*np.sinh(2*np.real(theta_in)))
    denom += (1 + 6*lambd**2 + lambd**4)*np.sinh(2*np.real(theta_tot))
    denom *= lambd
    return prefactor * np.real(num/denom)


ncf = 1001
nfrac = 101
cfs = np.linspace(400e6,800e6,ncf)
#Qes = np.zeros(nfrac,dtype=np.complex)

cmap = pl.get_cmap('jet')
all_Qes = []
all_f0s = []

L0 = 10e-9
Qi = 100000
#Cc = 0.1e-12
Qc = 40000.
C0 = 1/((2*pi*cfs)**2*L0)
Cc = np.sqrt(C0/(pi*cfs*Qc*50.))
C0 -= Cc
#Qc = C0/(pi*cfs*Cc**2*50.)

Zcs = np.r_[30:70:15j]
Qc_mean = np.zeros_like(Zcs)
Qc_std = np.zeros_like(Zcs)
Qc_dev = np.zeros_like(Zcs)



pl.figure(1, figsize=(10,10))
# We now have Qe as a function of the position along the feedline nfrac
for i, Zc in enumerate(Zcs):
    Qcs_pred = qc_prediction_monotonic(cfs,C0, Cc, Zc, xline_n)
    Qc_mean[i] = np.mean(Qcs_pred)
    Qc_std[i] = np.std(Qcs_pred)
    Qc_dev[i] = np.max(Qcs_pred) - np.min(Qcs_pred)

    pl.plot(cfs/1e6, Qcs_pred,ls='-', label='Zc = %dOhms'%Zc)

pl.figure(1)
pl.grid()
pl.legend(loc='upper right')
pl.ylabel('Qc')
pl.xlabel('Frequency [MHz]')
pl.show()
#pl.ylim((0, 50000))


pl.figure(figsize=(10,10))
pl.plot(Zcs, Qc_mean, 'bo')
pl.grid()
pl.xlabel('Zc [Ohms]')
pl.ylabel('Avg. Qc')
pl.show()

pl.figure(figsize=(10,10))
pl.plot(Zcs, Qc_std, 'bo')
pl.plot(Zcs, Qc_dev, 'rd')
pl.grid()
pl.xlabel('Zc [Ohms]')
pl.ylabel('Stdev. Qc')
pl.show()



