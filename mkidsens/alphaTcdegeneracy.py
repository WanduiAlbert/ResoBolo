import numpy as np
from scipy.constants import h,k,c,e
import scipy.special as special
import matplotlib.pyplot as plt
from matplotlib import cm
from math import pi

K0 = lambda x: special.kn(0, x)
I0 = lambda x: special.iv(0, x)

um = 1e-6
nm = 1e-9
ppm = 1e6
mK = 1e-3

Fd0 = 5.97e-5
Tc0 = 1.284
Delta0 = 1.763*k*Tc0
alphak0 = 0.375
N0 = 1.72e10/um**3/e
R = 10.7*um**3#/s
nqp_star = 518/um**3
Vsc = 16200*um**3*50*nm

def model_xMB(T, f0, Tc, alphak):
    eta = h*f0*1e6/(2*k*T)
    Delta = 1.763*k*Tc
    nqp = 2*N0*np.sqrt(2*Delta*pi*k*T)*np.exp(-Delta/k/T)
    S2 = 1 + np.sqrt(2*Delta/(pi*k*T))*np.exp(-eta)*I0(eta)
    print (np.mean(np.exp(-eta)*I0(eta)))
    xmb = -alphak*S2*nqp/(4*N0*Delta)
    return xmb

def model_dQMB(T, f0, Tc, alphak):
    eta = h*f0*1e6/(2*k*T)
    Delta = 1.763*k*Tc
    nqp = 2*N0*np.sqrt(2*Delta*pi*k*T)*np.exp(-Delta/k/T)
    S1 = (2/pi)*np.sqrt(2*Delta/(pi*k*T))*np.sinh(eta)*K0(eta)
    dQmb = alphak*S1*nqp/(2*N0*Delta)
    return dQmb

def dQdegeneracy_line(T, f0, Tc, alpha1, Tc1):
    eta = h*f0*1e6/(2*k*T)
    Delta = 1.763*k*Tc
    Delta1 = 1.763*k*Tc1
    C = Delta/(k*T)
    C1 = Delta1/(k*T)
    lna = C - C1 + np.log(alpha1)
    return lna

def xdegeneracy_line(T, f0, Tc, alpha1, Tc1):
    eta = h*f0*1e6/(2*k*T)
    Delta = 1.763*k*Tc
    Delta1 = 1.763*k*Tc1
    A = np.sqrt(2/(pi*k*T))*np.exp(-eta)*I0(eta)
    C = np.log(Delta**0.5/(1 + A*Delta**0.5)) + Delta/(k*T)
    C1 = np.log(Delta1**0.5/(1 + A*Delta1**0.5)) + Delta1/(k*T)
    lna = C - C1 + np.log(alpha1)
    return lna

#for Tbase in T:
#    #Tbase = 250e-3
#    x_fit = -model_xMB(Tbase, f0, X, Y)*ppm
#    dQ_fit = model_dQMB(Tbase, f0, X, Y)*ppm
#    beta_fit = dQ_fit/x_fit
#    lna = xdegeneracy_line(Tbase, f0, Tcs, 0.08, 1.1)
#    alpha_degen = np.exp(lna)
#    mask = (alpha_degen < 1)
#    lna = dQdegeneracy_line(Tbase, f0, Tcs, 0.08, 1.1)
#    dQalpha_degen = np.exp(lna)
#    dQmask = (dQalpha_degen < 1)
#
#    fig, ax = plt.subplots(figsize=(10,10))
#    cnt = ax.contour(X, Y, x_fit, 15, colors='black')
#    ax.clabel(cnt, inline=True, fontsize=10)
#    #mp = ax.pcolor(X, Y, x_fit, cmap=cm.viridis, alpha=0.6)
#    mp = ax.imshow(x_fit, extent=[X.min(), X.max(), Y.min(), Y.max()],
#            origin='lower', aspect='auto', cmap=cm.viridis, alpha = 0.6)
#    ax.plot(Tcs[mask], alpha_degen[mask], 'r-')
#    ax.set_xlabel('Tc [K]')
#    ax.set_ylabel('alphak')
#    ax.set_title('-x[ppm] at T = %1.1f mK'%(Tbase/mK))
#    fig.colorbar(mp, label='-x [ppm]')
#    plt.show()
#
#    fig, ax = plt.subplots(figsize=(10,10))
#    cnt2 = ax.contour(X, Y, dQ_fit, 15, colors='black')
#    ax.clabel(cnt2, inline=True, fontsize=10)
#    ##mp2 = ax.pcolor(X, Y, dQ_fit, cmap=cm.viridis, alpha=0.6)
#    mp2 = ax.imshow(dQ_fit, extent=[X.min(), X.max(), Y.min(), Y.max()],
#            origin='lower', aspect='auto', cmap=cm.viridis, alpha = 0.6)
#    ax.plot(Tcs[dQmask], dQalpha_degen[dQmask], 'r-')
#    ax.set_xlabel('Tc [K]')
#    ax.set_ylabel('alphak')
#    ax.set_title('dQ[ppm] at T = %1.1f mK'%(Tbase/mK))
#    fig.colorbar(mp2, label='dQ [ppm]')
#    plt.show()


T = np.r_[80:600:10j]*1e-3
Tcs = np.r_[1.283:1.287:1000j]
alphaks = np.r_[0.375:0.384:1000j] + 1e-3
X,Y = np.meshgrid(Tcs, alphaks)

f0 = 300 #MHz
# Simulate some data
Tcsim = 1.284
alphasim = 0.38
Tsim = np.r_[90:500:10j]*1e-3
trials = []
ysim = model_xMB(Tsim, f0, Tcsim, alphasim)
for i in range(10):
    ysimnoise = 1e-5*np.random.randn(ysim.size)
    trials.append(ysim + ysimnoise)

trials = np.array(trials)
yavg = np.mean(trials, axis=0)
ysig = np.std(trials, axis=0)

print (yavg)
print (ysig)

#plt.figure(figsize=(10,10))
#plt.errorbar(Tsim/1e-3, yavg/1e-6, yerr=ysig/1e-6, fmt='ko')
#plt.show()
#exit()
normalized_err = (yavg[:, np.newaxis, np.newaxis] - model_xMB(Tsim[:, np.newaxis,
    np.newaxis], f0, X, Y))/ysig[:, np.newaxis, np.newaxis]
print (normalized_err.shape)
chisqs = np.log(np.sum(normalized_err**2, axis=0))
print (chisqs.shape)


fig, ax = plt.subplots(figsize=(10,10))
cnt = ax.contour(X, Y, chisqs, 15, colors='black')
ax.clabel(cnt, inline=True, fontsize=10)
#mp = ax.pcolor(X, Y, x_fit, cmap=cm.viridis, alpha=0.6)
mp = ax.imshow(chisqs, extent=[X.min(), X.max(), Y.min(), Y.max()],
        origin='lower', aspect='auto', cmap=cm.viridis, alpha = 0.9)
ax.axhline(alphasim, 0, 1, color='r')
ax.axvline(Tcsim, 0, 1, color='r')
ax.set_xlabel('Tc [K]')
ax.set_ylabel('alphak')
fig.colorbar(mp, label='Log [chisq]')
plt.show()
    #fig, ax = plt.subplots(figsize=(10,10))
    #cnt3 = ax.contour(X, Y, beta_fit, 15, colors='black')
    #ax.clabel(cnt3, inline=True, fontsize=10)
    ##mp3 = ax.pcolor(X, Y, beta_fit, cmap=cm.viridis, alpha=0.6)
    #mp3 = ax.imshow(beta_fit, extent=[X.min(), X.max(), Y.min(), Y.max()],
    #        origin='lower', aspect='auto', cmap=cm.viridis, alpha = 0.6)
    #ax.set_xlabel('Tc [K]')
    #ax.set_ylabel('alphak')
    #ax.set_title('T = %1.1f mK'%(Tbase/mK))
    #fig.colorbar(mp3, label='beta [ppm]')
    #plt.show()
