import numpy as np
from scipy.constants import h,k,c,e
import scipy.special as special
import matplotlib.pyplot as plt
from matplotlib import cm
from math import pi
import scipy.optimize as opt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset




K0 = lambda x: special.kn(0, x)
I0 = lambda x: special.iv(0, x)

um = 1e-6
nm = 1e-9
ppm = 1e-6
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
    #print (np.mean(np.exp(-eta)*I0(eta)))
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


# Generate the "model set" of the bath temperature sweep data

f0 = 334 #MHz
# Simulate some data
Tctrue = 1.284
alphatrue = 0.38

Tmin = 80*mK
Tmaxs = np.arange(300, 850, 50)*mK

Npts = 10

Tcfit = np.zeros_like(Tmaxs)
alphafit = np.zeros_like(Tmaxs)
Tcsigma = np.zeros_like(Tmaxs)
alphasigma = np.zeros_like(Tmaxs)
rhofit = np.zeros_like(Tmaxs)

Tcsigmapred = np.zeros_like(Tmaxs)
alphasigmapred = np.zeros_like(Tmaxs)
rhopred = np.zeros_like(Tmaxs)

Nptsfine = 2000

for irun, Tmax in enumerate(Tmaxs):

    print (irun, Tmax)
    dT = (Tmax - Tmin)/Npts
    T = Tmin + dT*np.arange(Npts)
    dTfine = (Tmax - Tmin)/Nptsfine
    Tfine = Tmin + dTfine*np.arange(Nptsfine)

    ## Generate a dataset for fitting
    trials = []
    ysim = model_xMB(T, f0, Tctrue, alphatrue)
    for i in range(10):
        #ysimnoise = 1e-0*np.random.randn(ysim.size)
        #trials.append(ysim*(1 + ysimnoise))
        ysimnoise = 1e-5*np.random.randn(ysim.size)
        trials.append(ysim + ysimnoise)
    trials = np.array(trials)
    yavg = np.mean(trials, axis=0)
    ysig = np.std(trials, axis=0)

    p0 = [1.3, 0.5]
    fitxMB = lambda x, y, z: model_xMB(x, f0, y, z)
    chisq = lambda x,y: np.sum((yavg - fitxMB(T, x, y))**2/ysig**2)
    popt, pcov = opt.curve_fit(fitxMB, T, yavg, sigma=ysig, p0=p0, method='lm')
    psigma = np.sqrt(np.diag(pcov))
    rho = pcov/np.outer(psigma, psigma)
    #print (popt)
    #print (psigma)
    #print (rho[0,1])

    res = yavg - fitxMB(T, *popt)

    yfine = fitxMB(Tfine, *popt)
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, sharex=True, num=121 + irun,
        figsize=(10,12), constrained_layout=True, gridspec_kw={'height_ratios':[3,1], 'hspace':0})
    ax0.set_title(r'x vs bath temperature')

    #plt.plot(T/mK, yavg, 'ko', label='Data')
    ax0.errorbar(T/mK, yavg/ppm, yerr=ysig/ppm, fmt='ko', ms=12, label='Data')
    ax0.plot(Tfine/mK, yfine/ppm, 'r', label='Fit')
    ax0.set_xlabel('Bath Temp [mK]')
    ax0.set_ylabel('x [ppm]')
    ax0.grid()
    ax0.legend(loc='upper right')

    ax1.plot(T/mK, res/ppm, 'ko', ms=12)
    ax1.set_ylim(bottom=-20, top=20)
    ax1.set_xlabel('Bath Temp [mK]')
    ax1.set_ylabel('Residuals [ppm]')
    ax1.grid()

    plt.show()

    #Calculate the fisher matrix
    fisher = np.zeros((2,2))
    delta_x, delta_y = 1e-6, 1e-6
    fisher[0,0] = (chisq(Tctrue + delta_x, alphatrue)
            - 2*chisq(Tctrue, alphatrue)
            + chisq(Tctrue-delta_x, alphatrue))/delta_x**2
    fisher[0,1] = (chisq(Tctrue + delta_x, alphatrue + delta_y)
            - chisq(Tctrue + delta_x, alphatrue - delta_y)
            - chisq(Tctrue - delta_x, alphatrue + delta_y)
            + chisq(Tctrue - delta_x, alphatrue - delta_y))/(4*delta_x*delta_y)
    fisher[1,1] = (chisq(Tctrue, alphatrue + delta_y)
            - 2*chisq(Tctrue, alphatrue)
            + chisq(Tctrue, alphatrue - delta_y))/delta_y**2
    fisher[1,0] = fisher[0,1]

    cov_pr = np.linalg.inv(fisher)
    sig_pr = np.sqrt(np.diag(cov_pr))
    coh_pr = cov_pr/np.outer(sig_pr, sig_pr)
    rho_pr = coh_pr[0,1]
    Tcsigmapred[irun] = sig_pr[0]
    alphasigmapred[irun] = sig_pr[1]
    rhopred[irun] = rho_pr

    Tcfit[irun] = popt[0]
    alphafit[irun] = popt[1]
    Tcsigma[irun] = psigma[0]
    alphasigma[irun] = psigma[1]
    rhofit[irun] = rho[0,1]
    #exit()


plt.figure(figsize=(10,10))
#plt.plot(Tmaxs/mK, Tcfit, 'bo', ms=12)
plt.errorbar(Tmaxs/mK, Tcfit, yerr=Tcsigma, fmt='bo', ms=12)
plt.axhline(Tctrue, color='r')
plt.xlabel('Max Bath Temperature [mK]')
plt.ylabel('$T_c$ fitted [K]')
plt.grid()
plt.savefig('degeneracy_Tc_vs_bathtemp.png')
plt.show()

plt.figure(figsize=(10,10))
#plt.plot(Tmaxs/mK, alphafit, 'bo', ms=12)
plt.errorbar(Tmaxs/mK, alphafit, yerr=alphasigma, fmt='bo', ms=12)
plt.axhline(alphatrue, color='r')
plt.xlabel('Max Bath Temperature [mK]')
plt.ylabel(r'$\alpha_k$ fitted [K]')
plt.grid()
plt.savefig('degeneracy_alphak_vs_bathtemp.png')
plt.show()

plt.figure(figsize=(10,10))
plt.plot(Tmaxs/mK, rhofit, 'bo', ms=12, label='From Curve Fit')
#plt.plot(Tmaxs/mK, rhopred, 'rd', ms=12, label='From Fisher Matrix')
plt.xlabel('Max Bath Temperature [mK]')
plt.ylabel(r'$\rho_{\alpha_k,\ T_c}$ ')
plt.legend(loc='upper right')
plt.grid()
plt.savefig('degeneracy_rho_vs_bathtemp.png')
plt.show()

plt.figure(figsize=(10,10))
plt.plot(Tmaxs/mK, Tcsigma, 'bo', ms=12, label='From Curve Fit')
plt.plot(Tmaxs/mK, Tcsigmapred, 'rd', ms=12, label='From Fisher Matrix')
plt.xlabel('Max Bath Temperature [mK]')
plt.ylabel('$\sigma_{T_c}$ [K]')
plt.legend(loc='upper right')
plt.grid()
plt.show()

plt.figure(figsize=(10,10))
plt.plot(Tmaxs/mK, alphasigma, 'bo', ms=12, label='From Curve Fit')
plt.plot(Tmaxs/mK, alphasigmapred, 'rd', ms=12, label='From Fisher Matrix')
plt.xlabel('Max Bath Temperature [mK]')
plt.ylabel(r'$\sigma_{\alpha_k}$ [K]')
plt.legend(loc='upper right')
plt.grid()
plt.show()





exit()

T = np.r_[80:600:10j]*1e-3
Tcs = np.r_[1.283:1.287:1000j]
alphaks = np.r_[0.375:0.384:1000j] + 1e-3
X,Y = np.meshgrid(Tcs, alphaks)
#for Tbase in T:
#    #Tbase = 250e-3
#    x_fit = -model_xMB(Tbase, f0, X, Y)/ppm
#    dQ_fit = model_dQMB(Tbase, f0, X, Y)/ppm
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

plt.figure(figsize=(10,10))
plt.plot(Tsim, yavg, 'bo', ms=12)
plt.xlabel('Island Temperature [K]')
plt.ylabel('Fractional Frequency Shift')
plt.show()


p0 = [300, Tcsim, alphasim]
popt, pcov = opt.curve_fit(model_xMB, Tsim, yavg, sigma=ysig, p0=p0, method='lm')
print (popt)
psigma = np.sqrt(np.diag(pcov))
rho = pcov/np.outer(psigma, psigma)
print (psigma)
print (rho)
exit()


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
