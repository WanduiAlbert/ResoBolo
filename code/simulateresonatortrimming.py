#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


N0 = 600

Nfingers = np.array([470., 480., 462., 454., 408., 402., 370., 364., 326.,
    330., 320., 336., 316., 310., 306., 288., 284., 262., 248., 264., 234.,
    238., 244., 240., 226., 228.])#, 214.])
x = Nfingers/N0

old_fr = np.array([389.5, 384.3, 394.3, 399.8, 429.6, 434.6, 458., 462.9,
    495.9, 491.3, 500.6, 486.8, 505.4, 509.3, 514.3, 532., 536.3, 562.9,
    580.6, 558.1, 596.7, 592.6, 584.5, 588.5, 609.5, 605.2])#, 626.1])

fr_start = np.array([385.92196121, 387.30927425, 395.54482867, 405.4612828,
    439.64402861, 450.77711654, 479.00890532, 483.04826573, 485.93152532,
    494.97399303, 495.31168392, 498.65299833, 503.70840471, 516.57315632,
    519.13175432 , 547.636004, 554.03028855, 564.07021638, 565.99741343,
    569.89766088, 576.00749857, 578.25418412, 579.48527585, 580.21159329,
    587.81623511, 588.38744788])#, 611.54604793])


fr_target = np.array([389.8, 405.3, 420.8, 436.3, 451.8, 467.4, 483.8, 491.6,
    499.3, 507., 514.8, 522.5 , 530.2, 538., 545.7, 553.4, 561.2, 568.9, 574.2,
    579.4, 584.6, 589.9, 595.2, 600.4, 605.6, 611.1])#, 616.3])


fr_posttrim = np.array([386.84502242, 405.28719113, 422.80701906, 439.26758783,
    449.43517644, 465.13909717, 478.96436532, 486.53430194, 496.76567014,
    503.80791844, 512.88838599, 520.91491858, 531.20604609, 534.79983011,
    544.66718964, 546.81734932, 554.91864219, 561.37147932, 568.93245083,
    574.94220613, 580.15804762, 585.27208531, 591.28505198, 595.7805288,
    600.94231954, 607.72093532])

y_measured = fr_posttrim/fr_target

def freqmodelsimple(x, fa):
    return fa/np.sqrt(x)

#def freqmodel(x, fa, alpha, beta, gamma, epsilon):
#    return fa/np.sqrt(x)/np.sqrt(1 + alpha + beta*x + gamma*x**2 + epsilon/x)

def freqmodel(x, fa, a, b):
    return fa/np.sqrt(x)/np.sqrt(1 + a*b + a*x + b/x)


#def freqmodel(x, fa, a1, a2, b1):
#    alpha = a1*b1
#    beta = a1 + a2*b1
#    gamma = a2
#    epsilon = b1
#    return fa/np.sqrt(x)/np.sqrt(1 + alpha + beta*x + gamma*x**2 + epsilon/x)
model = freqmodel

falphas = fr_start*np.sqrt(x)

#p0 = [543.5284, 1.5670355, 0.198976]
p0 = [747.38, 2.8502, 0.3509]
xfine = np.r_[0.3:1.0:1000j]

y = old_fr
y = fr_start
popt, pcov = curve_fit(model, x, y, p0=p0)
print (popt)
#popt = [737.57, 2.1751, 0.4598]
yfine = model(xfine, *popt)
res = model(x, *popt) - y

#plt.figure()
#plt.plot(x, old_fr, 'bo', ms=12)
#plt.plot(xfine, yfine, 'k', ms=12)
#plt.plot(x, fr_start, 'rd', ms=12)
#plt.grid()
#plt.show()
#exit()

fa, a, b = popt

fr_bestguessstart = fa/np.sqrt(x)/np.sqrt(1 + a*b + a*x + b/x)
ratio = fr_start/fr_bestguessstart

#sigma = np.diag(np.sqrt(pcov))

f_i = fr_start
x_i = x

Ntries = 100
x_1s = np.zeros((Ntries, x.size))
y_1s = np.zeros((Ntries, x.size))
x_infs = np.zeros((Ntries, x.size))
y_infs = np.zeros((Ntries, x.size))
#for i in range(1):
#
#    x_i = x_i*(f_i/fr_target)**2
#    f_i = fa_shifted/np.sqrt(x_i)/np.sqrt(1 + a_shifted*x_i + b_shifted/x_i)
#    print (x_i/x)
#    print (f_i/fr_start)
#
#    #plt.figure(figsize=(10,10))
#    #plt.plot(fr_target, f_i/fr_target, 'o', ms=12)
#    #plt.grid()
#    #plt.ylim(bottom=0.88, top=1.12)
#    #plt.xlabel('Target Frequency [MHz]')
#    #plt.ylabel('Target Ratio')
#    #plt.show()

for itrial in range(Ntries):
    #p_shifted = 0.5*sigma*np.random.randn(x.size, popt.size)
    p_shifted = 0.1*np.random.multivariate_normal(np.zeros(popt.size), pcov, x.size)

    fa_shifted = fa*ratio + p_shifted[:, 0]
    a_shifted = a + p_shifted[:, 1]
    b_shifted = b + p_shifted[:, 2]

    #xinf = 1./(2*a)*(np.sqrt(1 - 4*a*b + 4*a*(fa*ratio/fr_target)**2) - 1)
    denom = 1#(1 + a*b)
    xinf = denom/(2*a)*(np.sqrt(1 - 4*a*b/denom**2 + 4*a/denom**2*(fa*ratio/fr_target)**2) - 1)
    #yinf = (fa*ratio/fr_target)/np.sqrt(b + xinf + a*xinf**2)
    yinf = (fa_shifted/fr_target)/np.sqrt(b_shifted + xinf*(1 +
        a_shifted*b_shifted) + a_shifted*xinf**2)
    #yinf = (fa*ratio/fr_target)/np.sqrt(xinf)/np.sqrt(1 + a*xinf +
    #        b/xinf)

    #yinf = (fa_shifted/fr_target)/np.sqrt(xinf)/np.sqrt(1 + a_shifted*xinf +
    #        b_shifted/xinf)

    x1 = x*(fr_start/fr_target)**2
    y1 = (fa_shifted/fr_target)/np.sqrt(b_shifted + x1*(1 +
        a_shifted*b_shifted) + a_shifted*x1**2)
    #x1 = 1./(2*a_shifted)*(fr_target/fr_start)**2*(np.sqrt(1 -
    #    4*a_shifted*b_shifted + 4*a_shifted*(fa_shifted/fr_target)**2) - 1)

    x_infs[itrial] = xinf
    y_infs[itrial] = yinf
    x_1s[itrial] = x1
    y_1s[itrial] = y1


y1_mean = np.mean(y_1s, axis=0)
yinf_mean = np.mean(y_infs, axis=0)
y1_std = np.std(y_1s, axis=0)
yinf_std = np.std(y_infs, axis=0)

plt.figure(figsize=(14,10))
#plt.plot(fr_target, yinf, 'o', ms=12, label='best solution')
plt.plot(fr_target, y_measured, 'kd', ms=12, label='CF200713 measured')
plt.errorbar(fr_target, y1_mean, y1_std, fmt='bo', ms=12, label='First approximation')
plt.errorbar(fr_target, yinf_mean, yinf_std, fmt='ro', ms=12, label='Best solution')
#for itrial in range(Ntries):
#    plt.plot(fr_target, y_1s[itrial], 'bo', ms=12, label='%d'%itrial)
#    plt.plot(fr_target, y_infs[itrial], 'ro', ms=12)
plt.grid()
plt.xlabel('$f_t$ [MHz]')
plt.ylabel('$f_r/f_t$')
#plt.ylim(bottom=0.94)
plt.title('Simulated resonator trimming using the full covariance matrix')
plt.legend(loc='upper right')
plt.savefig('simulated_resonator_shift_fullcov.png')
plt.show()

exit()

dx_cumulative = x - x_i
dx_singlestep = x*(1 - (fr_start/fr_target)**2)

plt.figure(figsize=(10,10))
plt.plot(fr_target, dx_cumulative*N0, 'o', ms=12, label='Cumulative')
plt.plot(fr_target, dx_singlestep*N0, 'o', ms=12, label='Single step')
plt.grid()
plt.legend(loc='upper right')
plt.xlabel('Target Frequency [MHz]')
plt.ylabel('Number of Fingers Trimmed')
plt.show()

plt.figure(figsize=(10,10))
plt.plot(fr_target, (dx_cumulative - dx_singlestep)*N0, 'o', ms=12, label='Cumulative')
plt.grid()
plt.legend(loc='upper right')
plt.xlabel('Target Frequency [MHz]')
plt.ylabel('Difference in the 2 solutions')
plt.show()
#plt.figure(figsize=(10,10))
#plt.plot(old_fr, falphas, 'bo', ms=12)
##plt.plot(x, falphas, 'bo', ms=12)
#plt.grid()
##plt.xlabel('N/N0')
#plt.xlabel('fr design [MHz]')
#plt.ylabel('f_alpha [MHz]')
#plt.show()
#
#plt.figure(figsize=(10,10))
#plt.plot(x, y, 'bo', ms=12, label='measured')
#plt.plot(xfine, yfine, 'r')
#plt.grid()
#plt.xlabel('N/N0')
#plt.ylabel('Frequency [MHz]')
#plt.show()
#
#plt.figure(figsize=(10,10))
#plt.plot(x, res, 'bo', ms=12, label='measured')
#plt.grid()
#plt.xlabel('N/N0')
#plt.ylabel('Frequency Residuals [MHz]')
#plt.show()




