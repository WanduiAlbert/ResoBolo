
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, sys
from scipy.special import iv, kn
from scipy.constants import h, k,c
import scipy.optimize as opt

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

nH = 1e-9
pF = 1e-12
cm = 1e-2
g = 1e-3
mJ = 1e-3

MHz = 1e6
mm = 1e-3
um = 1e-6

datadir = "../numerical_sims/"
spacing = 11*mm

warm_attenuation = 20
cold_attenuation = 50
total_attenuation = warm_attenuation + cold_attenuation
P = np.arange(-40,2,2) - total_attenuation

# Aluminum material properties
gamma = 1.35 * mJ#/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g#/u.mol
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k**2))

Z0 = 50
L = 10*nH

def get_MB_Qi(T, f):
	Q_int = 3e7
	Tc = 1.4
	alpha = 0.5
	Delta = 1.764 * k * Tc
	eta = (h * f / (2 * k * T))
	S_1 = ((2/pi)*np.sqrt(2*Delta/(pi*k*T))* np.sinh(eta)*K_0(eta))
	n_th = (2*N_0 * np.sqrt(2*pi* k * T* Delta)*np.exp(-Delta/(k*T)))
	Q_qp = ((2 * N_0 * Delta)/(alpha * S_1 * n_th))
	Q_i = 1./(1/Q_qp + 1./Q_int)

	return Q_i

# Initial design of the resonators
Nreso = 12
df = 5 #MHz
fstart = 300 #MHz
fr = fstart + np.arange(Nreso) * df
target_T = 380e-3
Qi = get_MB_Qi(target_T, fr*MHz)
wr = 2*pi*fr*MHz
C = 1./(wr**2*L)
C_branch = 0.11 * pF
C_load = C_branch * Nreso
Zin = Z0
Zout = 1/(1./Z0 + 1j*wr*C_load)
Cc = np.sqrt((2*C)/(Qi*Z0*wr))
CC = np.average(Cc) # use the same coupling cap for all the resonators
print (CC/pF)
Ca = 0.403*pF # from simulations
CC = Ca/2 # from simulations
La = 0.0*nH
w_star = 1./np.sqrt(L*CC)
f_star = w_star/2/pi/MHz

qc_prefactor = 2*w_star*L/Z0

Cpar = 0.1*pF
factor = 1 + (Ca*Z0*wr/4)**2*((Ca + 2*Cpar)/(Ca + Cpar))**2
#Qc = (8*C/(Ca**2*Z0*wr))*factor
#CC *= 1.5
y = wr * CC * Z0
Gprime = (wr*CC*Zin*y/(Zin + Zout)) - (1j*wr*CC*Zin**2*y**2)/(Zin + Zout)**2
dQe = Gprime/(wr*C)
Qe = 1/dQe
#Qc = 1./np.real(dQe)
Qc = qc_prefactor*(fr/f_star)**-3
Qc = (C/Ca)*(8/(Ca*Z0*wr) - 8*(La*wr/Z0) + wr*Ca*Z0/2 + 2*La**2*wr**3*Ca/Z0)
phi_c = np.arctan2(dQe.imag, dQe.real)
Qr = 1./(1./Qc + 1./Qi)
chi_c = 4*Qc*Qi/(Qi + Qc)**2#/np.cos(phi_c)


# Adjusting for the results from the sonnet simulations
#Nfingers = 502 - np.array([0, 17, 32, 47, 61, 74, 87, 100, 111, 122, 133,
#	153, 172, 189, 205, 219])
Nfingers = 502 - np.array([0, 17, 32, 47, 61, 74, 87, 100, 111, 122, 133,
	153])
C_true = (0.05*Nfingers + 0.788)*pF
Lpar = (6e-6*Nfingers**2 + 1e-3*Nfingers + 1.780)*nH
Lind = 10.057*nH #5.998*nH with no kinetic inductance
Ltot = Lind + Lpar
wr_true = 1./(Ltot*C_true)**0.5
fr_true = wr_true/2/pi/MHz
y_true = wr_true * CC * Z0
Zout_true = 1/(1./Z0 + 1j*wr_true*C_load)
Gprime_true = (wr_true*CC*Zin*y_true/(Zin + Zout_true)) -\
		(1j*wr_true*CC*Zin**2*y_true**2)/(Zin + Zout_true)**2
dQe_true = Gprime_true/(wr_true*C_true)
Qe_true = 1/dQe_true
Qc_true = 1./np.real(dQe_true)
#Qc_true = 2*C_true/(wr_true*Z0*CC**2)
factor = 1 + (Ca*Z0*wr_true/4)**2*((Ca + 2*Cpar)/(Ca + Cpar))**2
Qc_true = (8*C_true/(Ca**2*Z0*wr_true))*factor
Qc_true = (C_true/Ca)*(8/(Ca*Z0*wr_true) - 8*(La*wr_true/Z0) + wr_true*Ca*Z0/2 + 2*La**2*wr_true**3*Ca/Z0)
phi_c_true = np.arctan2(dQe_true.imag, dQe_true.real)
Qi_true = get_MB_Qi(target_T, fr_true*MHz)
Qr_true = 1./(1./Qc_true + 1./Qi_true)
chi_c_true = 4*Qc_true*Qi_true/(Qi_true + Qc_true)**2#/np.cos(phi_c)

freq_meas = np.mean(np.loadtxt(datadir + 'CF190822_f0_vs_P.txt'), axis=1)
Qc_meas = np.mean(np.loadtxt(datadir + 'CF190822_Qcs_vs_P.txt'), axis=1)
phic_meas = np.mean(np.loadtxt(datadir + 'CF190822_phics_vs_P.txt'), axis=1)
resoindices = np.array([0,1,2,3,4,6,8,9])
fr_fine1 = np.r_[230:290:1000j]
#Qc_guess = 2*C_true[resoindices]/(Z0*2*pi*freq_meas*MHz*CC**2)
p1 = np.polyfit(np.log(freq_meas/f_star), np.log(Qc_meas), 1)
Qc_fit = np.exp(np.polyval(p1, np.log(fr_fine1/f_star)))
print (p1)

freq_meas2 = np.mean(np.loadtxt(datadir + 'CF69C_f0_vs_P.txt')[:, :5], axis=1)
Qc_meas2 = np.mean(np.loadtxt(datadir + 'CF69C_Qcs_vs_P.txt')[:, :5], axis=1)
phic_meas2 = np.mean(np.loadtxt(datadir + 'CF69C_phics_vs_P.txt')[:, :5], axis=1)
resoindices2 = np.array([0,1,4,5,7,9])
#Qc_guess2 = 2*C_true[resoindices2]/(Z0*2*pi*freq_meas2*MHz*CC**2)
fr_fine2 = np.r_[270:330:1000j]
mask = np.ones_like(freq_meas2, dtype=bool)
#mask[0] = False
mask[1] = False
#mask[3] = False
#print (freq_meas2[mask])
p2 = np.polyfit(np.log(freq_meas2[mask]/f_star), np.log(Qc_meas2[mask]), 1)
Qc_fit2 = np.exp(np.polyval(p2, np.log(fr_fine2/f_star)))
fr_fine = np.r_[230:330:1000j]
print (p2)

#def get_Qc(fr):
#	wr = 2*pi*fr*MHz
#print (Qc_meas)
#print (Qc_meas2)
wr_meas = 2*pi*freq_meas*MHz
wr_meas2 = 2*pi*freq_meas2*MHz
actual_L = 1./wr_meas**2/C_true[resoindices]
actual_C = 1./wr_meas**2/Ltot[resoindices]
#print (C_true/pF)
#print (Lpar/nH)
#print (fr_true)
#print (actual_L/nH)
#print (Ltot[resoindices]/nH)
#print ((actual_L - Ltot[resoindices])/nH)
#print (actual_C/pF)
#print (C_true[resoindices]/pF)
#print ((actual_C - C_true[resoindices])/pF)
#exit()

wr_fine = 2*pi*fr_fine*MHz
L_fine = np.average(Ltot)
Qc_fine = 2./(L_fine*Z0*wr_fine**3*CC**2)


#def Qc_model(f, wstar, fref):
#	w = 2*pi*(f-fref)
#	return 2*L*wstar*MHz/Z0*(w/wstar)**-3
def Qc_model(f, wstar, fref):
	w = 2*pi*f
	return (2*L*wstar*MHz/Z0*(w/wstar)**-3)*(1 - (f/fref)**2)
popt, pcov = opt.curve_fit(Qc_model, freq_meas, Qc_meas, p0=[w_star/MHz, 1000])
sigma = np.sqrt(np.diag(pcov))
fstar = popt[0]/2/pi
sigfstar = sigma[0]/2/pi
fref = popt[1]
sigfref = sigma[1]
popt2, pcov2 = opt.curve_fit(Qc_model, freq_meas2[mask], Qc_meas2[mask],
		p0=[w_star/MHz, 1000])
sigma2 = np.sqrt(np.diag(pcov2))
fstar2 = popt2[0]/2/pi
sigfstar2 = sigma2[0]/2/pi
fref2 = popt2[1]
sigfref2 = sigma2[1]

Ca1 = 1./((popt[0]*MHz)**2*L/2)/pF
Ca2 = 1./((popt2[0]*MHz)**2*L/2)/pF

print ("\n\n\n\n\n")
print (fstar, sigfstar)
print (fref, sigfref)
print (Ca1, 2*Ca1/popt[0]*sigma[0])
print ("\n\n\n\n\n")
print (fstar2, sigfstar2)
print (fref2, sigfref2)
print (Ca2, 2*Ca2/popt[0]*sigma2[0])

Qc_meas_fit = Qc_model(fr_fine1, *popt)
Qc_meas_fit2 = Qc_model(fr_fine2, *popt2)


# Maybe getting more desperate. What happens if line mismatch is actually
# important in this case.
la = 2000*um
lb = 8754*um
lc = 3580*um
ld = 17508*um
le = 11672*um
Delta = 4200*um
delta = 2281*um
gamma = 4554*um
lt = 2*(la + lb + lc + le) + 3*ld

pos = np.zeros(Nreso)
pos[0] = la + lb + lc + ld - delta
pos[1] = pos[0] - Delta
pos[2] = pos[0] - (Delta + gamma)
pos[3] = pos[0] - (2*Delta + gamma)
pos[4] = pos[0] + le + 2*delta
pos[5] = pos[4] + Delta
pos[6] = pos[4] + (Delta + gamma)
pos[7] = pos[4] + (2*Delta + gamma)
pos[8] = pos[7] + le + ld
pos[9] = pos[8] - Delta
pos[10] = pos[8] - (Delta + gamma)
pos[11] = pos[8] - (2*Delta + gamma)

# For reverse mapping x -> lt - x
#pos = lt - pos

er = 7.95
n = er**0.5
wt = c/(lt*n)
lin = np.linspace(0.1,0.9,wr_fine.size)*lt + 1e-9
win = c/(lin*n)
#win = c/(pos*n)
k = pos/lt
#print (2*k-1)
#print (wr_fine/MHz)
#print (win/MHz)
#print ((2*wt - win)/win)
zlambda = 0.98

num = 4*zlambda**2*np.cos(wr_fine/wt)**2 + (1 + zlambda**2)**2*np.sin(wr_fine/wt)**2
denom = 2*zlambda**2*((1+zlambda**2) -
		(-1+zlambda**2)*np.cos(wr_fine/wt)*np.cos(wr_fine/wt*(2*wt-win)/win))
Qcfeedline = (8/(L*wr_fine**3*Z0*Ca**2))*num/denom
Qcfeedline = 1/(wr_fine**2*L*Ca)*(8/(Ca*Z0*wr_fine) - 8*(La*wr_fine/Z0) + wr_fine*Ca*Z0/2 + 2*La**2*wr_fine**3*Ca/Z0)


plt.figure(figsize=(10,10))
#plt.plot(fr_fine, Qcfeedline, 'r-')
#plt.plot(fr_fine1, Qc_fit, 'k--')
#plt.plot(fr_fine2, Qc_fit2, 'k--')
plt.plot(fr_fine1, Qc_meas_fit, 'k--')
plt.plot(fr_fine2, Qc_meas_fit2, 'k--')
#plt.plot(fr_true, Qc, 'ks', label='Design')
#plt.plot(freq_meas, Qc_true[resoindices], 'bo', label='Modified Design')
plt.plot(freq_meas, Qc_meas, 'rh', label='CF190822 Measured')
plt.plot(freq_meas2, Qc_meas2, 'bh', label='CF69C Measured')
plt.grid()
plt.legend(loc='best')
plt.xlabel('Freq [MHz]')
plt.ylabel('Qc')
plt.savefig('modifiedQc_vs_measuredQc.png')

plt.show()
exit()

plt.figure(figsize=(10,10))
plt.plot(freq_meas, Qc_meas/Qc_true[resoindices], 'ro', label='CF190822')
plt.plot(freq_meas2, Qc_meas2/Qc_true[resoindices2], 'bh', label='CF69C')
plt.grid()
plt.legend(loc='best')
plt.xlabel('Freq [MHz]')
plt.ylabel('scaling of Qc between measured and design')
plt.savefig('measuredQc_vs_scaling.png')

plt.show()
exit()


predanalytical_Qc = np.loadtxt('analytical_Qc_31.txt', unpack=True)
_, frac_pos, pred_Qc = np.loadtxt('predicted_mismatchedQc.txt', unpack=True)
_, frac_pos30, pred_Qc30 = np.loadtxt('predicted_mismatchedQc_30.txt', unpack=True)
_, frac_pos50, pred_Qc50 = np.loadtxt('predicted_mismatchedQc_50.txt', unpack=True)
fr = freq[:, -1]

deltaQi = Qi[:, -1]
x = (freq[:, :-1] - freq[:, -1][:, np.newaxis])/freq[:, -1][:, np.newaxis] * 1e6

z = np.zeros_like(x.T)
z[x.T == np.min(x.T, axis=0)] = 1

reso2led = np.argmax(z, axis=0)
# some manual adjustments
reso2led[16] = 35
reso2led[21] = 10

#print (fr/MHz)
#print (reso2led)

ordered_resos = np.zeros(NLEDs)
ordered_resos[reso2led] = fr/MHz

ordered_Qis = np.zeros(NLEDs)
ordered_Qis[reso2led] = deltaQi
ordered_Qis[ordered_Qis == 0] = np.nan
arrayed_Qis = ordered_Qis.reshape((6,-1)).T
Qis = arrayed_Qis.flatten()

arrayed_resos = ordered_resos.reshape((6, -1)).T
design_freqs = np.array([306, 321, 318, 471, 468, 483, 339, 324, 327, 462, 465,
	450, 336, 333, 330, 459, 456, 453, 351, 390, 393, 396, 399, 438, 354, 387,
	384, 405, 402, 435, 369, 372, 381, 408, 417, 420])
design_freqs2 = np.array([420, 417, 408, 381, 372, 369, 354, 387, 384, 405, 402,
	435, 438, 399, 396, 393, 390, 351, 336, 333, 330, 459, 456, 453, 450, 465,
	462, 327, 324, 339, 306, 321, 318, 471, 468, 483])

pos_indices = np.array([31, 32, 33, 34, 35, 36, 30, 29, 28, 27, 26, 25, 19, 20, 21,
	22, 23, 24, 18, 17, 16, 15, 14, 13, 7, 8, 9, 10, 11, 12, 6, 5, 4, 3, 2,
	1])-1

# The data from the Qc predictions is sorted according to the position along the
# feedline. I want to remap this to being ordered according to the indexing
# scheme of the array
mapper = sorted(zip(np.arange(pos_indices.size), pos_indices))
mapping_order = np.array([x for _, x in mapper])

termination = 25.5*mm
pos = termination + (pos_indices + 0.5)*spacing
total_len = 2*termination + (design_freqs.size + 5)*spacing
frac_position = pos/total_len

#design_freqs = 1.0*design_freqs - 38.9
arrayed_designfreqs = design_freqs.reshape((6, -1))

arrayed_delta = (arrayed_resos - arrayed_designfreqs)/arrayed_designfreqs
arrayed_delta[arrayed_delta == -1] = np.nan
delta = arrayed_delta.flatten()
mask = np.isfinite(delta)
print (np.mean(delta[mask]))
print (np.std(delta[mask]))

Qes = np.zeros_like(design_freqs, dtype=np.complex64)
Qes[reso2led] = np.mean(Qe, axis=1)
Qe = Qes.reshape((6, -1)).T
dQe = 1./Qe
Qc = 1./np.real(dQe)
phic = -np.arctan2(dQe.imag, dQe.real)
Qc[Qc == 0] = np.nan
phic[Qc == 0] = np.nan

X, Y = np.meshgrid(np.arange(Nrows+1), np.arange(Ncols+1))
#plt.figure(figsize=(12,10))
#plt.pcolormesh(X, Y, Qc)
#plt.title(r'$Q_c$ spread across the wafer')
##plt.xticks([])
#plt.xlabel('X')
##plt.yticks([])
#plt.ylabel('Y')
#plt.colorbar(label=r'$Q_c$')
#plt.savefig('waferuniformity_in_Qc.png')
#
#plt.figure(figsize=(12,10))
#plt.pcolormesh(X, Y, phic)
#plt.title(r'$\phi_c$ spread across the wafer')
##plt.xticks([])
#plt.xlabel('X')
##plt.yticks([])
#plt.ylabel('Y')
#plt.colorbar(label=r'$\phi_c$')
#plt.savefig('waferuniformity_in_phic.png')


fullmask = np.logical_and(mask, delta < 1.0)

pat = np.hstack([[0,1]*3, [1,0]*3])
release = np.hstack([pat, pat, pat]).astype(bool)

masked_released = np.logical_and(fullmask, release)
masked_unreleased = np.logical_and(fullmask, np.logical_not(release))

row, col = np.meshgrid(np.arange(Nrows), np.arange(Ncols), indexing='ij')
X0, Y0 = 1.52*mm, 0*mm
#X0, Y0 = -spacing/2, -spacing/2
X = X0 - (Ncols - 2*col - 1)/2.0 * spacing
Y = Y0 - (Nrows - 2*row - 1)/2.0 * spacing

R = np.sqrt((X)**2 + (Y)**2)
radial_pos = R.flatten()
x_pos = X.flatten()
y_pos = Y.flatten()

#plt.figure()
#plt.plot(x_pos[masked_released]/mm, Qc.flatten()[masked_released], 'bo',
#		label='Released', ms=12)
#plt.plot(x_pos[masked_unreleased]/mm, Qc.flatten()[masked_unreleased], 'rh',
#		label='Unreleased', ms=12)
#plt.xlabel('X [mm]')
#plt.ylabel(r'$Q_c$')
#plt.grid()
#plt.title('X Dependence of Qc')
#plt.legend(loc='upper left')
#plt.savefig('Qc_xdependence.png')
#
#plt.figure()
#plt.plot(y_pos[masked_released]/mm, Qc.flatten()[masked_released], 'bo',
#		label='Released', ms=12)
#plt.plot(y_pos[masked_unreleased]/mm, Qc.flatten()[masked_unreleased], 'rh',
#		label='Unreleased', ms=12)
#plt.xlabel('Y [mm]')
#plt.ylabel(r'$Q_c$')
#plt.grid()
#plt.title('Y Dependence of Qc')
#plt.legend(loc='upper left')
#plt.savefig('Qc_ydependence.png')

plt.figure()
#plt.plot(frac_position[masked_released], Qc.flatten()[masked_released], 'bo',
#		label='Released', ms=12)
#plt.plot(frac_position[masked_unreleased], Qc.flatten()[masked_unreleased], 'rh',
#		label='Unreleased', ms=12)
plt.plot(frac_pos, predanalytical_Qc, 'gP', label='Analytical Prediction: 31 Ohm line impedance',
		ms=12)
#plt.plot(frac_pos, pred_Qc, 'gP', label='Predicted: 38 Ohm line impedance',
#		ms=12)
#plt.plot(frac_pos30, pred_Qc30, 'bP', label='Predicted: 30 Ohm line impedance',
#		ms=12)
#plt.plot(frac_pos50, pred_Qc50, 'kP', label='Predicted: 50 Ohm line impedance',
#		ms=12)
plt.plot(frac_position[fullmask], Qc.flatten()[fullmask], 'rh',
		label='data', ms=12)
plt.xlabel('Fractional Position along Feedline')
plt.ylabel(r'$Q_c$')
plt.grid()
plt.xlim((0,1))
plt.title('Qc variation along the feedline')
lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.savefig('Qc_feedlinevariation.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

plt.figure()
#plt.plot(frac_position[masked_released], Qc.flatten()[masked_released], 'bo',
#		label='Released', ms=12)
#plt.plot(frac_position[masked_unreleased], Qc.flatten()[masked_unreleased], 'rh',
#		label='Unreleased', ms=12)
plt.plot(design_freqs2, predanalytical_Qc, 'gP', label='Predicted: 31 Ohm line impedance',
		ms=12)
#plt.plot(design_freqs2, pred_Qc30, 'bP', label='Predicted: 30 Ohm line impedance',
#		ms=12)
#plt.plot(design_freqs2, pred_Qc50, 'kP', label='Predicted: 50 Ohm line impedance',
#		ms=12)
plt.plot(design_freqs[fullmask], Qc.flatten()[fullmask], 'rh',
		label='data', ms=12)
plt.xlabel('Design Frequency [MHz]')
plt.ylabel(r'$Q_c$')
plt.grid()
#plt.xlim((0,1))
plt.title('Qc variation along the feedline')
lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.savefig('Qc_vsdesignfreq.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

plt.figure()
#plt.plot(frac_position[masked_released], Qc.flatten()[masked_released], 'bo',
#		label='Released', ms=12)
#plt.plot(frac_position[masked_unreleased], Qc.flatten()[masked_unreleased], 'rh',
#		label='Unreleased', ms=12)
plt.plot(Qc.flatten()[fullmask], np.roll(predanalytical_Qc[mapping_order][fullmask], 2), 'rh', label='Predicted: 31 Ohm line impedance', ms=12)
#plt.plot(Qc.flatten()[fullmask], np.roll(pred_Qc30[mapping_order][fullmask], 2), 'bh', label='Predicted: 30 Ohm line impedance', ms=12)
#plt.plot(Qc.flatten()[fullmask], pred_Qc50[mapping_order][fullmask], 'kP', label='Predicted: 50 Ohm line impedance', ms=12)
plt.xlabel('Measured Qc')
plt.ylabel(r'Predicted Qc')
plt.grid()
plt.title('Qc variation along the feedline')
lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.savefig('Qc_measured_vs_pred.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

plt.show()
#print (Qc.flatten()[fullmask])
#print (np.roll(pred_Qc30[mapping_order][fullmask], 2))
#plt.close('all')

