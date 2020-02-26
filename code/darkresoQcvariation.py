import numpy as np
import matplotlib.pyplot as pl
import reso_fit
from scipy.interpolate import RectBivariateSpline
pi = np.pi
from scipy import optimize
from scipy.constants import c, epsilon_0

# transmission line
# Nb/Nb, 0.3um ILD, 20um MS, 80um CPW 20um wide slot

Z0 = 50
xline_eps = 13.4 - 0.0081j
xline_n = np.sqrt(xline_eps)
xline_Z0 = 31.
xline_ltot = 0.502
pF = 1e-12
La = 0.0e-6

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
	beta = gamma.imag
	prefactor = 2*(C)/(Cc**2*Z0*w0)*(1 - 2*Cc*La*w0**2)
	print ((1 - 2*Cc*La*w0**2))
	cosh_tot = np.cos(beta*xline_ltot)
	sinh_tot = np.sin(beta*xline_ltot)
	cosh_in = np.cos(2*beta*l)
	sinh_in = np.sin(2*beta*l)
	cosh_out = np.cos(2*beta*(xline_ltot-l))
	sinh_out = np.sin(2*beta*(xline_ltot-l))
	lambd = Zc/Z0
	num = (4*lambd**2*cosh_tot**2 + (1+lambd**2)**2*sinh_tot**2)
	denom = lambd**2*(2*(1+lambd**2) + (lambd**2 - 1)*cosh_in +\
			(lambd**2-1)*sinh_out)
	return prefactor * num/denom

ncf = 11
nfrac = 101
fracs = np.linspace(0,1,nfrac)
cfs = np.linspace(290e6,500e6,ncf)
#Qes = np.zeros(nfrac,dtype=np.complex)

cmap = pl.get_cmap('jet')
design_freqs = np.array([420, 417, 408, 381, 372, 369, 354, 387, 384, 405, 402,
	435, 438, 399, 396, 393, 390, 351, 336, 333, 330, 459, 456, 453, 450, 465,
	462, 327, 324, 339, 306, 321, 318, 471, 468, 483])
argsort = np.argsort(design_freqs)
num_indices = design_freqs.size
pos_indices = np.arange(num_indices)
pos_indices[6:] += 1
pos_indices[2*6:] += 1
pos_indices[3*6:] += 1
pos_indices[4*6:] += 1
pos_indices[5*6:] += 1

L = 10e-9
Cs = 1./(L*(2*pi*design_freqs*1e6)**2)
# Ccs in frequency sorted order
Ccs = np.array([0.23960004, 0.23960004, 0.23960004, 0.23960004, 0.23960004, 0.21205023
, 0.21205023, 0.21205023, 0.21205023, 0.21205023, 0.21205023, 0.21205023
, 0.21205023, 0.21205023, 0.18450043, 0.18450043, 0.18450043, 0.18450043
, 0.18450043, 0.18450043, 0.18450043, 0.18450043, 0.18450043, 0.18450043
, 0.18450043, 0.18450043, 0.15695062, 0.15695062, 0.15695062, 0.15695062
, 0.15695062, 0.15695062, 0.15695062, 0.15695062, 0.15695062, 0.15695062
, 0.15695062, 0.15695062, 0.15695062, 0.15695062, 0.15695062, 0.15695062
, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082
, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082
, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082, 0.12940082
, 0.12940082, 0.12940082, 0.12940082, 0.10185101])*pF/2

Ccs[0] = np.nan
Ccs[1] = np.nan
Ccs[3] = np.nan
Ccs[4] = np.nan
Ccs[5] = np.nan
Ccs[14] = np.nan
Ccs[15] = np.nan
Ccs[16] = np.nan
Ccs[19] = np.nan
Ccs[20] = np.nan
Ccs[21] = np.nan
Ccs[22] = np.nan
Ccs[25] = np.nan
Ccs[26] = np.nan
Ccs[37] = np.nan
Ccs[38] = np.nan
Ccs[41] = np.nan
Ccs[42] = np.nan
Ccs[43] = np.nan
Ccs[44] = np.nan
Ccs[47] = np.nan
Ccs[48] = np.nan
Ccs[49] = np.nan
Ccs[58] = np.nan
Ccs[59] = np.nan
Ccs[60] = np.nan
Ccs[62] = np.nan
Ccs[63] = np.nan
mask = np.logical_not(np.isnan(Ccs))
Ccs = Ccs[mask]

mapper = dict(zip(design_freqs[argsort], Ccs))
#print (mapper)

ordered_Ccs = np.array(list(map(lambda x: mapper[x], design_freqs)))
#fig, ax = pl.subplots(figsize=(10,10))
#ax.plot(design_freqs, ordered_Ccs/pF, 'ro')
#ax.plot(design_freqs[argsort], Ccs/pF, 'bh')
#pl.show()
#exit()

mm = 1e-3
termination = 25.5*mm
reso_segment = 11*mm
pos = termination + (pos_indices + 0.5)*reso_segment
total_len = 2*termination + (num_indices + 5)*reso_segment
frac_position = pos/total_len

print (frac_position)
Qcs_pred = qc_prediction(design_freqs*1e6,Cs, ordered_Ccs, xline_Z0, xline_n,
		xline_ltot*(1 - frac_position))
data = np.vstack([design_freqs, frac_position, Qcs_pred]).T

fig,ax = pl.subplots(figsize=(10,10))
ax.plot(design_freqs, Qcs_pred, 'ro')
ax.grid()
ax.set_xlabel('Design Frequency')
ax.set_ylabel('Predicted Qc')
pl.show()

np.savetxt('analytical_Qc_%d_0p2nH.txt'%(xline_Z0), data)

