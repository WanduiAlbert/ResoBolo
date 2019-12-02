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
xline_Z0 = 30.
xline_ltot = 0.104086
pF = 1e-12

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
	prefactor = 2*(C)/(Cc**2*Z0*w0)#*(1 + (Cc*Z0*w0/2)**2)
	cosh_tot = np.real(np.cosh(gamma*xline_ltot))
	sinh_tot = -np.imag(np.sinh(gamma*xline_ltot))
	cosh_in = np.real(np.cosh(2*gamma*l))
	sinh_in = -np.imag(np.sinh(2*gamma*l))
	cosh_out = np.real(np.cosh(2*gamma*(xline_ltot-l)))
	sinh_out = -np.imag(np.sinh(2*gamma*(xline_ltot-l)))
	lambd = Zc/Z0
	num = -(4*lambd**2*cosh_tot**2 + (1+lambd**2)**2*sinh_tot**2)
	denom = lambd**2*(-2*(1+lambd**2) + (-1+lambd**2)*cosh_in +\
			(1-lambd**2)*cosh_out)
	return prefactor * num/denom

cmap = pl.get_cmap('jet')
design_freqs = np.array([266.1, 271.7, 276.9, 282.4, 287.6, 292.7, 298.0, 303.5,
	308.3, 313.3, 318.5, 323.4])
num_indices = design_freqs.size
pos_indices = np.arange(num_indices)

arrayed_pos = pos_indices[::-1].reshape((3, -1))
arrayed_pos[0] = arrayed_pos[0][::-1]
arrayed_pos[1] = arrayed_pos[1][::-1]
arrayed_pos[-1] = arrayed_pos[-1][::-1]
arrayed_pos = arrayed_pos.flatten()

dist = np.zeros(design_freqs.size)
Cs = np.array([28.14, 27.23, 26.36, 25.53, 24.74, 23.98, 23.26, 22.57, 21.91,
	21.28, 20.68, 20.10])*pF
L = 1./(Cs*(2*pi*design_freqs*1e6)**2)
Cc = 0.1883*pF

mm = 1e-3
termination = (0.25 + 1.505+8.754 + 3.580)*mm
hor_segment = 17.508*mm
vert_segment = 11.692*mm
pair_spacing = 4.200*mm
interpair_spacing = 8.754*mm
gap = (hor_segment - interpair_spacing - pair_spacing)/2
dist += termination + gap
dist[1:] += pair_spacing
dist[2:] += (interpair_spacing - pair_spacing)
dist[3:] += pair_spacing
dist[4:] += 2*gap + vert_segment
dist[5:] += pair_spacing
dist[6:] += (interpair_spacing - pair_spacing)
dist[7:] += pair_spacing
dist[8:] += 2*gap + vert_segment
dist[9:] += pair_spacing
dist[10:] += (interpair_spacing - pair_spacing)
dist[11:] += pair_spacing
frac_position = dist/xline_ltot

print (frac_position)
Qcs_pred = qc_prediction(design_freqs[arrayed_pos]*1e6,Cs[arrayed_pos], Cc, xline_Z0, xline_n,
		xline_ltot*(frac_position))

freqmeas = np.array([298.0, 303.5, 313.3, 323.4])
Qcmeas = np.array([31017, 28325, 14796, 10527])
fig,ax = pl.subplots(figsize=(10,10))
ax.plot(design_freqs[arrayed_pos], Qcs_pred, 'ro')
ax.plot(freqmeas, Qcmeas, 'bh')
ax.grid()
ax.set_xlabel('Design Frequency')
ax.set_ylabel('Predicted Qc')
pl.show()
