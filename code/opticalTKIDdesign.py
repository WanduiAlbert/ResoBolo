#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.special import iv, kn
from scipy.constants import h, k
from waffleTKIDcapacitordesign import IDC

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

nH = 1e-9
pF = 1e-12
MHz = 1e6
cm = 1e-2
g = 1e-3
mJ = 1e-3

makeplots = False
# Aluminum material properties
gamma = 1.35 * mJ#/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g#/u.mol
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k**2))

Z0 = 50
L = 10*nH

# Rounds numbers to the nearest base
def roundto(num, base):
	return ((num // base) + ((num % base) > (base//2))) * base

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


def get_simulated_LC(Nfingers):
	C_sim = 0.05*Nfingers + 0.788
	L_sim = 6.1e-6*Nfingers**2 + 1e-3*Nfingers + 1.780
	return C_sim, L_sim


if __name__=="__main__":
	df = 5 #MHz
	N = 12
	fstart = 300 #MHz
	fr = fstart + np.arange(N) * df
	target_T = 380e-3
	Qi = get_MB_Qi(target_T, fr*MHz)
	#print (Qi)
	#print (Qi)
	wr = 2*pi*fr*MHz
	C = 1./(wr**2*L)
	C_branch = 0.11 * pF
	C_load = C_branch * N
	L_branch = 0.0 * nH
	L_load = L_branch * N
	Zin = Z0
	Zout = 1/(1./Z0 + 1j*wr*C_load + 1./(1j*wr*L_load))
	Cc = np.sqrt((2*C)/(Qi*Z0*wr))
	CC = np.average(Cc) # use the same coupling cap for all the resonators
	#CC = 0.1945 * pF # Using a number from actual calculations
	y = wr * CC * Z0
	Gprime = (wr*CC*Zin*y/(Zin + Zout)) - (1j*wr*CC*Zin**2*y**2)/(Zin + Zout)**2
	dQe = Gprime/(wr*C)
	Qe = 1/dQe
	Qc = 1./np.real(dQe)
	#print (Qc)
	#print (Qc)
	phi_c = np.arctan2(dQe.imag, dQe.real)
	#print (phi_c*180/pi)
	#Qc = (C*pF)/(0.5*Z0*wr*(CC*pF/2)**2)
	Qr = 1./(1./Qc + 1./Qi)
	#chi_c = 4*Qr**2/Qi/Qc/np.cos(phi_c)
	#chi_c = 4*np.abs(Qe)*Qi/(Qi*np.cos(phi_c) + np.abs(Qe))**2
	chi_c = 4*Qc*Qi/(Qi + Qc)**2#/np.cos(phi_c)

	plt.figure()
	plt.scatter(np.arange(N), chi_c)
	plt.xlabel('Resonator Index')
	plt.ylabel('Coupling Efficiency chi_c')
	plt.grid()
	plt.savefig('opticalTKID_coupling_efficiency.png')
	plt.close()
	dQr = 1/Qr

	#exit()
	if makeplots:
		for i in range(N):
			frequency = np.linspace(-1.0, 1.0, 1000) + fr[i]
			x = (frequency - fr[i])/fr[i]
			S21 = 1 - dQe[i]/(dQr[i] + 1j * 2 * x)
			label = 'Qr=%d Qi=%d\nQc=%d phi_c=%1.3fdeg'%(Qr[i], Qi[i], Qc[i],
					phi_c[i]*180/pi)
			S21dB = 10*np.log10(np.abs(S21))


			plt.figure(i)
			plt.plot(frequency, S21dB, 'b', label=label)
			plt.xlabel('Frequency [MHz]')
			plt.ylabel('|S21|')
			plt.legend(loc='center right')
			plt.grid()
			plt.savefig("reso_%1.3fMHz_S21.png"%fr[i])
			plt.close()

			plt.figure(i + 250)
			plt.plot(S21.real, S21.imag, 'b', label=label)
			plt.xlabel('I')
			plt.ylabel('Q')
			plt.legend(loc='center right')
			plt.grid()
			plt.axis('square')
			plt.savefig("reso_%1.3fMHz_IQ.png"%fr[i])
			plt.close()
	#exit()

	print ("Total Coupling Capacitance Needed ", CC/pF)
	CC *= 2 #Needed because we have 2 coupling capacitors in series
	# I want to pick the largest capacitor as the base IDC and then design blade
	# cells to create the other capacitors
	Cmain = np.max(C)
	print (C/pF)
	cap = IDC(1.0)
	trace_width = 2.
	gap_width = 2.
	finger_length = 1000.
	finger_gap = 2.
	contact_width = 2
	nfinger = 10 # For now
	cap.set_dimensions(trace_width, gap_width,
			finger_length, finger_gap, nfinger, contact_width)

	# Want to determine the number of fingers for this capacitor structure.
	Nfingers, capfrac = cap.getnumfingers(Cmain)
	print ("Expected number of fingers is ", Nfingers)
	print ("Expected number of fractional fingers is ", capfrac)

	cap.set_dimensions(trace_width, gap_width,
			finger_length, finger_gap, Nfingers, contact_width)
	C_actual = cap.capacitance()/pF
	print ("The expected capacitance of the structure is %1.3f pF"%C_actual)
	print ("Width/height: (%1.3f, %1.3f) um"%(cap.width, cap.height))

