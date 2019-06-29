#! /usr/bin/env python3

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.special import iv, kn
from scipy.constants import h, k

K0 = lambda x: kn(0, x)
I0 = lambda x: iv(0, x)

nH = 1e-9
pF = 1e-12
MHz = 1e6
cm = 1e-2
g = 1e-3
mJ = 1e-3
pW = 1e-12
aW = 1e-18

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
	S_1 = ((2/pi)*np.sqrt(2*Delta/(pi*k*T))* np.sinh(eta)*K0(eta))
	n_th = (2*N_0 * np.sqrt(2*pi* k * T* Delta)*np.exp(-Delta/(k*T)))
	Q_qp = ((2 * N_0 * Delta)/(alpha * S_1 * n_th))
	Q_i = 1./(1/Q_qp + 1./Q_int)

	return Q_i

def get_responsivity(fr, P, Tb, K, n, Qi, Tc=1.284):
	"""
	function: get_responsivity(fr, P, Tb, K, n, Qi, Tc=1.385)

	args: fr - resonance frequency in Hz
		   P - Power dissipated in the heater in W
		  Tb - Temperature of the bath
		   K - Conductance of the bolometer legs. In W/K^(n+1)
		   n - Leg exponent. Defined by P = K (T^(n+1) - Tb^(n+1))
		  Qi - Quality factor of the resonator
		  Tc - Superconducting transition temperature

	returns: responsivity in Hz/W
	"""
	MHz = 1e6
	pW = 1e-12
	eV = 1.6e-19 # Joules
	um = 1e-6
	N0 = 1.72e10/um**3/eV
	alpha_k = 0.378
	nu = fr

	T = (P/K + Tb**(n+1))**(1./(n + 1)) # island temperature
	G = K*(n+1)*T**n
	delta = 1.763*k*Tc
	kappa = 0.5 + delta/(k*T)
	q = (h*nu/2/k/T)
	nqp = 2*N0*np.sqrt(2*pi*delta*k*T)*np.exp(-delta/(k*T))
	S1 = 2/pi * np.sqrt(2*delta/(pi*k*T))*np.sinh(q)*K0(q)
	S2 = 1 + np.sqrt(2*delta/(pi*k*T))*np.exp(-q)*I0(q)
	beta = S2 / S1

	x = alpha_k * S2 * nqp / (4*N0*delta)
	S = (beta * kappa)/(2*Qi*G * T)

	return S



def get_simulated_LC(Nfingers):
	C_sim = 0.05*Nfingers + 0.788
	L_sim = 6.1e-6*Nfingers**2 + 1e-3*Nfingers + 1.780
	return C_sim, L_sim


if __name__=="__main__":

	fr = 337.4
	n = 1.862
	K = 122*pW#/K^(n+1)
	nu = fr*MHz
	Tbath = 80e-3

	Ntemps = 115
	Npowers = 110

	temps = np.linspace(Tbath*1e3, 500, 115)*1e-3
	#Pread = -110 + 2*np.arange(16)
	Pread = np.linspace(-110, -80, 100)

	Ts, Pg = np.meshgrid(temps, Pread)
	Ps = 1e-3*10**(Pg/10)

	#T = 380e-3
	Pthermal = K*(Ts**(n+1) - Tbath**(n+1))
	Qi = get_MB_Qi(Ts, fr*MHz)
	wr = 2*pi*nu
	Qc = 15346
	Qr = 1./(1./Qc + 1./Qi)
	E = (2*Qr**2/Qc)*(Ps/wr)

	Ptotal = Pthermal + Ps*(2*Qr**2/Qc/Qi)
	N = E/(h*nu) #Microwave photon number

	responsivity = get_responsivity(nu, Ptotal, temps, K, n, Qi, Tc=1.284)

	Atls = 1.2e-17
	# I want the noise at 1 Hz
	nu1 = 1
	nu2 = 1e3
	T2 = 120e-3
	beta = 1.8
	Stls = Atls/np.sqrt(N)*(nu1/nu2)**(-0.5)*(temps/T2)**(-beta)
	NEFtls = nu*np.sqrt(Stls)
	NEPtls = np.sqrt(Stls)/responsivity

	fig, ax = plt.subplots(figsize=(10,10))
	#img = ax.pcolormesh(Pthermal/pW, Pg, np.log(NEPtls/aW), cmap="viridis")
	img = ax.pcolormesh(Pthermal/pW, Pg, NEPtls/aW, cmap="viridis")
	ax.set_xlabel('Popt [pW]')
	ax.set_ylabel('Pread [dBm]')
	ax.grid()
	ax.set_title("Tbath = %1.1f mK"%(Tbath*1e3))
	fig.colorbar(img, label='NEP TLS(1 Hz) [aW/rtHz]')
	plt.show()

	fig, ax = plt.subplots(figsize=(10,10))
	#img = ax.pcolormesh(Pthermal/pW, Pg, np.log(NEPtls/aW), cmap="viridis")
	img = ax.pcolormesh(Pthermal/pW, Pg, NEFtls, cmap="viridis")
	ax.set_xlabel('Popt [pW]')
	ax.set_ylabel('Pread [dBm]')
	ax.grid()
	ax.set_title("Tbath = %1.1f mK"%(Tbath*1e3))
	fig.colorbar(img, label='NEF TLS(1 Hz) [Hz/rtHz]')
	plt.show()

	exit()
	# Let's make the TLS NEP histogram just for kicks
	fig, ax = plt.subplots(figsize=(10,10))
	ax.hist((NEPtls/aW).flatten(), bins=20, color='k',\
		density=True, histtype='stepfilled')
	ax.set_xlabel('NEP TLS [aW/rtHz]')
	ax.grid()
	plt.show()
	exit()



	#C = 1./(wr**2*L)
	#C_branch = 0.11 * pF
	#N = 10
	#C_load = C_branch * N
	#Zin = Z0
	#Zout = 1/(1./Z0 + 1j*wr*C_load)
	#CC = np.sqrt((2*C)/(Qi*Z0*wr))
	#y = wr * CC * Z0
	#Gprime = (wr*CC*Zin*y/(Zin + Zout)) - (1j*wr*CC*Zin**2*y**2)/(Zin + Zout)**2
	#dQe = Gprime/(wr*C)
	#Qe = 1/dQe
	#Qc = 1./np.real(dQe)
	#print (Qc)
	#phi_c = np.arctan2(dQe.imag, dQe.real)
	#print (phi_c)
	exit()
	#print (Qc)
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

	print ("Total Coupling Capacitance Needed ", CC/pF)
	CC *= 2 #Needed because we have 2 coupling capacitors in series
	# I want to pick the largest capacitor as the base IDC and then design blade
	# cells to create the other capacitors
	Cmain = np.max(C)
	print (C/pF)

	# Want to determine the number of fingers for this capacitor structure.
	Nfingers, capfrac = cap.getnumfingers(Cmain)
	print ("Expected number of fingers is ", Nfingers)
	print ("Expected number of fractional fingers is ", capfrac)

	cap.set_dimensions(trace_width, gap_width,
			finger_length, finger_gap, Nfingers, contact_width)
	C_actual = cap.capacitance()/pF
	print ("The expected capacitance of the structure is %1.3f pF"%C_actual)
	print ("Width/height: (%1.3f, %1.3f) um"%(cap.width, cap.height))

	cap.layer = 0
	cap.cellname = 'Cap'
	capcell = cap.draw()
	cap_cell = gdspy.Cell('Cap_%1.0fMHz'%fr[0])
	cref = gdspy.CellReference(capcell, rotation=90)
	cap_cell.add(cref)

	print ("Target coupling capacitance", CC/pF)
	print ("Actual coupling capacitance", coupcap.capacitance()/pF)

	# Quick comparison between the design here and the predictions from sonnet
	# simulations of the capacitor structures.
	C_sim, L_par = get_simulated_LC(Nfingers - np.array(coarse_blades))
	L_total = L + L_par*nH
	C_sim *= pF
	wr = 1./np.sqrt(L_total * C_sim)
	fr = wr/2/pi/MHz
	print (fr)
	C_branch = 0.11 * pF
	C_load = C_branch * N
	Zin = Z0
	Zout = 1/(1./Z0 + 1j*wr*C_load)
	Cc = np.sqrt((2*C)/(Qi*Z0*wr))
	y = wr * Cc * Z0
	Gprime = (wr*Cc*Zin*y/(Zin + Zout)) - (1j*wr*Cc*Zin**2*y**2)/(Zin + Zout)**2
	dQe = Gprime/(wr*C_sim)
	Qe = 1/dQe
	Qc = 1./np.real(dQe)
	phi_c = np.arctan2(dQe.imag, dQe.real)
	print (Qc)
	print (phi_c*180/pi)
	#Qc = (C*pF)/(0.5*Z0*wr*(CC*pF/2)**2)
	Qr = 1./(1./Qc + 1./Qi)
	#chi_c = 4*Qr**2/Qi/Qc/np.cos(phi_c)
	#chi_c = 4*np.abs(Qe)*Qi/(Qi*np.cos(phi_c) + np.abs(Qe))**2
	chi_c = 4*Qc*Qi/(Qi + Qc)**2#/np.cos(phi_c)
	print (chi_c)

	exit()












