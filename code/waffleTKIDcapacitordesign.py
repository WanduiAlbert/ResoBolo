"""
2019/03/14:
-----------------------------------------------------------------------------
I'm retroactively working with this code to see if I can improve the predictions
for the resonance frequencies of the waffle TKID devices.

"""
import numpy as np
from scipy import special, constants
from math import pi
import gdspy
import pprint
import astropy.units as u
from scipy.special import iv, kn
from scipy.constants import h, k

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

nH = 1e-9
pF = 1e-12
MHz = 1e6
cm = 1e-2
g = 1e-3
mJ = 1e-3
pp = pprint.PrettyPrinter(indent=5)

mUnits = 1e-6
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

def load_data(fn, nports=2, paramtype='Y'):
    dataset = np.loadtxt(fn, skiprows=9, delimiter=',')
    p = paramtype
    if nports==2:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "21", np.complex128), (p + "22", np.complex128)])
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        p12 = dataset[:, 3] + 1j*dataset[:, 4]
        p21 = dataset[:, 5] + 1j*dataset[:, 6]
        p22 = dataset[:, 7] + 1j*dataset[:, 8]
        tup = list(zip(dataset[:, 0], p11, p12, p21, p22))
        return np.array(tup, dtype=dtypes)
    elif nports==3:
        dtypes = np.dtype([("frequency", np.float64), (p + "11", np.complex128),\
            (p + "12", np.complex128), (p + "13", np.complex128), (p + "21", np.complex128),\
            (p + "22", np.complex128), (p + "23", np.complex128), (p +"31", np.complex128),\
            (p + "32", np.complex128), (p  +"33", np.complex128)])
        p11 = dataset[:, 1] + 1j*dataset[:, 2]
        p12 = dataset[:, 3] + 1j*dataset[:, 4]
        p13 = dataset[:, 5] + 1j*dataset[:, 6]
        p21 = dataset[:, 7] + 1j*dataset[:, 8]
        p22 = dataset[:, 9] + 1j*dataset[:,10]
        p23 = dataset[:,11] + 1j*dataset[:,12]
        p31 = dataset[:,13] + 1j*dataset[:,14]
        p32 = dataset[:,15] + 1j*dataset[:,16]
        p33 = dataset[:,17] + 1j*dataset[:,18]
        tup = lipt(zip(dataset[:, 0], p11, p12, p13, p21, p22, p23, p31, p32, p33))
        return np.array(tup, dtype=dtypes)

class IDC():
	def __init__(self,capfrac):
		self.trace_width = 2.
		self.gap_width = 2.
		self.finger_length = 900.
		self.finger_gap = 2.
		self.nfinger = 64
		self.contact_width = self.trace_width
		self.width = self.contact_width + self.finger_gap + \
			self.finger_length + self.contact_width
		self.height = self.nfinger*(self.gap_width+self.trace_width) + \
			self.trace_width
		self.capfrac = capfrac
		self.layer = 0
		self.C = self.capacitance()

	def set_dimensions(self, trace_width, gap_width, finger_length, finger_gap,\
			nfinger, contact_width=24):
		self.trace_width = trace_width
		self.gap_width = gap_width
		self.finger_length = finger_length
		self.finger_gap = finger_gap
		self.nfinger = nfinger
		self.contact_width = contact_width
		self.width = self.contact_width + self.finger_gap + \
			self.finger_length + self.contact_width
		self.height = self.nfinger*(self.gap_width+self.trace_width) + \
			self.trace_width
		self.C = self.capacitance()

	def capacitance2(self):
		p = mUnits*(self.finger_length - self.finger_gap)
		q = mUnits*self.height
		s = mUnits*self.gap_width
		w = mUnits*self.trace_width
		a = s + w
		e = constants.epsilon_0 * 11.7

		c = 0.0
		for n in range(1,100):
			j = special.j0((2*n-1)*pi*s/(2*a))
			c += (1.0/(2*n-1)) * j*j

		c *= p*q*(4.*e/(pi*a))
		return c

	def capacitance(self):
		e0 = constants.epsilon_0
		e1 = 1.0
		es = 11.7

		w = mUnits*self.trace_width
		g = mUnits*self.gap_width
		h = mUnits*self.finger_length
		N = int(self.nfinger*self.capfrac) - 1
		l = 2. * (w + g)
		eta = 2.*w/l
		kiinf = np.sin(np.pi*eta / 2.)
		kiinfp = np.sqrt(1 - kiinf**2)
		keinf = 2.*np.sqrt(eta) / (1 + eta)
		keinfp = np.sqrt(1 - keinf**2)

		Kkiinf = special.ellipk(kiinf)
		Kkiinfp = special.ellipk(kiinfp)
		Kkeinf = special.ellipk(keinf)
		Kkeinfp = special.ellipk(keinfp)

		Ci = e0 * (1+es) * h * Kkiinf/Kkiinfp
		Ce = e0 * (1+es) * h * Kkeinf/Kkeinfp
		Ct = (N-3) * Ci/2 + 2*(Ci*Ce)/(Ci+Ce)
		return Ct

	def getnumfingers(cap, target_C):
		e0 = constants.epsilon_0
		e1 = 1.0
		es = 11.7

		w = mUnits*cap.trace_width
		g = mUnits*cap.gap_width
		h = mUnits*cap.finger_length
		l = 2. * (w + g)
		eta = 2.*w/l
		kiinf = np.sin(np.pi*eta / 2.)
		kiinfp = np.sqrt(1 - kiinf**2)
		keinf = 2.*np.sqrt(eta) / (1 + eta)
		keinfp = np.sqrt(1 - keinf**2)

		Kkiinf = special.ellipk(kiinf)
		Kkiinfp = special.ellipk(kiinfp)
		Kkeinf = special.ellipk(keinf)
		Kkeinfp = special.ellipk(keinfp)

		Ci = e0 * (1+es) * h * Kkiinf/Kkiinfp
		Ce = e0 * (1+es) * h * Kkeinf/Kkeinfp
		n = (target_C -  2*(Ci*Ce)/(Ci+Ce)) * 2/Ci + 4
		nfingers = int(n)
		frac_fingers = n - nfingers

		return nfingers, frac_fingers

	def draw(self, less_one=False):
		dx = (self.finger_length + self.gap_width )/2
		dy = self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
		self.make_fingers(less_one)
		self.make_contacts()
		self.C = self.capacitance()

		self.left_contact.layers = [self.layer]
		self.left_contact.translate(-dx, -dy/2)
		self.right_contact.translate(-dx, -dy/2)
		self.right_contact.layers = [self.layer]
		self.cell = gdspy.Cell(self.cellname)
		for f in self.fingers:
			f.layers = [self.layer]
			f.translate(-dx, -dy/2)
			self.cell.add(f)
		self.cell.add(self.left_contact)
		self.cell.add(self.right_contact)

		return self.cell

	def make_fingers(self, less_one):
		self.fingers = []
		xcur = 0.
		ycur = 0.

		left_fingers = []
		right_fingers = []

		nrf = (self.nfinger + 1) / 2
		nrfx = nrf * self.capfrac
		#print "nrfx: ",nrfx
		nrff = int(nrfx)	# Number of right fingers to leave fully intact
		nrfr = nrfx - nrff	# How much to truncate the last right finger

		minx = xcur + self.finger_gap
		maxx = xcur + self.finger_length
		partialx = nrfr * minx + (1-nrfr) * maxx

		range_val = self.nfinger + 1
		if less_one:
			range_val = self.nfinger

		for i in range(range_val):
			if i % 2 == 0:
				lower_leftx = xcur
				upper_rightx = xcur + self.finger_length
			else:
				if i/2 == nrff:
					lower_leftx = partialx
					upper_rightx = xcur + self.finger_length + self.finger_gap
				else:
					lower_leftx = xcur+self.finger_gap
					upper_rightx = lower_leftx + self.finger_length

			assert lower_leftx >= 0.
			assert upper_rightx <= self.finger_length + self.finger_gap

			lower_lefty = ycur
			upper_righty = lower_lefty + self.trace_width

			ycur += self.trace_width + self.gap_width
			box = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])
			if i % 2 == 0:
				self.fingers.append(box)
			elif (i/2) <= nrff:
				self.fingers.append(box)

	def make_contacts(self):
		lower_leftx = -self.contact_width
		lower_lefty = 0.
		upper_rightx = lower_leftx + self.contact_width
		upper_righty = lower_lefty + self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
		self.left_contact = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])
		lower_leftx = self.finger_gap + self.finger_length
		upper_rightx = lower_leftx + self.contact_width
		self.right_contact = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])

def get_simulated_LC(Nfingers):
	C_sim = 0.05*Nfingers + 0.788
	L_sim = 6.1e-6*Nfingers**2 + 1e-3*Nfingers + 1.780
	return C_sim, L_sim

def admittance_model(x, C, L, R):
    w = 2*pi*x
    return w*C/np.sqrt((1 - w**2*L*C)**2 + (w*R*C)**2)
    #return np.log(np.abs(1j*w*C/(1 - w**2*L*C + 1j*w*R*C)))

def tline_model(x, Zc, w0, eps):
	w = 2*pi*x
	return np.abs(-1/Zc/(eps*np.cos(w/w0) + 1j*np.sin(w/w0)))

def chisq(theta, x, y):
    return np.sum((y - admittance_model(x, *theta))**2)

def get_S21(Y):
    Y0 = 1/Z0
    DY = (Y['Y11'] + Y0)*(Y['Y22'] + Y0) -Y['Y12']*Y['Y21']
    return -2 * Y['Y21']*Y0/DY

if __name__=="__main__":
	tracew = 2
	gapw = 2
	finger_length = 75
	finger_gap = 1.5
	#nfingers = [39.5*2, 42.5*2, 45.5*2, 48.5*2, 52.5*2]
	#main_nfingers = [336.5,  365.5,  399.5,  437.5,  482.5]
	nfingers = [39.5*2, 42.5*2, 45.5*2, 48.5*2]
	main_nfingers = [336.5,  365.5,  399.5,  437.5]
	N = len(nfingers)
	Cs = []
	ccs = []
	Qcs = []
	Z0 = 50
	L0 = 5.998*nH
	x = list(map(lambda x: get_simulated_LC(x-0.5), main_nfingers))
	C_cap, L_cap = list(zip(*x))
	L_cap = np.array(L_cap)*nH
	L_geom = L0 + L_cap
	C_cap = np.array(C_cap)*pF
	fr_meas = np.array([351.5, 337.4, 318.4, 305.8])
	target_T = 380e-3
	Qi = []
	for fr in fr_meas:
		Qi.append(get_MB_Qi(target_T, fr*MHz))
	Qi = np.array(Qi)
	for i in range(N):
		cmain = IDC(1.0)
		cmain.contact_width=25
		cmain.set_dimensions(tracew, gapw, 996, finger_gap,\
			main_nfingers[i])
		Cs.append(cmain.capacitance())
	Cs = np.array(Cs)
	print (Cs/pF)
	wr = 2*pi*fr_meas*MHz
	L_total = 1./(wr**2*Cs)
	L_al = 10.057*nH
	print (Cs/pF - C_cap/pF)
	#print (L_cap/nH)
	#print (L_geom/nH)
	print (L_total/nH)
	print ((L_total - L_al)/nH)
	C_load = 7 * pF
	Zin = Z0
	Zout = 1/(1./Z0 + 1j*wr*C_load)
	#Cc = np.sqrt((2*Cs)/(Qi*Z0*wr))
	Cc = np.array([0.3478, 0.3731, 0.3984, 0.4322])*pF
	CC = Cc/2
	print (CC/pF)
	#CC = np.average(Cc) # use the same coupling cap for all the resonators
	#CC = 0.1945 * pF # Using a number from actual calculations
	y = wr * CC * Z0
	Gprime = (wr*CC*Zin*y/(Zin + Zout)) - (1j*wr*CC*Zin**2*y**2)/(Zin + Zout)**2
	dQe = Gprime/(wr*Cs)
	Qe = 1/dQe
	Qc = 1./np.real(dQe)
	print (Qc)
	phi_c = np.arctan2(dQe.imag, dQe.real)
	print (phi_c)
	exit()
	#print (Qc)
	#print (phi_c*180/pi)
	#Qc = (C*pF)/(0.5*Z0*wr*(CC*pF/2)**2)
	Qr = 1./(1./Qc + 1./Qi)
	#chi_c = 4*Qr**2/Qi/Qc/np.cos(phi_c)
	#chi_c = 4*np.abs(Qe)*Qi/(Qi*np.cos(phi_c) + np.abs(Qe))**2
	chi_c = 4*Qc*Qi/(Qi + Qc)**2#/np.cos(phi_c)
	for i in range(N):
		c = IDC(1.0)
		c.contact_width=25
		c.set_dimensions(tracew, gapw, finger_length, finger_gap, nfingers[i]-1)
		ccs.append(c.capacitance()*pF)
		# Now for the main capacitor
		cmain = IDC(1.0)
		cmain.contact_width=25
		cmain.set_dimensions(tracew, gapw, 1000, finger_gap,\
			main_nfingers[i]-0.5)
		Cs.append(cmain.capacitance()*pF)

	Cs = np.array(Cs)*u.pF
	pp.pprint (Cs)
	ccs = np.array(ccs)*u.pF
	pp.pprint (ccs)
	wrs = (1/(L*Cs)**0.5).to(1/u.s)
	frs = (wrs/2/pi).to(u.MHz)
	pp.pprint (frs)
	Qcs = (8*Cs/(wrs*ccs**2*Z0)).to(1)
	pp.pprint (Qcs)
	#pp.pprint (Cs)

	# Next we do the same for the larger spacing capacitors

	tracew = 4
	gapw = 4
	finger_length = 77
	finger_gap = 3
	nfingers = [60.5, 64.5, 68.5, 74.5, 80.5]
	main_nfingers = [168.5,  184.5,  200.5,  220.5,  242.5]
	N = len(nfingers)
	uCs = []
	uccs = []
	uQcs = []
	for i in range(N):
		c = IDC(1.0)
		c.contact_width=25
		c.set_dimensions(tracew, gapw, finger_length, finger_gap, nfingers[i]-0.5)
		# print (c.width, c.height)
		uccs.append(c.capacitance()*pF)
		# Now for the main capacitor
		cmain = IDC(1.0)
		cmain.contact_width=25
		cmain.set_dimensions(tracew, gapw, 997, finger_gap,\
			main_nfingers[i]-0.5)
		# print (cmain.width, cmain.height)
		uCs.append(cmain.capacitance()*pF)

	uCs = np.array(uCs)*u.pF
	pp.pprint (uCs)
	uccs = np.array(uccs)*u.pF
	pp.pprint (uccs)
	uwrs = (1/(L*uCs)**0.5).to(1/u.s)
	ufrs = (uwrs/2/pi).to(u.MHz)
	pp.pprint (ufrs)
	uQcs = (8*uCs/(uwrs*uccs**2*Z0)).to(1)
	pp.pprint (uQcs)
	#pp.pprint (Cs)
