#!/usr/bin/env python3

import numpy as np
from math import pi
import gdspy
import sys
import matplotlib.pyplot as plt
from scipy import constants,special
from scipy.special import iv, kn
from scipy.constants import h, k
import astropy.units as u
import patches
from openpyxl import Workbook
import pdb
import pprint
pp = pprint.PrettyPrinter(indent=4)

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

number_of_points = 32

mUnits = 1e-6
wafer_width = 101e3
wafer_len = 101e3
inv_margin = 200
edge_margin = 1000
finger_length = 502
cap_margin = 4
indboxwidth = 325.
indboxheight = 127.
coup_cap_finger_length = 50
bondpad_size = 140
island_halfwidth = 800
#Qi = 40000
Z0 = 50 * u.Ohm
gnd_box_margin = 150
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def_layers = {"Wafer_footprint":1, "Reticle_Outline":2, "GP":3, "CIF_LSNSUB":4, "CIF_ILD":5,
		"Cap_Nb_120nm":6, "Capacitor_Etch":7, "Al":8, "Au_Heater":9, "RES_PRO":10,
		"MS":11, "Au_Thermal_Sink":12, 'CIF_LSN1':13, 'XeF2':14,
		"Filler_btwn_cells":15, "Anodize":23}
layer_order = [12, 9, 6, 3, 8 ]
feed_cell_length = 100
feed_cell_width = 20
feed_main_width = 8
feed_ms_fraction = 0.2
feedres_dist = 209
mask_width = 22000
mask_length = 26000


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



def moveabove_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_top = fixed_ury
	free_bottom = free_lly

	yshift = (fixed_top - free_bottom)
	return yshift

def moveabove(object_fixed,object_free,spacing=0):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	yshift = moveabove_get_shift(box_fixed,box_free)
	object_free.translate(0,yshift-spacing)

def movebelow_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_bottom = fixed_lly
	free_top = free_ury

	yshift = (fixed_bottom - free_top)
	return yshift

def movebelow(object_fixed,object_free,spacing=0):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	yshift = movebelow_get_shift(box_fixed,box_free)
	object_free.translate(0,yshift-spacing)

def moveleft_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_left = fixed_llx
	free_right = free_urx

	xshift = (fixed_left - free_right)
	return xshift

def moveleft(object_fixed,object_free,spacing=0):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	xshift = moveleft_get_shift(box_fixed,box_free)
	object_free.translate(xshift-spacing,0)

def moveright_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_right = fixed_urx
	free_right = free_llx

	xshift = (fixed_right - free_right)
	return xshift

def moveright(object_fixed,object_free,spacing=0):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	xshift = moveright_get_shift(box_fixed,box_free)
	object_free.translate(xshift-spacing,0)

def centerx_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_mid = 0.5*(fixed_llx + fixed_urx)
	free_mid = 0.5*(free_llx + free_urx)

	xshift = (fixed_mid - free_mid)
	return xshift

def centerx(object_fixed,object_free):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	xshift = centerx_get_shift(box_fixed,box_free)
	object_free.translate(xshift,0)

def centery_get_shift(box_fixed,box_free):
	(fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
	(free_llx,free_lly),(free_urx,free_ury) = box_free

	fixed_mid = 0.5*(fixed_lly + fixed_ury)
	free_mid = 0.5*(free_lly + free_ury)

	yshift = (fixed_mid - free_mid)
	return yshift

def centery(object_fixed,object_free):
	box_fixed = object_fixed.get_bounding_box()
	box_free = object_free.get_bounding_box()
	yshift = centery_get_shift(box_fixed,box_free)
	object_free.translate(0,yshift)

def recenter(obj):
	(xmin, ymin), (xmax, ymax) = obj.get_bounding_box()
	midx = (xmax - xmin)//2
	midy = (ymax - ymin)//2
	x0, y0 = obj.origin
	dx = (x0 - midx)
	dy = (y0 - midy)

	obj.translate(dx, dy)

class InductorMeander():
	def __init__(self):
		self.tracewidth = 1.0
		self.gap = 1.0
		self.boxwidth = indboxwidth
		self.boxheight = indboxheight
		self.x0 = 0
		self.y0 = 0
		self.layer = 0
		self.cellname = 'Inductor'
	def draw(self):
		path = gdspy.Path(self.tracewidth,(self.x0,self.y0))
		specseg = {'layer':self.layer}
		specturn = {'layer':self.layer,'number_of_points':number_of_points}
		radius = 0.5*(self.gap + self.tracewidth)
		turnwidth = 0.5*self.gap + self.tracewidth
		path.segment(turnwidth,**specseg)
		seglen = self.boxwidth - 2*turnwidth
		while path.y - self.y0 + 4*(self.gap + self.tracewidth) < self.boxheight:
			path.segment(seglen,**specseg)
			path.turn(radius,'ll',**specturn)
			path.segment(seglen,**specseg)
			path.turn(radius,'rr',**specturn)
		path.segment(seglen,**specseg)
		path.turn(radius,'ll',**specturn)
		path.segment(seglen,**specseg)
		path.segment(turnwidth,**specseg)
		self.x1,self.y1 = path.x,path.y

		self.cell = gdspy.Cell(self.cellname)
		self.cell.add(path)
		return self.cell

class IDC():
	def __init__(self,capfrac):
		self.trace_width = 2.
		self.gap_width = 2.
		self.finger_length = 500.
		self.finger_gap = 2.
		self.nfinger = 64
		self.contact_width = cap_margin#self.trace_width
		self.width = self.contact_width + self.finger_gap + self.finger_length + self.contact_width
		self.height = self.nfinger*(self.gap_width+self.trace_width) + self.trace_width
		self.capfrac = capfrac
		self.cellname = 'Capacitor'
		self.layer= 6 
		self.C = self.capacitance()

	def draw(self, less_one=False):
		dx = (self.finger_length + self.gap_width )/2
		dy = self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
		self.make_fingers(less_one)
		self.make_contacts()
		self.C = self.capacitance()

		#self.left_contact.layer = self.layer
		self.left_contact.translate(-dx, -dy/2)
		self.right_contact.translate(-dx, -dy/2)
		#self.right_contact.layer = self.layer
		self.cell = gdspy.Cell(self.cellname)
		contacts = gdspy.boolean(self.left_contact, self.right_contact, 'or',
				layer=self.layer)
		#print ("adding fingers")
		Nfingers = len(self.fingers)//2
		for i, f in enumerate(self.fingers):
			f.translate(-dx, -dy/2)
			contacts = gdspy.boolean(contacts, f, 'or', max_points=Nfingers*4,
					layer=self.layer)
			#print (self.cellname, i, len(contacts.polygons))
		self.cell.add(contacts)
		#print (self.cellname, len(contacts.polygons))

		return self.cell

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
			box = gdspy.Rectangle([lower_leftx,lower_lefty],
					[upper_rightx,upper_righty], layer=self.layer)
			if i % 2 == 0:
				self.fingers.append(box)
			elif (i/2) <= nrff:
				self.fingers.append(box)

	def make_contacts(self):
		lower_leftx = -self.contact_width
		lower_lefty = 0.
		upper_rightx = lower_leftx + self.contact_width
		upper_righty = lower_lefty + self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
		self.left_contact = gdspy.Rectangle([lower_leftx,lower_lefty],
				[upper_rightx,upper_righty], layer=self.layer)
		lower_leftx = self.finger_gap + self.finger_length
		upper_rightx = lower_leftx + self.contact_width
		self.right_contact = gdspy.Rectangle([lower_leftx,lower_lefty],
				[upper_rightx,upper_righty], layer=self.layer)

def Rectangle(x0,y0,x1,y1):
	r = gdspy.Rectangle((x0,y0),(x1,y1))
	def get_bounding_box():
		pt0,pt1,pt2,pt3 = r.points
		return pt0,pt2
	r.get_bounding_box = get_bounding_box
	return r

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
	capfrac = n - nfingers

	return nfingers, capfrac

def placeuniquecaps(caps_inv, mask, mcols, nrows, ncols):
	# mrows, mcols are the number of rows and columns on the mask.
	N = nrows * ncols
	# All the capacitors have the same size
	cap_dx, cap_dy = get_size(caps_inv[0])

	# startx = -mask_width/2 + cap_dx/2
	# starty = mask_length/2 - cap_dy/2
	startx = -mask_width/2 + cap_dx/2
	starty = (N // mcols + 1) * cap_dy + cap_dy/2

	caps_inv_ref = [None]*N
	caps_inv_ref[0] = gdspy.CellReference(caps_inv[0],\
			[startx, starty])
	for icap in range(1, N):
		cap_ref = gdspy.CellReference(caps_inv[icap])
		if icap % mcols == 0:
			movebelow(caps_inv_ref[icap - mcols], cap_ref)
			centerx(caps_inv_ref[icap - mcols], cap_ref)
		else:
			moveright(caps_inv_ref[icap - 1], cap_ref)
			centery(caps_inv_ref[icap - 1], cap_ref)
		caps_inv_ref[icap] = cap_ref

	mask.add(caps_inv_ref)
	return caps_inv_ref

def get_wafer_edges():
	edge_l = 125
	edge_w = 15e3
	hor = gdspy.Rectangle([-edge_w/2, edge_l/2], [edge_w/2, -edge_l/2],
			layer=def_layers['Cap_Nb_120nm'])
	vert = gdspy.Rectangle([-edge_l/2, edge_w/2], [edge_l/2, -edge_w/2],
			layer=def_layers['Cap_Nb_120nm'])
	hor_cell = gdspy.Cell('50umX15mm_Hline')
	vert_cell = gdspy.Cell('50umX15mm_Vline')
	hor_cell.add(hor)
	vert_cell.add(vert)

	return hor_cell, vert_cell

# Makes a rectangle on a default layer
def make_rectangle(width,length):
	return gdspy.Rectangle([-width/2, length/2], [width/2, -length/2], layer=def_layers['Cap_Nb_120nm'])

def fill_empty_space(cell, width, length):
	filler = make_rectangle(width, length)
	subcells = cell.references
	for subcell in subcells:
		refname = subcell.ref_cell.name
		print (refname)
		dx, dy = get_size(subcell)
		subrect = make_rectangle(dx, dy)
		subrect.translate(*subcell.origin)
		filler = gdspy.fast_boolean(filler, subrect,
				'xor',layer=def_layers['Filler_btwn_cells'])
		# print (subcell)
	return filler

def get_mask_lens():
	diameter = 31112
	radius = diameter/2
	lens = gdspy.Round([0,0], radius, number_of_points=2000, layer=def_layers['Anodize'])
	return gdspy.fast_boolean(lens, None, 'or', layer=def_layers['Anodize'])

def get_wafer_outline():
	o_diameter = 150000
	o_radius = o_diameter/2
	i_diameter = o_diameter - 6000
	i_radius = i_diameter/2
	o_outline= gdspy.Round([0,0], o_radius, number_of_points=2000, layer=def_layers['Wafer Outline'])
	i_outline= gdspy.Round([0,0], i_radius, number_of_points=2000, layer=def_layers['Wafer Outline']) 

	o_outline = gdspy.fast_boolean(o_outline, None, 'or', layer=def_layers['Wafer Outline'])
	i_outline = gdspy.fast_boolean(i_outline, None, 'or', layer=def_layers['Wafer Outline'])
	return o_outline, i_outline   

def get_size(cell):
	(xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
	return (xmax - xmin), (ymax - ymin)

def flip_cell(cell, layer, max_points=4000):
	cell_ref = gdspy.CellReference(cell, rotation=180)
	cell_polys = cell_ref.get_polygonsets()
	cell_polys_merged = cell_polys[0]
	cell_polys = gdspy.boolean(cell_polys,None,'or',
			max_points=max_points, layer=layer)
	cell.remove_polygons(lambda pts, lyr, datatype: lyr == layer)
	cell.add(cell_polys)
	return cell


def get_cap_tank(index, cap, coupcap, nfingers, commoncap):
	#pdb.set_trace()
	assembled_cap = gdspy.Cell('assembled_' + cap.ref_cell.name)
	Dx, Dy = get_size(cap)
	cDx, cDy = get_size(commoncap)
	dx, dy = get_size(coupcap)
	coupoffset = 10
	overlap = 2
	ysize = Dy + cDy
	xsize = Dx
	if index % 4 == 0 or index % 4 == 1:
		# Up case
		cap.translate(0, -(Dy+2*cDy-2*overlap-ysize)/2)
		commoncap.translate(0, (ysize-cDy-coupoffset)/2)
		coupcap.translate(0, (ysize-dy)/2-coupoffset)
	else:
		# Down case
		cap.translate(0, (Dy+2*cDy-2*overlap-ysize)/2)
		commoncap.translate(0, -(ysize-cDy-coupoffset)/2)
		coupcap.translate(0, -(ysize-dy)/2+coupoffset)
	assembled_cap.add(cap)
	assembled_cap.add(commoncap)
	assembled_cap.add(coupcap)

	return assembled_cap


def make_capacitor_connector():
	cap_edge_connector = gdspy.Cell('capacitor_edge_connector')
	small_dim = 2
	conn_length = 253 + small_dim/2
	conn_width = 8
	main_connector1 = gdspy.Rectangle([-conn_length/2, conn_width/2],
			[conn_length/2, -conn_width/2], layer=def_layers['Cap_Nb_120nm'])
	main_connector2 = gdspy.Rectangle([-conn_length/2, conn_width/2],
			[conn_length/2, -conn_width/2], layer=def_layers['Cap_Nb_120nm'])
	small_connector1 = gdspy.Rectangle([-small_dim, small_dim], [small_dim, -small_dim],
			layer=def_layers['Cap_Nb_120nm'])
	small_connector2 = gdspy.Rectangle([-small_dim, small_dim], [small_dim, -small_dim],
			layer=def_layers['Cap_Nb_120nm'])
	cap_edge_connector.add(main_connector1.translate(-conn_length/2-small_dim, 0))
	cap_edge_connector.add(main_connector2.translate(conn_length/2+small_dim, 0))
	cap_edge_connector.add(small_connector1.translate(-conn_length,
		-conn_width/2-small_dim))
	cap_edge_connector.add(small_connector2.translate(conn_length,
		-conn_width/2-small_dim))
	main_lib.add(cap_edge_connector)
	return cap_edge_connector

def make_common_capacitors(nfinger, capfrac):
	connector = make_capacitor_connector()
	common_idc_down = IDC(1.0)
	common_idc_down.finger_length = finger_length
	common_idc_down.nfinger = nfinger
	common_idc_down.cellname = "Capacitor_common_down"
	common_cap_down = common_idc_down.draw(less_one=True)
	common_cap_down_ref = gdspy.CellReference(common_cap_down)
	Dx, Dy = get_size(common_cap_down_ref)
	cap_edge_ref = gdspy.CellReference(connector, rotation=180)
	cap_edge_ref.translate(0, -(Dy/2 + 6))
	common_cap_down.add(cap_edge_ref)
	common_cap_down.flatten()

	common_cap_up = gdspy.Cell("Capacitor_common_up")
	common_cap_down_ref2 = gdspy.CellReference(common_cap_down, rotation=180)
	common_cap_up.add(common_cap_down_ref2)
	common_cap_up.flatten()

	return common_cap_up, common_cap_down


def generate_caps_from_freqs(nrows, ncols):
	N_res = 2*nrows*ncols
	L = 10 #nH
	fstart = 324 # MHz
	df = 3 #MHz
	fs = fstart + df * np.arange(N_res)

	capfrac = 1/(fs/fs[0])**2
	Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
	idcs = [IDC(capfrac[i]) for i in range(N_res)]

	T_target = 380e-3#K
	Qc = get_MB_Qi(T_target, fs*MHz)

	# Let's make the coupling capacitors as well
	Ccs = 2 * np.sqrt((Cs * u.F)/((Z0/2) * Qc * 2*np.pi*fs*u.MHz)).to(u.F).value
	coupcaps = [IDC(1.0) for i in range(len(Cs))]

	caps = []
	coup_Nfingers = []
	Nfingers = []
	cfracs = []
	coupcap_counter = 3
	for idc, c, cs, coupcap, freq in zip(idcs, Cs, Ccs, coupcaps, fs):
		idc.layer = def_layers['Cap_Nb_120nm']
		idc.finger_length = finger_length
		idc.cellname = "Capacitor_{0:3d}MHz".format(freq)
		nfinger, cfrac = getnumfingers(idc, c)
		Nfingers.append(nfinger)
		cfracs.append(cfrac)

		coupcap.layer = def_layers['Cap_Nb_120nm']
		coupcap.finger_length = coup_cap_finger_length - 2
		nfingers, _ = getnumfingers(coupcap, cs)
		nfingers = int(round(nfingers/coupcap_counter)*coupcap_counter)
		coup_Nfingers.append(nfingers)


	coup_Nfingers = np.array(coup_Nfingers)
	unq_nfingers, unq_indices, unq_inv = np.unique(coup_Nfingers,\
			return_index=True, return_inverse=True)
	unq_coupcaps = []
	for index in unq_indices:
		coupcaps[index].nfinger = coup_Nfingers[index]
		coupcaps[index].cellname =\
				"coupling_capacitor_{0:d}".format(coup_Nfingers[index])
		coupcap_cell = coupcaps[index].draw()
		dx, dy = get_size(coupcap_cell)
		gap = 2
		overlap = cap_margin
		xspacing = dx + finger_length + gap
		coupcap_array = gdspy.CellArray(coupcap_cell, columns=2, rows=1,
				spacing=[xspacing, 0], origin=[-xspacing/2, 0])
		capdoublet = gdspy.Cell('doublet_' + coupcap_cell.name)
		capdoublet.add(coupcap_array)
		capdoublet.flatten()
		unq_coupcaps.append(capdoublet)

	print ("I've made all the coupling capacitors!")
	Nfinger = int(round(np.max(Nfingers)/10)*10)
	print ("The maximum number of fingers is, ", Nfinger)
	common_Nfinger = np.min(Nfingers)
	common_capfrac = capfrac[-1]
	common_up, common_down = make_common_capacitors(nfinger, capfrac)
	obtained_Cs = []
	obtained_Ccs = []
	for ires, coup_index in zip(range(N_res), unq_inv):
		idcs[ires].nfinger = Nfinger - common_Nfinger
		idcs[ires].capfrac = (idcs[ires].capfrac - common_capfrac)/(1 - common_capfrac)
		unqcap = idcs[ires].draw()
		#idcs[ires].capfrac = idcs[ires].capfrac
		print ("Working on idc for resonator, ", ires)
		if ires % 4 == 0 or ires % 4 == 1:
			commoncap_ref = gdspy.CellReference(common_up)
			unqcap = flip_cell(unqcap, max_points=4*idcs[ires].nfinger,
					layer=def_layers['Cap_Nb_120nm'])
			unqcap_ref = gdspy.CellReference(unqcap)
		else:
			commoncap_ref = gdspy.CellReference(common_down)
			unqcap_ref = gdspy.CellReference(unqcap)
		caps.append(get_cap_tank(ires, unqcap_ref,
			gdspy.CellReference(unq_coupcaps[coup_index]),
			Nfinger, commoncap_ref))
		obtained_Cs.append(idcs[ires].C)
		obtained_Ccs.append(coupcaps[ires].C)

	# Make the capacitor cover cell
	length, width = 700, 3500
	capcover = gdspy.Cell('cover_over_capacitor')
	rect = gdspy.Rectangle([-length/2, -width/2], [length/2, width/2],
			layer=def_layers['Cap_Nb_120nm'])
	capcover.add(rect)

	# Make the capacitor cover cell for the dark resonators
	length, width = 2500, 1200
	capcoverd = gdspy.Cell('cover_over_capacitor_darkresonators')
	rect = gdspy.Rectangle([-length/2, -width/2], [length/2, width/2],
			layer=def_layers['Cap_Nb_120nm'])
	capcoverd.add(rect)
	# Add the cells on the mask


def generate_mask():
	print ("Generating the mask....\n\n")
	all_cells = main_lib.cell_dict
	#pp.pprint(all_cells)
	make_inverted_cells()
	filler = fill_empty_space(mask, mask_width, mask_length)
	mask.add(filler)

	print ("\n\nMask Generation Completed.\n\n")
	return mask

def cellIsPresent(cellname):
	allcellnames = main_lib.cells.keys()
	return cellname in allcellnames

def get_inv_cellname(cellname):
	if cellname.endswith('_inv'):
		inv_cell_name = cellname
	else:
		inv_cell_name = cellname + '_inv'
	return inv_cell_name

def inverter(cell, rotation=0):
	cell_name = get_inv_cellname(cell.name)
	print (cell_name)
	if cellIsPresent(cell_name): return
	inv_cell = gdspy.Cell(cell_name)
	cell = cell.flatten()
	cell_ref = gdspy.CellReference(cell, rotation=rotation)
	dx, dy = get_size(cell_ref)
	dx += 2*inv_margin
	dy += 2*inv_margin
	dx = roundto(dx, 100)
	dy = roundto(dy, 100)
	layer = cell.get_layers().pop()
	polyset = cell_ref.get_polygonsets()
	print (len(polyset))
	#polyset = gdspy.PolygonSet(polys, layer=layer)
	bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)
	new_polyset = gdspy.boolean(polyset, bbox, 'xor', max_points=4000, layer=layer)
	inv_cell.add(new_polyset)
	# main_lib.add(inv_cell)
	return inv_cell

def invert_cell(cell, rotation=0):
	layers = cell.get_layers()

	if not cell.get_dependencies():
		inverter(cell)

	icells = []
	for subcell in cell.get_dependencies():
		icells.append(invert_cell(subcell, rotation))
	return icells

def make_inverted_cells():
	allcells = main_lib.cells
	cellnames = list(allcells.keys())
	# print (cellnames)
	top = allcells['Wafer_Layout_new']
	for cell in top.get_dependencies():
		if cell.name == 'WaferOutline': continue
		if not cell.name.endswith('_inv'):
			invert_cell(cell)
	# Also make inverted cells of all the capacitor cells
	for cellname in cellnames:
		if cellname.startswith('Capacitor') or cellname.startswith('doublet')\
				or cellname.startswith('cover'):
			print (cellname)
			invert_cell(allcells[cellname])

def make_inverted_cell(cellref, mask_components):
	cellname = cellref.ref_cell.name
	invcellname = get_inv_cellname(cellname)
	if invcellname in mask_components:
		invcellref = gdspy.CellReference(invcellname)
	#elif cellname in mask_components:
	#	invcellref = gdspy.CellReference(cellname)
	else:
		try:
			invcell = main_lib.cell_dict[invcellname]
			invcellref = gdspy.CellReference(invcell)
		except KeyError:
			invcell = gdspy.Cell(invcellname)
			subcells = cellref.ref_cell.elements
			for a_cell in subcells:
				a_cellref = make_inverted_cell(a_cell, mask_components)
				if type(a_cell) == gdspy.CellArray:
					a_cellarr = gdspy.CellArray(a_cellref.ref_cell, columns=a_cell.columns, rows=a_cell.rows,\
						spacing=a_cell.spacing,origin=a_cell.origin)
					invcell.add(a_cellarr)
				else:
					a_cellref.origin = a_cell.origin
					invcell.add(a_cellref)
			invcellref = gdspy.CellReference(invcell)
	invcellref.origin = cellref.origin


	return invcellref

def generate_inverted_overlay(wafer, mask):
	mask_components = {x.ref_cell.name for x in mask.elements if type(x) in patches.allowed_element_types}

	# invcomponent = make_inverted_cell(component, mask_components)
	invwafer = gdspy.Cell('Global_Overlay_Inverted')
	gcomponents = wafer.elements
	for component in gcomponents:
		# print (component)
		ctype = type(component)
		if ctype not in patches.allowed_element_types: continue
		invcomponent = make_inverted_cell(component, mask_components)
 
		if ctype == gdspy.CellArray:
			# continue
			invcomponent = gdspy.CellArray(invcomponent.ref_cell, columns=component.columns, rows=component.rows,\
			spacing=component.spacing)
		invcomponent.origin = component.origin
		invwafer.add(invcomponent)


	return invwafer

def check_cell_for_overlaps(cell):
	maxdepth = 100
	all_polygons = cell.get_polygons(depth=0)
	all_polygons = [gdspy.Polygon(x) for x in all_polygons]
	N = len(all_polygons)
	print (N)
	counter = 0
	for i in range(N):
		for j in range(i + 1, N):
			overlap = gdspy.fast_boolean(all_polygons[i], all_polygons[j], 'and')
			if overlap is not None:
				return True
			# print (counter)
			counter += 1

	return False

def collapse_cell(cell, layer):
	new_cellname = cell.name + '_singlelayer'
	new_cell = cell.copy(name=new_cellname,deep_copy=True)
	new_cell = new_cell.flatten(single_layer=layer)
	return new_cell


def makechipoutline(width, length, outlinename):
	specs = {'layer':10}
	outline = gdspy.Path(2, (0,0))
	#outline.layers = [0]
	outline.segment(width, '+x', layer=def_layers['Filler_btwn_cells'])
	outline.segment(length, '+y', layer=def_layers['Filler_btwn_cells'])
	outline.segment(width, '-x', layer=def_layers['Filler_btwn_cells'])
	outline.segment(length, '-y', layer=def_layers['Filler_btwn_cells'])
	outline.translate(-width/2, -length/2)
	outlinecell = gdspy.Cell(outlinename)
	outlinecell.add(outline)
	return outlinecell

def get_inverted_cells():
	cells = main_lib.cells
	cellnames = cells.keys()

	return [cells[name] for name in cellnames if name.endswith('_inv')]

def get_cell_area(cell):
	dx, dy = get_size(cell)
	return dx * dy

def generate_main_mask():
	print ("Generating the main mask....\n\n")
	all_cells = main_lib.cells

	mask = all_cells['TKIDModule_Main_Reticle']
	filler = fill_empty_space(mask, mask_width, mask_length)
	mask.add(filler)

	print ("\n\nCapacitor Mask Generation Completed.\n\n")
	return mask



def generate_capacitor_mask():
	print ("Generating the capacitor mask....\n\n")
	all_cells = main_lib.cells
	mask = all_cells['TKIDModule_Capacitor_Reticle']
	return mask

	mask = gdspy.Cell('TKIDModule_Capacitor_Reticle')
	maskoutline = makechipoutline(mask_width, mask_length, 'MaskOutline')
	outline_width = 2
	# mask.add(gdspy.CellReference(maskoutline))

	default_spacing=200 # A suitable spacing between the cells on the reticle
	intercap_spacing=50

	# Placement of the unique capacitors on the mask. I'll fit a set of 30 caps above the common pixel
	# another 30 below and 2 flanking each side of the common pixel
	inv_cell_list = get_inverted_cells()
	not_yet = set(inv_cell_list)
	caps_inv = [gdspy.CellReference(cell) for cell in inv_cell_list \
			if cell.name.startswith('Capacitor') and cell.name.find('common')<0]
	caps_inv.sort(key=lambda x:x.ref_cell.name)
	total_mask_area = mask_length * mask_width
	total_area_needed = np.sum(list(map(lambda x: get_cell_area(x),\
			caps_inv)))
	mrows = 8
	mcols = 16

	u_dx, u_dy = get_size(caps_inv[0])

	u_dx += intercap_spacing
	u_dy += intercap_spacing

	# Place all the caps above the common pixel
	for index, cap in enumerate(caps_inv):
		irow, icol = index // mcols, index % mcols
		y_displacement = u_dy * (2*irow - mrows + 1)/2
		x_displacement = u_dx * (2*icol - mcols + 1)/2
		cap.translate(x_displacement , y_displacement)
		mask.add(cap)

	icommoncap_up = all_cells['Capacitor_common_up_inv']
	icommoncap_up_ref = gdspy.CellReference(icommoncap_up)
	moveright(cap, icommoncap_up_ref, spacing=-intercap_spacing)
	o_dx, o_dy = icommoncap_up_ref.origin
	u_dx, u_dy = get_size(icommoncap_up_ref)
	icommoncap_up_ref.translate(0, -o_dy - intercap_spacing/2 - u_dy/2)


	icommoncap_down = all_cells['Capacitor_common_down_inv']
	icommoncap_down_ref = gdspy.CellReference(icommoncap_down)
	moveabove(icommoncap_up_ref, icommoncap_down_ref, spacing=-intercap_spacing)
	icommoncap_down_ref.translate(o_dx, 0)

	mask.add(icommoncap_up_ref)
	mask.add(icommoncap_down_ref)

	capcover = all_cells['cover_over_capacitor_inv']
	capcover_ref = gdspy.CellReference(capcover)
	c_dx, c_dy = get_size(capcover_ref)
	capcover_ref.translate(o_dx + c_dx/2 - u_dx/2, 0)
	moveabove(icommoncap_down_ref, capcover_ref, spacing=-intercap_spacing)
	mask.add(capcover_ref)
	capcoverd = all_cells['cover_over_capacitor_darkresonators_inv']
	capcoverd_ref = gdspy.CellReference(capcoverd)
	c_dx, c_dy = get_size(capcoverd_ref)
	capcoverd_ref.translate(o_dx + c_dx/2 - u_dx/2, 0)
	moveabove(capcover_ref, capcoverd_ref, spacing=-intercap_spacing)
	mask.add(capcoverd_ref)


	coupcaps_inv = [gdspy.CellReference(cell) for cell in inv_cell_list \
			if cell.name.startswith('doublet')]
	coupcaps_inv.sort(key=lambda x:x.ref_cell.name)

	cu_dx, cu_dy = get_size(caps_inv[0])
	u_dx, u_dy = get_size(coupcaps_inv[0])

	u_dx += intercap_spacing
	u_dy += intercap_spacing
	cmcols = mcols
	cmrows = mrows
	mrows = 15
	mcols = 2
	starting_x = -(cmcols + 1)/2*(cu_dx + intercap_spacing) - intercap_spacing
	starting_y= -(mrows+1)/2*u_dy
	#print (starting_x, starting_y)
	#exit()


	# Place all the caps above the common pixel
	for index, coupcap in enumerate(coupcaps_inv):
		icol, irow = index // mrows, index % mrows
		#print (icol, irow)
		if icol == 0:
			if irow == 0:
				coupcap.translate(starting_x, starting_y)
			else:
				coupcap.translate(starting_x, 0)
				moveabove(coupcaps_inv[index-1], coupcap, spacing=-intercap_spacing)
		elif icol == 1:
			if irow == 0:
				coupcap.translate(starting_x - u_dx, starting_y)
			else:
				coupcap.translate(starting_x - u_dx, 0)
				moveabove(coupcaps_inv[index-1], coupcap, spacing=-intercap_spacing)
		#print (coupcap.origin)

		mask.add(coupcap)
	#exit()
	filler = fill_empty_space(mask, mask_width, mask_length)
	mask.add(filler)


	print ("\n\nCapacitor Mask Generation Completed.\n\n")
	return mask



def main():
	# Wafer organization all dimensions in microns
	nrows = 8
	ncols = 8
	###########################################################################
	#																		 #
	#				   CAPACITOR GENERATION.							#
	#																		 #
	###########################################################################
	final_fn = 'TKID_Module_20200423_AW.gds'
	main_lib.read_gds(final_fn)
	cell_dict = list(main_lib.cells.keys())


	#for cell in cell_dict:
	#	if cell.startswith('TKIDModule_Capacitor_Reticle'): main_lib.remove(cell)
	#	if cell.startswith('MaskOutline'): main_lib.remove(cell)
	#	if cell.startswith('assembled'): main_lib.remove(cell)
	#	if cell.startswith('doublet'): main_lib.remove(cell)
	#	if cell.startswith('Capacitor'): main_lib.remove(cell)
	#	if cell.startswith('capacitor_edge_connector'): main_lib.remove(cell)
	#	if cell.startswith('coupling'): main_lib.remove(cell)
	#	if cell.startswith('cover'): main_lib.remove(cell)

	#generate_caps_from_freqs(nrows, ncols)
	#make_inverted_cells()
	main_mask = generate_main_mask()
	#cap_mask = generate_capacitor_mask()
	#canon_lens = get_mask_lens()
	#cap_mask.add(canon_lens)
	#main_mask.add(canon_lens)
	#main_lib.add(cap_mask)
	main_lib.write_gds(final_fn)#,unit=1e-6,precision=1e-9)



if __name__=='__main__':
	main()


