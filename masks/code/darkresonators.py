#!/usr/bin/env python3

import numpy as np
import gdspy
import matplotlib.pyplot as plt
from scipy import constants,special
import astropy.units as u
number_of_points = 32

mUnits = 1e-6
wafer_width = 101e3
wafer_len = 101e3
edge_margin = 500
finger_length = 800
indboxwidth = 200.
indboxheight = 280.

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
    self.finger_length = 300.
    self.finger_gap = 2.
    self.nfinger = 64
    self.contact_width = self.trace_width
    self.width = self.contact_width + self.finger_gap + self.finger_length + self.contact_width
    self.height = self.nfinger*(self.gap_width+self.trace_width) + self.trace_width
    self.capfrac = capfrac
    self.cellname = 'Capacitor'
    self.layer= 1

  def draw(self, less_one=False):
    self.make_fingers(less_one)
    self.make_contacts()
    self.C = self.capacitance()

    self.left_contact.layer = self.layer
    self.right_contact.layer = self.layer
    self.cell = gdspy.Cell(self.cellname)
    for f in self.fingers:
      f.layer = self.layer
      self.cell.add(f)
    self.cell.add(self.left_contact)
    self.cell.add(self.right_contact)

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
    N = self.nfinger - 1
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
    nrff = int(nrfx)    # Number of right fingers to leave fully intact
    nrfr = nrfx - nrff    # How much to truncate the last right finger

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

def Rectangle(x0,y0,x1,y1):
  r = gdspy.Rectangle((x0,y0),(x1,y1))
  def get_bounding_box():
    pt0,pt1,pt2,pt3 = r.points
    return pt0,pt2
  r.get_bounding_box = get_bounding_box
  return r

def makechipoutline(wafer_w, wafer_l):
  specs = {'layer':2}
  outline = gdspy.Path(1, (0,0))
  #outline.layers = [0]
  outline.segment(wafer_w, '+x', **specs)
  outline.segment(wafer_l, '+y', **specs)
  outline.segment(wafer_w, '-x', **specs)
  outline.segment(wafer_l, '-y', **specs)
  outlinecell = gdspy.Cell('WaferOutline')
  outlinecell.add(outline)
  return outlinecell

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

# returns the number x rounded to the nearest ten
rounded = lambda x: int(x/10)*10

def makerestankmodel(N_fingers):
  model_cap = IDC(1.0)
  min_n = np.min(N_fingers)
  model_cap.nfinger = rounded(min_n) if min_n >= 10 else min_n
  model_cap.finger_length = finger_length
  model_cap.cellname = 'Model_Res_Tank'
  return model_cap

get_size = lambda x: (x[1,0] - x[0,0], x[1,1]-x[0,1])


def get_cap_array(cap, model, model_fingers):
  N = cap.nfinger // model_fingers
  res = cap.nfinger - (N * model_fingers)
  # Make the array as a row of the ref cells lined up
  dx, dy = get_size(model.get_bounding_box())
  array = gdspy.CellArray(model, 1, N, (dx , dy - cap.gap_width), rotation=-90)
  return array


def get_cap_position(cap, index, N_res, nrows, ncols, w_spacing, l_spacing,\
    ind_start):
  dx, dy = get_size(cap.get_bounding_box())
  index += 1
  row = (N_res - index) // nrows
  col = (N_res - index) - (row * ncols)
  space = 6 #small gap between cap and ind
  # First let's locate the inductor associated with the capacitor
  x_ind = ind_start[0] + row * (w_spacing)
  y_ind = ind_start[1] - col * (l_spacing)

  x = x_ind - space - dx
  y = y_ind + (dy - indboxheight)/2

  return (x,y)

def make_connector(len_conn):
  connector = gdspy.Path(2, (0,0))
  connector.segment(9, '+x', layer=1)
  connector.turn(2, 'r', layer=1)
  connector.segment(len_conn, '-y', layer=1)
  conn = gdspy.Cell('connectors')
  conn.add(connector)
  return conn

def add_connector_array(ncols, nrows, indspacing, ind_start):
  len_conn = (finger_length - indboxheight + 5 )/2
  connector = make_connector(len_conn)
  origin = (ind_start[0] - 10, ind_start[0] + indboxheight+ len_conn )
  return gdspy.CellArray(connector, ncols, nrows, indspacing, origin)

def main():
  # Wafer organization all dimensions in microns
  nrows = 36
  ncols = 36
  N_res = nrows * ncols
  width_spacing = (wafer_width - 2*edge_margin)/ncols
  len_spacing = (wafer_len - 2*edge_margin)/nrows
  wafer = gdspy.Cell('6_inch_wafer')
  chip_outline = makechipoutline(wafer_width, wafer_len)
  wafer.add(gdspy.CellReference(chip_outline, (0,0)))

  # Make the array of inductors first since it is simple enough
  ind = InductorMeander()
  indcell = ind.draw()
  indspacing = [width_spacing, len_spacing]
  indarray = gdspy.CellArray(indcell,ncols, nrows, indspacing)
  indarray.translate(edge_margin, edge_margin)
  indbbox = indarray.get_bounding_box()
  ind_start = (indbbox[0,0], indbbox[1,1])
  # Let's think about the frequency spacing
  df = 2 #MHz
  fstart = 300 #MHz
  fs = (fstart + df * np.arange(N_res))
  L = 22 #nH
  Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
  caps = [IDC(1.0) for i in range(len(Cs))]

  N_fingers = [0]*N_res
  for i in range(N_res):
    caps[i].finger_length = finger_length
    nfingers, capfrac = getnumfingers(caps[i], Cs[i])
    caps[i].nfinger = nfingers
    caps[i].capfrac = capfrac
    N_fingers[i] = nfingers

  model_cap = makerestankmodel(N_fingers)
  model_cap_cell = model_cap.draw(less_one=True)

  cap_cells = []
  for i in reversed(range(N_res)):
    # Using the model capacitor construct the full capacitor as an array
    cap_array = get_cap_array(caps[i], model_cap_cell, model_cap.nfinger)
    origin= get_cap_position(cap_array, i, N_res, nrows, ncols,\
        width_spacing, len_spacing, ind_start)
    cap_array.origin = origin
    cap_cells.append(cap_array)
  #wafer.add(model_cap_ref)
  wafer.add(indarray)
  wafer.add(cap_cells)
  top_connectors = add_connector_array(ncols, nrows, indspacing, ind_start)
  top_connectors.translate(0, -4)

  bot_connectors = gdspy.copy(top_connectors, 0, -indboxheight)
  bot_connectors.x_reflection=True
  #bot_connectors.rotation = 180
  ind_dx, ind_dy = get_size(indbbox)
  bot_connectors.translate(0, ind_dy - finger_length +1 )
  wafer.add(top_connectors)
  wafer.add(bot_connectors)
  gdspy.write_gds('darkres.gds',unit=1e-6,precision=1e-9)

  # Let's make a map of the number of unit caps needed along the array
  #grid = np.arange(N_res)[::-1].reshape((nrows, ncols))
  #z = np.array(N_fingers)[::-1].reshape((nrows, ncols)) // model_cap.nfinger
  #img = plt.imshow(z, cmap='jet', aspect='equal', interpolation=None)
  ##img.set_cmap('jet')
  #plt.colorbar()
  #plt.show()


def main_old():
  specturn = {'number_of_points':number_of_points}
  ind = InductorMeander()
  indcell = ind.draw()
  cap = IDC(1.0)
  cap.nfinger = 64
  capcell = cap.draw()
  print("cap.C: ",cap.C*1e12,"pF")

  filtcap = IDC(1.0)
  filtcap.nfinger = 128
  filtcap.cellname = 'FilterCapacitor'
  filtcapcell = filtcap.draw()
  print("filtcap.C: ",filtcap.C*1e12,"pF")

  filtind = InductorMeander()
  filtind.cellname = 'FilterInductor'
  filtind.boxheight = 400
  filtind.boxwidth = 140
  filtindcell = filtind.draw()

  tankcap = gdspy.CellReference(capcell,rotation=90)
  tankind = gdspy.CellReference(indcell)
  moveright(tankcap,tankind,-10)
  centery(tankcap,tankind)

  coupcap = IDC(1.0)
  coupcap.nfinger = 15
  coupcap.finger_length = 25
  coupcap.cellname = 'CouplingCap'
  coupcapcell = coupcap.draw()
  print("coupcap: ",coupcap.C*1e12,"pF")

  coupcaphigh = gdspy.CellReference(coupcapcell)
  coupcaplow = gdspy.CellReference(coupcapcell)
  moveleft(tankcap,coupcaphigh,10)
  moveleft(tankcap,coupcaplow,10)

  (x0,y0),(x1,y1) = coupcaphigh.get_bounding_box()
  (x2,y2),(x3,y3) = tankcap.get_bounding_box()
  coupcaphigh.translate(0,y3-y1)
  (x0,y0),(x1,y1) = coupcaplow.get_bounding_box()
  coupcaplow.translate(0,y2-y0)

  ccw1 = gdspy.Rectangle((x2-10,y3-2),(x2,y3))
  ccw2 = gdspy.Rectangle((x2-10,y2),(x2,y2+2))
  #gndtrace = gdspy.Rectangle((x0,y2-60),(x0+2160,y2-10))
  gndtrace = gdspy.Path(2,(x0+1,y2))
  gndtrace.segment(35,direction='-y')
  gndtrace.segment(1,direction='-x')
  gndtrace.segment(0,final_width=40,direction='+x')
  gndtrace.segment(2950)

  rotrace = gdspy.Rectangle((x0-25,y3-2),(x0,y3))

  #movebelow(coupcaphigh,coupcaplow,10)

  (cx0,cy0),(cx1,cy1) = tankcap.get_bounding_box()
  (ix0,iy0),(ix1,iy1) = tankind.get_bounding_box()

  cap2indcell = gdspy.Cell('Cap2Tankind')
  wire1 = gdspy.Path(2,(cx1,cy0+1))
  wire1.segment(9,direction='+x')
  wire1.turn(2,'l',**specturn)
  wire1.segment(iy0-cy0-3.0)
  wire2 = gdspy.Path(2,(cx1,cy1-1))
  wire2.segment(9,direction='+x')
  wire2.turn(2,'r',**specturn)
  wire2.segment(cy1-iy1-3.0)
  cap2indcell.add(wire1)
  cap2indcell.add(wire2)

  tankcap2tankind = gdspy.CellReference(cap2indcell)

  # Pair two inductors
  filtind1 = gdspy.CellReference(filtindcell,rotation=-90)
  filtind2 = gdspy.CellReference(filtindcell,rotation=90)
  movebelow(filtind1,filtind2,10)
  centerx(filtind1,filtind2)
  indpaircell = gdspy.Cell('InductorPair')
  indpaircell.add(filtind1)
  indpaircell.add(filtind2)

  indpair1 = gdspy.CellReference(indpaircell)
  centery(tankcap,indpair1)
  
  ind2indpaircell = gdspy.Cell('ind2indpair')
  (px0,py0),(px1,py1) = indpair1.get_bounding_box()
  (ix0,iy0),(ix1,iy1) = tankind.get_bounding_box()
  wire1 = gdspy.Path(2,(ix1-2.5,iy0))
  wire1.segment(iy0-py0-3.0,direction='-y')
  wire1.turn(2,'l',**specturn)
  wire1.segment(10.5)

  wire2 = gdspy.Path(2,(ix1-2.5,iy1))
  wire2.segment(py1-iy1-3.0,direction='+y')
  wire2.turn(2,'r',**specturn)
  wire2.segment(10.5)
  ind2indpaircell.add(wire1)
  ind2indpaircell.add(wire2)

  tankind2indpair1 = gdspy.CellReference(ind2indpaircell)

  filtcap1 = gdspy.CellReference(filtcapcell,rotation=90)
  moveright(indpair1,filtcap1,-10)
  centery(indpair1,filtcap1)

  filterstagecell = gdspy.Cell('FilterStage')
  filterstagecell.add(indpair1)
  filterstagecell.add(filtcap1)

  (px0,py0),(px1,py1) = indpair1.get_bounding_box()
  (cx0,cy0),(cx1,cy1) = filtcap1.get_bounding_box()
  wire1 = gdspy.Path(1,(px1-0.5,py0))
  wire1.segment(py0-cy0-3.0,direction='-y',final_width=2)
  wire1.turn(2,'l',**specturn)
  wire1.segment(8.5)
  wire2 = gdspy.Path(1,(px1-0.5,py1))
  wire2.segment(cy1-py1-3.0,direction='+y',final_width=2)
  wire2.turn(2,'r',**specturn)
  wire2.segment(8.5)
  wire3 = gdspy.Path(2,(cx1,cy1-1.))
  wire3.segment(8.5)
  wire3.turn(2,'r',**specturn)
  wire3.segment(py0-cy0-3.0,final_width=1)
  wire4 = gdspy.Path(2,(cx1,cy0+1.))
  wire4.segment(8.5)
  wire4.turn(2,'l',**specturn)
  wire4.segment(py0-cy0-3.0,final_width=1)

  filterstagecell.add(wire1)
  filterstagecell.add(wire2)
  filterstagecell.add(wire3)
  filterstagecell.add(wire4)

  filterstage1 = gdspy.CellReference(filterstagecell)
  moveright(tankind,filterstage1,-10)
  centery(tankind,filterstage1)

  filterstage2 = gdspy.CellReference(filterstagecell)
  moveright(filterstage1,filterstage2,1.5)
  centery(filterstage1,filterstage2)

  finalind = gdspy.CellReference(indpaircell)
  moveright(filterstage2,finalind,1.5)
  centery(filterstage2,finalind)

  pad = gdspy.Rectangle((0,0),(140,140))
  padcell = gdspy.Cell('WirebondPad')
  padcell.add(pad)

  pad1 = gdspy.CellReference(padcell)
  pad2 = gdspy.CellReference(padcell)
  movebelow(pad1,pad2,10)

  stubcell = gdspy.Cell('PadStub')
  stub1 = gdspy.Rectangle((0,0),(10,2))
  stubcell.add(stub1)

  padstackcell = gdspy.Cell('Padstack')
  padstackcell.add(pad1)
  padstackcell.add(pad2)
  stub1 = gdspy.CellReference(stubcell)
  stub2 = gdspy.CellReference(stubcell)
  moveabove(padstackcell,stub1,2)
  moveleft(padstackcell,stub1)
  movebelow(padstackcell,stub2,-2)
  moveleft(padstackcell,stub2)
  padstackcell.add(stub1)
  padstackcell.add(stub2)

  padstack1 = gdspy.CellReference(padstackcell)
  moveright(finalind,padstack1,0)
  centery(finalind,padstack1)

  cells = [indcell,indpaircell,filtcapcell,padcell,padstackcell]
  cells += [filterstagecell,filtindcell,cap2indcell,ind2indpaircell]
  cells += [coupcapcell]
  cells += [stubcell]

  kpupcell = gdspy.Cell('kpup')
  kpupcell.add(tankind)
  kpupcell.add(filterstage1)
  kpupcell.add(filterstage2)
  kpupcell.add(finalind)
  kpupcell.add(padstack1)
  kpupcell.add(tankcap2tankind)
  kpupcell.add(tankind2indpair1)
  kpupcell.add(coupcaphigh)
  kpupcell.add(coupcaplow)
  kpupcell.add(ccw1)
  kpupcell.add(ccw2)
  kpupcell.add(gndtrace)
  kpupcell.add(rotrace)
  cells += [kpupcell]

  #fs = np.linspace(1,1.066,33)
  fs = np.arange(200, 800, 2) # Aiming for 300 resonators with a 2 MHz spacing.
  capfracs = 1 / np.sqrt(fs)
  #capfracs = np.linspace(0.6,1,33)
  ncap = len(capfracs)

  kpuparray = gdspy.CellArray(kpupcell,1,ncap,[0,380])
  chipcell = gdspy.Cell('Chip')
  for i in range(ncap):
    cap = IDC(capfracs[i])
    cap.nfinger = 64
    cap.cellname = 'cap%03d'%i
    capcell = cap.draw()
    origin = (0,380*i)
    #kpupref = gdspy.CellReference(kpup,origin)
    capref = gdspy.CellReference(capcell,origin,rotation=90)
    chipcell.add(capref)
    #chipcell.add(kpupref)
    cells += [capcell]
  chipcell.add(kpuparray)


  chipref = gdspy.CellReference(chipcell)
  (x0,y0),(x1,y1) = chipref.get_bounding_box()

  rolinecell = gdspy.Cell('roline')
  roline = gdspy.Path(40,(x0,y0-40))
  roline.segment(y1-y0+80,direction='+y')
  rolinecell.add(roline)
  cells += [rolinecell]

  rflinecell = gdspy.Cell('signalline')
  roline = gdspy.CellReference(rolinecell)
  moveleft(kpupcell,roline)
  rflinecell.add(roline)
  pad1 = gdspy.CellReference(padcell)
  pad2 = gdspy.CellReference(padcell)
  centerx(roline,pad1)
  centerx(roline,pad2)
  pad1.translate(50,0)
  pad2.translate(50,0)
  moveabove(roline,pad1)
  movebelow(roline,pad2)
  rflinecell.add(pad1)
  rflinecell.add(pad2)
  rfline = gdspy.CellReference(rflinecell)
  chipcell.add(rfline)

  gndlinecell = gdspy.Cell('gndline')
  roline = gdspy.CellReference(rolinecell)
  moveright(kpupcell,roline)
  gndlinecell.add(roline)
  pad1 = gdspy.CellReference(padcell)
  pad2 = gdspy.CellReference(padcell)
  centerx(roline,pad1)
  centerx(roline,pad2)
  pad1.translate(-50,0)
  pad2.translate(-50,0)
  moveabove(roline,pad1)
  movebelow(roline,pad2)
  gndlinecell.add(pad1)
  gndlinecell.add(pad2)
  gndline = gdspy.CellReference(gndlinecell)
  chipcell.add(gndline)
  cells += [gndlinecell]

  chipref = gdspy.CellReference(chipcell)
  (x0,y0),(x1,y1) = chipref.get_bounding_box()
  outlinecell = gdspy.Cell('ChipOutline')
  path = gdspy.Path(1,(x0-100,y0-100))
  w = x1 - x0 + 200
  h = y1 - y0 + 200
  path.segment(w,direction='+x',layer=2)
  path.segment(h,direction='+y',layer=2)
  path.segment(w,direction='-x',layer=2)
  path.segment(h,direction='-y',layer=2)
  outlinecell.add(path)
  outline = gdspy.CellReference(outlinecell)
  chipcell.add(outline)
  cells += [chipcell,rflinecell,outlinecell]

  mstxtcell = gdspy.Cell('MSTXT')
  mstxt = gdspy.Text('MS',140)
  mstxtcell.add(mstxt)
  gndtxtcell = gdspy.Cell('GNDTXT')
  gndtxt = gdspy.Text('GND',140)
  gndtxtcell.add(gndtxt)

  labeltxtcell = gdspy.Cell('LABELTXT')
  labeltxt = gdspy.Text('BAS20180114',140)
  labeltxtcell.add(labeltxt)

  mstxt1 = gdspy.CellReference(mstxtcell)
  moveleft(chipcell,mstxt1)
  moveabove(chipcell,mstxt1)
  mstxt1.translate(470,-220)
  chipcell.add(mstxt1)
  mstxt2 = gdspy.CellReference(mstxtcell)
  moveleft(chipcell,mstxt2)
  movebelow(chipcell,mstxt2)
  mstxt2.translate(470,220)
  chipcell.add(mstxt2)

  gndtxt1 = gdspy.CellReference(gndtxtcell)
  moveright(chipcell,gndtxt1)
  moveabove(chipcell,gndtxt1)
  gndtxt1.translate(-610,-220)
  chipcell.add(gndtxt1)
  gndtxt2 = gdspy.CellReference(gndtxtcell)
  moveright(chipcell,gndtxt2)
  movebelow(chipcell,gndtxt2)
  gndtxt2.translate(-610,220)
  chipcell.add(gndtxt2)

  labeltxt1 = gdspy.CellReference(labeltxtcell)
  centerx(chipcell,labeltxt1)
  movebelow(chipcell,labeltxt1,-220)
  chipcell.add(labeltxt1)

  cells += [mstxtcell,gndtxtcell,labeltxtcell]

  gdspy.write_gds('kpup.gds',cells=cells,unit=1e-6,precision=1e-9)

if __name__=='__main__':
  main()


