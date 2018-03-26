#!/usr/bin/env python3

import numpy as np
import gdspy
import matplotlib.pyplot as plt
from scipy import constants,special
import astropy.units as u
from stepper import Shot
number_of_points = 32

mUnits = 1e-6
wafer_width = 101e3
wafer_len = 101e3
inv_margin = 200
edge_margin = 1000
finger_length = 900
indboxwidth = 325.
indboxheight = 127.
coup_cap_finger_length = 100
bondpad_size = 140
island_halfwidth = 400
Qi = 40000
Z0 = 50 * u.Ohm
gnd_box_margin = 200
main_lib = gdspy.GdsLibrary('main')
gdspy.current_library = main_lib
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
    "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9}
feed_cell_length = 100
feed_cell_width = 20
feed_main_width = 8
feed_ms_fraction = 0.2
feedres_dist = 210


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
  midx = (xmax - xmin)/2
  midy = (ymax - ymin)/2
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
    self.finger_length = 900.
    self.finger_gap = 2.
    self.nfinger = 64
    self.contact_width = self.trace_width
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

    self.left_contact.layer = self.layer
    self.left_contact.translate(-dx, -dy/2)
    self.right_contact.translate(-dx, -dy/2)
    self.right_contact.layer = self.layer
    self.cell = gdspy.Cell(self.cellname)
    for f in self.fingers:
      f.layer = self.layer
      f.translate(-dx, -dy/2)
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

def makechipoutline(wafer_w, wafer_l):
  specs = {'layer':10}
  outline = gdspy.Path(1, (0,0))
  #outline.layers = [0]
  outline.segment(wafer_w, '+x', layer=def_layers['PRO1'])
  outline.segment(wafer_l, '+y', layer=def_layers['PRO1'])
  outline.segment(wafer_w, '-x', layer=def_layers['PRO1'])
  outline.segment(wafer_l, '-y', layer=def_layers['PRO1'])
  outlinecell = gdspy.Cell('WaferOutline')
  outlinecell.add(outline)
  return outlinecell


def cellIsPresent(cellname):
  allcellnames = main_lib.cell_dict.keys()
  return cellname in allcellnames

def inverter(cell, rotation=0):
  cell_name = cell.name + '_r'
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
  polys = cell_ref.get_polygons(depth=1)
  polyset = gdspy.PolygonSet(polys, layer=layer)
  bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)

  new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
  inv_cell.add(new_polyset)

def invert_cell(cell, rotation=0):
  layers = cell.get_layers()

  if len(layers) == 1:
    inverter(cell, rotation)

  for cell in cell.get_dependencies():
    invert_cell(cell, rotation)

  #cell_name = cell.name + '_r'
  #inv_cell = gdspy.Cell(cell_name)
  #cell = cell.flatten()
  #cell_ref = gdspy.CellReference(cell, rotation=rotation)
  #dx, dy = get_size(cell_ref)
  #dx += 2*inv_margin
  #dy += 2*inv_margin
  #polys = cell_ref.get_polygons(depth=1)
  #polyset = gdspy.PolygonSet(polys, layer=layer)
  #bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)

  #new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
  #inv_cell.add(new_polyset)
  #return inv_cell

def make_cap_to_ind_lines():
  layer = 8
  cell_name = 'connector'
  ms_width = 2
  len_conn = (finger_length + 3 * ms_width -ms_width)/2
  cap2ind_conn = gdspy.Cell(cell_name)
  #box = gdspy.Rectangle(layer=layer)
  conn1 = gdspy.Path(ms_width, (0,0))
  conn1.segment(5*ms_width, '+x', layer=layer)
  conn2 = gdspy.Path(ms_width, (conn1.x - ms_width/2, conn1.y + ms_width/2))
  conn2.segment(len_conn , '-y', layer=layer)
  conn3 = gdspy.Path(ms_width, (conn2.x - ms_width/2, conn2.y + ms_width/2))
  conn3.segment(island_halfwidth - 5*ms_width + ms_width + 1.25*ms_width, '+x', layer=layer)
  conn4 = gdspy.Path(ms_width, (conn3.x - ms_width/2, conn3.y - ms_width/2))
  conn4.segment(12.5, '+y', layer=layer)
  conn5 = gdspy.Path(2.5*ms_width, (conn4.x, conn4.y - ms_width/2))
  conn5.segment(51, '+y', layer=layer)

  dx = island_halfwidth/2
  dy = len_conn/2 - ms_width/2
  conn1.translate(-dx, dy)
  conn2.translate(-dx, dy)
  conn3.translate(-dx, dy)
  conn4.translate(-dx, dy)
  conn5.translate(-dx, dy)

  cap2ind_conn.add(conn1)
  cap2ind_conn.add(conn2)
  cap2ind_conn.add(conn3)
  cap2ind_conn.add(conn4)
  cap2ind_conn.add(conn5)


  return cap2ind_conn

def get_inductor():
  fn = '../resobolo_files/new_inductor.gds'
  gdsii = gdspy.GdsLibrary()
  gdsii.read_gds(fn)
  ind = gdsii.extract('Al_Inductor_Islands_right')
  polys = ind.get_polygons()
  (xmin, ymin), (xmax, ymax) = ind.get_bounding_box()
  #ind_view = gdspy.CellReference(ind, [-xmin, -ymin])
  inductor = gdspy.Cell('Al_inductor')
  for poly in polys:
    polygon = gdspy.Polygon(poly, layer=3)
    polygon = polygon.translate(-xmin, -ymin)
    polygon = polygon.translate(-(xmax - xmin)/2, -(ymax - ymin)/2)
    inductor.add(polygon)
  inductor.flatten()
  return inductor

def get_size(cell):
  (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
  return (xmax - xmin), (ymax - ymin)

def get_cap_tank(cap, coupcap):
  Dx, Dy = get_size(cap)
  tofeed = gdspy.CellReference(coupcap)
  dx, dy = get_size(tofeed)
  #print (Dx, Dy, dx, dy)
  overlap = 2
  displacement_y = -(Dy - dy)/2
  displacement_x = (Dx + dx)/2 - overlap
  tofeed.origin = [0, 0]
  tofeed.translate(displacement_x, displacement_y)
  cap.add(tofeed)
  tognd = gdspy.CellReference(coupcap)
  tognd.translate(-displacement_x, displacement_y)
  cap.add(tognd)
  cap.flatten()
  return cap

def rot(n, i, j, ri, rj):
  if (rj == 0):
    if (ri == 1):
      i = (n-1) - i
      j = (n-1) - j
    i, j = j, i
  return (i, j)

def double_range(start, stop):
  while start < stop:
    yield start
    start <<= 1

# Maps from a position in an array d, into (i,j) coordinates of an n by n
# square. Assumes we are using a Hilbert curve to fill out the space
# d in {0, n**2 - 1}
def d2ij(d, n):
  i, j = 0, 0
  t = d
  for s in double_range(1,n):
    ri = 1 & (t//2)
    rj = 1 & (t ^ ri)
    i,j = rot(s, i, j, ri,rj)
    i += s*ri
    j += s*rj
    t //= 4
  return (i, j)

uint = np.frompyfunc(int, 1,1)

def get_feedline():
  n =63
  s_gndsub = gdspy.Cell('unit_GP_sub')

  dx = n * feed_cell_length
  dy = feed_main_width
  mainline = gdspy.Rectangle([-dx/2, dy/2],\
      [dx/2, -dy/2], layer=def_layers['400nm_NbWiring'])
  main_cell = gdspy.Cell('feedline_main')
  main_cell.add(mainline)

  dy = feed_cell_width
  ild = gdspy.Rectangle([-dx/2, dy/2],\
      [dx/2, -dy/2], layer=def_layers['ILD'])
  ild_cell = gdspy.Cell('feedline_ILD')
  ild_cell.add(ild)

  # To make each cell more symmetric, I'll divide the ground sub region into 2
  sub_len = (1 - feed_ms_fraction)*feed_cell_length/2
  sub1 = gdspy.Rectangle([-sub_len/2, feed_cell_width/2],\
      [sub_len/2, -feed_cell_width/2], layer=def_layers['LSNSUB'])
  sub2 = gdspy.Rectangle([-sub_len/2, feed_cell_width/2],\
      [sub_len/2, -feed_cell_width/2], layer=def_layers['LSNSUB'])
  sdx = (sub_len + feed_ms_fraction * feed_cell_length)/2
  sub1.translate(-sdx, 0)
  sub2.translate(sdx, 0)
  s_gndsub.add(sub1)
  s_gndsub.add(sub2)
  gndsub_arr = gdspy.CellArray(s_gndsub, n, 1, [feed_cell_length,0])
  gndsub_arr.translate(-(dx - feed_cell_length)/2, 0)
  gndsub_cell = gdspy.Cell('feedline_GP_sub')
  gndsub_cell.add(gndsub_arr)
  gndsub_cell.flatten()

  feedcell = gdspy.Cell('MainFeedline')
  feedcell.add(gdspy.CellReference(main_cell))
  feedcell.add(gdspy.CellReference(ild_cell))
  feedcell.add(gdspy.CellReference(gndsub_cell))
  # feedcell.add(gdspy.CellReference(main_cell))
  # feedcell.add(gdspy.CellReference(ild_cell))
  # feedcell.add(gdspy.CellReference(gndsub_cell))


  # feedline = gdspy.Cell('large_feedline_sec')
  # feedline.add(feedarr)
  # # feedline.flatten
  # print (get_size(feedarr))

  return feedcell

# Rounds numbers to the nearest base
def roundto(num, base):
  return ((num // base) + ((num % base) > (base//2))) * base

# cres - common part of resonator
# ures - unique part of resonator
def get_resonator_GPsub(cres_bbox, ures_size):
  u_dx, u_dy = ures_size
  (xmin, ymin), (xmax, ymax) = cres_bbox
  c_dx, c_dy = (xmax - xmin, ymax - ymin)
  ymargin =  200
  xmargin = 300
  overlap = 2
  dx = c_dx + u_dx + 2 * xmargin - overlap
  dy = u_dy + 2 * ymargin

  # I want the gnd sub region to be a rounded to nearest 100
  dx = roundto(dx, 100)
  dy = roundto(dy, 100)
  # The margins have changed. Adjust them accordingly
  xmargin = (dx + overlap - c_dx - u_dx)//2
  excess = (dx + overlap - c_dx - u_dx) % 2
  ymargin = (dy - u_dy)//2
  #offset = dx/2 - (xmax + xmargin + overlap/2)
  gnd = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=def_layers['LSNSUB'])
  #gnd.translate(-offset, 0)

  c_xoffset = dx/2 - (xmargin + excess + c_dx//2 + c_dx % 2)
  delta = c_dx//2 - (overlap + c_xoffset)
  u_xoffset = delta + u_dx//2 + u_dx % 1

  # Need to make a small GND plane sub region for connecting the coupling caps to the feedline
  sdx = 20
  sdy = feedres_dist
  small = gdspy.Rectangle([-sdx/2, sdy/2], [sdx/2, -sdy/2], layer=def_layers['LSNSUB'])
  small_xoffset = delta + sdx/2 - overlap + sdx/2 + feed_cell_length/2
  small.translate(-small_xoffset, (dy + sdy)/2) #position relative to common arr

  gndsub = gdspy.Cell('reso_GP_sub')
  gndsub.add(gnd)
  gndsub.add(small)
  return gndsub, c_xoffset, u_xoffset, small_xoffset

def make_vertical_feed(feed_cell):
  vfeed = gdspy.CellReference(feed_cell, rotation=90)
  vdx, vdy = get_size(vfeed)
  vfeed.translate(0, vdy/2)
  vfeed2 = gdspy.CellReference(feed_cell, rotation=90)
  movebelow(vfeed, vfeed2)
  l = 20
  s = 8
  # Small adjustments to make the proper connections between
  # horizontal and vertical feedlines
  feed = gdspy.Rectangle([-s/2,s/2], [s/2, -s/2],\
      layer=def_layers['400nm_NbWiring'])
  feed_corner = gdspy.Cell('feed_corner')
  feed_corner.add(feed)
  feed_box = gdspy.CellReference(feed_corner)
  moveabove(vfeed, feed_box, spacing=s/2)
  ild = gdspy.Rectangle([-l/2, l/2], [l/2, -l/2],\
      layer=def_layers['ILD'])
  ild_corner = gdspy.Cell('ILD_corner')
  ild_corner.add(ild)
  ild_box = gdspy.CellReference(ild_corner)
  moveabove(vfeed, ild_box, spacing=l/2)
  gndsub = gdspy.Rectangle([-l/2, l/2], [l/2, -l/2],\
      layer=def_layers['LSNSUB'])
  gndsub_corner = gdspy.Cell('gndsub_corner')
  gndsub_corner.add(gndsub)
  gndsub_box = gdspy.CellReference(gndsub_corner)
  moveabove(vfeed, gndsub_box, spacing=l/2)
  # Same corrections at the bottom side
  feed_box_l = gdspy.CellReference(feed_corner)
  movebelow(vfeed2, feed_box_l, spacing=-s/2)
  ild_box_l = gdspy.CellReference(ild_corner)
  movebelow(vfeed2, ild_box_l, spacing=-l/2)
  gndsub_box_l = gdspy.CellReference(gndsub_corner)
  movebelow(vfeed2, gndsub_box_l, spacing=-l/2)
  vertfeedcell = gdspy.Cell('MainFeedline_vert')
  vertfeedcell.add(vfeed)
  vertfeedcell.add(vfeed2)
  vertfeedcell.add(feed_box)
  vertfeedcell.add(ild_box)
  vertfeedcell.add(gndsub_box)
  vertfeedcell.add(feed_box_l)
  vertfeedcell.add(ild_box_l)
  vertfeedcell.add(gndsub_box_l)

  return vertfeedcell

def get_coupling_conns(g_dy, u_dy):
  ymargin = (g_dy - u_dy)//2
  w = 8
  l2gnd = ymargin + 2
  l2feed = l2gnd + feedres_dist
  l2gnd += 4 # A little extra

  c2gnd = gdspy.Rectangle([-w/2,l2gnd/2], [w/2, -l2gnd/2],\
      layer=def_layers['400nm_NbWiring'])
  c2feed = gdspy.Rectangle([-w/2,l2feed/2], [w/2, -l2feed/2],\
      layer=def_layers['400nm_NbWiring'])
  cap2feed = gdspy.Cell('cap_to_feed')
  cap2feed.add(c2feed)

  cap2gnd = gdspy.Cell('cap_to_gnd')
  cap2gnd.add(c2gnd)

  return cap2feed, cap2gnd

def get_common_resonator(ind, common_cap):
  common_resonator = gdspy.Cell('common_resonator')

  ms_width = 2
  connector = make_cap_to_ind_lines()
  conn_dx, conn_dy = get_size(connector)
  dy = conn_dy/2  + ms_width/2
  top_conn_ref = gdspy.CellReference(connector, [0, dy] )
  bot_conn_ref = gdspy.CellReference(connector, [0, -dy], x_reflection=True)

  resonator_connector = gdspy.Cell('Cap_to_Ind_lines')
  resonator_connector.add(top_conn_ref)
  resonator_connector.add(bot_conn_ref)
  resonator_connector.flatten()
  resonator_connector_ref = gdspy.CellReference(resonator_connector)
  ind_ref = gdspy.CellReference(ind)
  cap_ref = gdspy.CellReference(common_cap, rotation = 90)

  conn_dx, conn_dy = get_size(resonator_connector)
  ind_dx, ind_dy = get_size(ind)
  cap_dx, cap_dy = get_size(cap_ref)

  capconn_overlap = 0
  indconn_overlap = 4
  dx = conn_dx + ind_dx + cap_dx - indconn_overlap - capconn_overlap

  # Currently the center is at the connector.
  offset = dx//2 - (cap_dx + (conn_dx- indconn_overlap)/2 - capconn_overlap )
  resonator_connector_ref.translate(-offset,0)
  common_resonator.add(resonator_connector_ref)
  #ind_origin = [island_halfwidth/2  + indboxwidth/2, 0]
  #ind_ref = gdspy.CellReference(ind, origin = ind_origin)
  moveright(resonator_connector_ref, ind_ref, spacing=indconn_overlap)
  common_resonator.add(ind_ref)
  moveleft(resonator_connector_ref, cap_ref, spacing=capconn_overlap)
  common_resonator.add(cap_ref)

  return common_resonator

def get_bondpads():
  size = 400
  feedpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['400nm_NbWiring'])
  feedpad_cell = gdspy.Cell('MSfeed_bondpad')
  feedpad_cell.add(feedpad)

  ildpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['ILD'])
  ildpad_cell = gdspy.Cell('ILDfeed_bondpad')
  ildpad_cell.add(ildpad)

  gndpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['120nm_NbWiring'])
  gndpad_cell = gdspy.Cell('GNDfeed_bondpad')
  gndpad_cell.add(gndpad)

  pad_cell = gdspy.Cell('Feed_Bonding_Pad')
  pad_cell.add(gdspy.CellReference(feedpad_cell))
  pad_cell.add(gdspy.CellReference(ildpad_cell))
  pad_cell.add(gdspy.CellReference(gndpad_cell))
  return pad_cell

def make_inverted_cells():
  allcells = main_lib.cell_dict
  cellnames = allcells.keys()
  top = allcells['Global_Overlay']
  for cell in top.get_dependencies():
    inv_cellname = cell.name + '_r'
    if inv_cellname not in cellnames:
      if inv_cellname.startswith('Capacitor'):
        invert_cell(cell, rotation=90)
      else:
        invert_cell(cell)



def get_inverted_cells():
  cellnames = main_lib.cell_dict.keys()

  return [name for name in cellnames if name.endswith('_r') and name !=\
      'WaferOutline_r']

def get_cell_area(cell):
  dx, dy = get_size(cell)
  return dx * dy

def getunitfeedline():
  dx = feed_cell_length
  dy = feed_main_width
  mainline = gdspy.Rectangle([-dx/2, dy/2],\
      [dx/2, -dy/2], layer=def_layers['400nm_NbWiring'])
  main_cell = gdspy.Cell('unit_main_feed')
  main_cell.add(mainline)

  dy = feed_cell_width
  ild = gdspy.Rectangle([-dx/2, dy/2],\
      [dx/2, -dy/2], layer=def_layers['ILD'])
  ild_cell = gdspy.Cell('unit_ILD_feed')
  ild_cell.add(ild)

  unitcell = gdspy.Cell('unit_feedline')
  unitcell.add('unit_GP_sub')
  unitcell.add(main_cell)
  unitcell.add(ild_cell)
  return unitcell

def getlinetopad(nrows):
  cells = main_lib.cell_dict
  unitfeed = getunitfeedline()
  dx, dy = get_size(unitfeed)
  n = roundto(edge_margin/dx, 10)
  hor_section = gdspy.CellArray(unitfeed, n, 1, [dx, dy])
  vfeed = gdspy.CellReference('MainFeedline', rotation=90)
  vdx, vdy = get_size(vfeed)
  vert_section = gdspy.CellArray(vfeed, 1, nrows - 1, [vdx, vdy])
  #feed_corner = gdspy.CellReference(cells['feed_corner'])
  #ild_corner = gdspy.CellReference(cells['ILD_corner'])
  #gndsub_corner = gdspy.CellReference(cells['gndsub_corner'])

  #moveabove(vert_section, gndsub_corner)
  #moveabove(vert_section, ild_corner)
  #moveabove(vert_section, feed_corner)

  #moveleft(ILD_corner, hor_section)

  terminus = gdspy.Cell('feedline_to_pad')
  terminus.add(hor_section)
  #terminus.add(vert_section)
  #terminus.add(gndsub_corner)
  #terminus.add(ILD_corner)
  #terminus.add(feed_corner)
  return terminus

def main():
  # Wafer organization all dimensions in microns
  nrows = 8
  ncols = 8
  N_res = nrows * ncols

  wafer = gdspy.Cell('Global_Overlay')
  ind = get_inductor()

  # Generate all the capacitors
  L = 10 #nH
  fstart = 300 # MHz
  df = 3 #MHz
  fs = fstart + df * np.arange(N_res)

  capfrac = 1/(fs/fs[0])**2
  Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
  idcs = [IDC(capfrac[i]) for i in range(N_res)]

  # Let's make the coupling capacitors as well
  Qc = Qi
  Ccs = 2 * np.sqrt((Cs * u.F)/((Z0/2) * Qc * 2*np.pi*fs*u.MHz)).to(u.F).value
  coupcaps = [IDC(1.0) for i in range(len(Cs))]

  caps = []
  coup_Nfingers = []
  Nfingers = []
  cfracs = []
  for idc, c, cs, coupcap, freq in zip(idcs, Cs, Ccs, coupcaps, fs):
    idc.finger_length = finger_length
    idc.cellname = "Capacitor_{0:3d}MHz".format(freq)
    nfinger, cfrac = getnumfingers(idc, c)
    Nfingers.append(nfinger)
    cfracs.append(cfrac)

    coupcap.finger_length = coup_cap_finger_length - 2
    nfingers, _ = getnumfingers(coupcap, cs)
    nfingers = int(round(nfingers/5)*5)
    coup_Nfingers.append(nfingers)


  coup_Nfingers = np.array(coup_Nfingers)
  unq_nfingers, unq_indices, unq_inv = np.unique(coup_Nfingers,\
      return_index=True, return_inverse=True)
  unq_coupcaps = []
  for index in unq_indices:
    coupcaps[index].nfinger = coup_Nfingers[index]
    coupcaps[index].cellname =\
        "coupling_capacitor_{0:d}".format(coup_Nfingers[index])
    unq_coupcaps.append(coupcaps[index].draw())

  Nfinger = int(round(np.max(Nfingers)/10)*10)
  common_Nfinger = np.min(Nfingers)
  common_capfrac = capfrac[-1]

  common_idc = IDC(1.0)
  common_idc.nfinger = common_Nfinger
  common_idc.cellname = "Capacitor_common"
  common_cap = common_idc.draw(less_one=True)
  obtained_Cs = []
  obtained_Ccs = []
  for ires, coup_index in zip(range(N_res), unq_inv):
    idcs[ires].nfinger = Nfinger - common_Nfinger
    idcs[ires].capfrac = (idcs[ires].capfrac - common_capfrac)/(1 - common_capfrac)
    caps.append(get_cap_tank(idcs[ires].draw(), unq_coupcaps[coup_index]))
    obtained_Cs.append(idcs[ires].C)
    obtained_Ccs.append(coupcaps[ires].C)

  # We will need the size of the capacitor when positioning it finally
  cap_ref = gdspy.CellReference(caps[0],rotation=90)
  u_dx, u_dy = get_size(cap_ref)

  # Make the common portion of the resonator
  common_resonator = get_common_resonator(ind, common_cap)

  cres_dx, cres_dy = get_size(common_resonator)
  cres_x_uneven = 1 #xmax - (-xmin) for the resonator

  # Make an array from the common resonator
  xspacing = roundto(wafer_width/ncols, 100)
  yspacing = roundto(wafer_len/nrows, 100)
  arr_xshift = (xspacing/2)*(ncols-1)
  arr_yshift = (yspacing/2)*(nrows-1)
  common_res_arr = gdspy.CellArray(common_resonator, ncols, nrows, [xspacing, yspacing] )
  common_res_arr.translate(-arr_xshift, -arr_yshift)
  # Make the GND plane subtract region around each resonator
  cres_bbox = common_resonator.get_bounding_box()

  # Generate the region around the resonator where the ground plane is to be
  # removed. Also return the offsets in the x direction from the center of this
  # region where the common resonator and the unique capacitor are to be placed.
  reso_gnd_sub, c_xoffset, u_xoffset, feed_xoffset =\
      get_resonator_GPsub(cres_bbox, get_size(cap_ref))
  g_dx, g_dy = get_size(reso_gnd_sub)
  g_dy -= feedres_dist

  reso_gnd_sub_arr = gdspy.CellArray(reso_gnd_sub, ncols, nrows, [xspacing, yspacing] )
  reso_gnd_sub_arr.translate(-arr_xshift, -arr_yshift)
  wafer.add(reso_gnd_sub_arr)
  cu_overlap = 2
  common_res_arr.translate(c_xoffset, 0)
  wafer.add(common_res_arr)

  # The position of the bottom left corner of the resonator array which is the
  # origin of the i, j coordinate system
  x_origin = -arr_xshift
  y_origin = -arr_yshift
  # Need to now map the variable capacitor pieces relative to the full array
  for ires in range(N_res):
    i, j = d2ij(ires, nrows)
    xpos = x_origin + i*xspacing - u_xoffset
    ypos = y_origin + j*yspacing
    cap_ref = gdspy.CellReference(caps[ires], [xpos, ypos], rotation=90)
    #cap2feed_ref = gdspy.CellReference(cap2feed)
    #moveabove(cap_ref, cap2feed_ref, spacing=2)
    #cap2gnd_ref = gdspy.CellReference(cap2gnd)
    #movebelow(cap_ref, cap2gnd_ref, spacing=-2)
    wafer.add(cap_ref)
    #wafer.add(cap2feed_ref)
    #wafer.add(cap2gnd_ref)



  #Make the feedline
  feed_cell = get_feedline()
  f_dx, f_dy = get_size(feed_cell)
  fncols = int((xspacing//f_dx)*ncols) - 1
  fnrows = nrows
  fxspacing = f_dx
  fyspacing = yspacing
  # Need a small correction to center on the whole resonator instead of common
  feed_array = gdspy.CellArray(feed_cell, fncols, fnrows, [fxspacing, fyspacing])
  feed_array.translate(-arr_xshift, -arr_yshift + g_dy/2 + feedres_dist )
  wafer.add(feed_array)
  (fxmin, fymin), (fxmax, fymax) = feed_array.get_bounding_box()

  vertfeedcell = make_vertical_feed(feed_cell)
  leftconns = gdspy.CellArray(vertfeedcell, 1, nrows//2, [0, 2*fyspacing])
  rightconns = gdspy.CellArray(vertfeedcell, 1, nrows//2-1,[0, 2*fyspacing])

  linetopad = getlinetopad(nrows)
  pad = get_bondpads()
  p_dx, p_dy = get_size(pad)
  top_pad_ref = gdspy.CellReference(pad, [fxmax + p_dx/2, fymax - f_dy/2])
  bot_pad_ref = gdspy.CellReference(pad, [fxmax + p_dx/2, fymin - f_dy/2])

  leftconns.translate(fxmin , fymin + yspacing/2 + f_dy/2)
  rightconns.translate(fxmax , fymin + 1.5* yspacing + f_dy/2)

  wafer.add(leftconns)
  wafer.add(rightconns)
  wafer.add(top_pad_ref)
  wafer.add(bot_pad_ref)

  cap2feed, cap2gnd = get_coupling_conns(g_dy, u_dy)
  cap2feed_dx, cap2feed_dy = get_size(cap2feed)
  cap2gnd_dx, cap2gnd_dy = get_size(cap2gnd)
  overlap = 2
  cap2feed_xoffset = u_xoffset - (u_dx//2 + u_dx % 2) + 8

  cap2feed_arr = gdspy.CellArray(cap2feed, ncols, nrows, [xspacing, yspacing])
  cap2feed_arr.translate(-arr_xshift - feed_xoffset, -arr_yshift +\
      u_dy/2+ cap2feed_dy/2 - overlap)
  wafer.add(cap2feed_arr)
  cap2gnd_arr = gdspy.CellArray(cap2gnd, ncols, nrows, [xspacing, yspacing])
  cap2gnd_arr.translate(-arr_xshift - feed_xoffset, -arr_yshift -\
      u_dy/2 - cap2gnd_dy/2 + overlap)
  wafer.add(cap2gnd_arr)
  chip_outline = makechipoutline(wafer_width, wafer_len)
  wafer.add(gdspy.CellReference(chip_outline, (-wafer_width/2, -wafer_len/2)))


  #make_left_connectors(vertfeedcell, )
  # Generate the patch and shot table

  all_cells = main_lib.cell_dict
  make_inverted_cells()

  inv_cell_list = get_inverted_cells()

  mask_len = 22000
  mask_width = 26000

  total_mask_area = mask_len * mask_width
  total_area_needed = np.sum(list(map(lambda x: get_cell_area(all_cells[x]),\
      inv_cell_list)))
  main_lib.write_gds('sscale_darkres.gds',unit=1e-6,precision=1e-9)








if __name__=='__main__':
  main()


