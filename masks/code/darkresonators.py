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
finger_length = 900
indboxwidth = 200.
indboxheight = 280.
coup_cap_finger_length = 100
bondpad_size = 140
island_halfwidth = 400
Qi = 40000
Z0 = 50 * u.Ohm

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
    self.C = self.capacitance()

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
  specs = {'layer':10}
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

# Lambda functions
rounded_five = lambda x: int(round(x/5)*5)
rounded_ten = lambda x: int(round(x/10)*10)
rounded_even = lambda x: int(round(x/2)*2)
get_size = lambda x: (x[1,0] - x[0,0], x[1,1]-x[0,1])
get_coords = lambda x: (x[0,0], x[1,1])

def makerestankmodel(nfingers):
  model_cap = IDC(1.0)
  model_cap.nfinger = nfingers
  if model_cap.nfinger == 0:
    model_cap.nfinger = 1
  model_cap.finger_length = finger_length
  model_cap.cellname = 'Model_Res_Tank_%d' %model_cap.nfinger
  return model_cap

# I'm defining a series of template capacitors used to build up all the
# capacitors in the resonator array. The size of the largest capacitor is the
# minimum number of fingers of the capacitors in the array, rounded down to the
# nearest 10, if larger than 10. The smaller model capacitors are constructed by
# dividing the length of the largest capacitor by 2 until we construct a
# capacitor of length 1.
def make_captank_models(num, binpow):
  if num == 1:
    cap = makerestankmodel(num)
    return [cap.draw(less_one=True)], [num]

  if num == 0:
    return make_captank_models(1, binpow)

  cap = makerestankmodel(num)
  cell = cap.draw(less_one=True)

  caps, nfingers = make_captank_models(int(num // binpow), binpow )
  return [cell] + caps, [cap.nfinger] + nfingers

def update_origins(cap_array):
  curr = (0, 0)
  for i, arr in enumerate(cap_array):
    if i == 0:
      curr = arr.origin
      continue
    dx, dy = get_size(arr.get_bounding_box())
    curr = (curr[0] - dx +2, curr[1])
    arr.origin = curr

def make_capacitor(cap, model, model_fingers):
  cap_array = []
  origin = (0,0)
  # Loop over each model capacitor constructing the full capacitor
  res = cap.nfinger
  for m, nf in zip(model, model_fingers):
    N = res // nf
    if N == 0:
      continue
    res %= nf#Number of fingers remaining
    # Make the array as a row of the ref cells lined up
    dx, dy = get_size(m.get_bounding_box())
    array = gdspy.CellArray(m, 1, N, (dx , dy - cap.gap_width), rotation=-90)
    array.origin = origin
    #print (N, nf)
    cap_array.append(array)
    if res == 0:
      break
  update_origins(cap_array)
  return cap_array


def get_cap_position(cap, index, N_res, nrows, ncols, w_spacing, l_spacing,\
    ind_start, space):
  dx, dy = get_size(cap.get_bounding_box())
  row = index  // ncols
  col = index - (row * ncols)
  #index += 1
  #row = (N_res - index) // nrows
  #col = (N_res - index) - (row * ncols)
  # First let's locate the inductor associated with the capacitor
  x_ind = ind_start[0] + col* (w_spacing)
  y_ind = ind_start[1] - row * (l_spacing)
  x = x_ind - space - dx #space is the distance between the cap and ind
  y = y_ind + (dy - indboxheight)/2 - 2
  return (x,y)

def make_connector(len_conn):
  connector = gdspy.Path(2, (0,0))
  connector.segment(6, '+x', layer=1)
  connector.turn(2, 'r', layer=1)
  connector.segment(len_conn, '-y', layer=1)
  connector.turn(2, 'l', layer=1)
  connector.segment(island_halfwidth -1.5 * 2, '+x', layer=1)
  connector.turn(2, 'l', layer=1)
  connector.segment(11, '+y', layer=1)
  conn = gdspy.Cell('connectors')
  conn.add(connector)
  return conn

def add_connector_array(ncols, nrows, indspacing, start):
  ms_width = 2
  len_conn = (finger_length - indboxheight )/2 -ms_width + 61.5#magic number
  connector = make_connector(len_conn)
  (xmin, ymin), (xmax, ymax) = connector.get_bounding_box()
  dx = xmax - xmin - 0.25 * ms_width
  dy = ymax - ymin - ms_width
  origin = (start[0] - dx + ms_width, start[1] + dy +(indboxheight/2 + ms_width))
  return gdspy.CellArray(connector, ncols, nrows, indspacing, origin)

def add_row_feeds(ms_width, sf_start, sf_len, nrows, ncols, indspacing):
  # Now let's make the actual side feedline
  rfeed = gdspy.Path(ms_width)
  rfeed.segment(sf_len, '+x', layer=3)
  row_feed = gdspy.Cell('ms_feed_row')
  row_feed.add(rfeed)
  row_feed_arr = gdspy.CellArray(row_feed, columns=1,\
      rows=nrows,spacing=[sf_len, indspacing[1]])
  row_feed_arr.origin = sf_start

  return row_feed_arr

def get_left_connector_feeds(feed_conn, sf_start, num_lrows, spacing):
  left_conn_feeds = gdspy.CellArray(feed_conn, columns = 1, rows=num_lrows,\
      spacing= spacing)
  left_conn_feeds.origin = sf_start
  return left_conn_feeds

def get_right_connector_feeds(feed_conn, sf_start, conn_len, sf_len,\
    num_rrows, spacing):
  right_conn_feeds = gdspy.CellArray(feed_conn,columns = 1, rows=num_rrows ,\
      spacing= spacing)
  right_conn_feeds.origin = [sf_start[0] + sf_len, sf_start[1] + spacing[1]//2]
  return right_conn_feeds

def add_corner_connector_feeds(ms_width, sf_start, sf_len, nrows, ncols,\
    indspacing, conn_len, num_lrows, num_rrows):
  radius = ms_width * 1.5
  specturn = {'layer':3,'number_of_points':number_of_points}
  turn  = gdspy.Path(ms_width)
  turn.segment(ms_width, '-x', layer=3)
  turn.turn(radius, 'r', **specturn)
  turn.segment(ms_width, '+y', layer=3)
  turn_cell = gdspy.Cell('feed_corner')
  turn_cell.add(turn)
  turn_arr_bl = gdspy.CellArray(turn_cell, columns=1, rows = num_lrows,\
      spacing=[0, 2*indspacing[1]])
  turn_arr_bl.origin = [sf_start[0] + ms_width, sf_start[1] ]
  corner_spacing = conn_len + 3 * ms_width
  turn_arr_tl = gdspy.copy(turn_arr_bl, 0, (num_lrows * 2 - 1)*corner_spacing)
  turn_arr_tl.x_reflection = True

  # I'll make a second array out of the rotated turn_cell
  turn_cell_rot = gdspy.CellReference(turn_cell, rotation=180)
  turn_arr_tr = gdspy.CellArray(turn_cell, columns=1, rows = num_rrows,\
      spacing=[0, 2*indspacing[1]], rotation=180)
  dy = num_rrows * 2 *indspacing[1]
  turn_arr_tr.origin = [sf_start[0] + sf_len - 4*ms_width, sf_start[1] ]
  turn_arr_tr.translate(0, dy)
  turn_arr_br = gdspy.copy(turn_arr_tr, 0,  0)
  turn_arr_br.x_reflection = True
  turn_arr_br.translate(0, -dy + corner_spacing)

  corners = gdspy.Cell('all_corners')
  corners.add([turn_arr_bl, turn_arr_tl, turn_arr_br, turn_arr_tr])

  return gdspy.CellReference(corners)

def make_feedline_cells(sf_len, ms_width):
  bondpad = gdspy.Rectangle([0,0], [bondpad_size, -bondpad_size], layer=3)
  bondpad_cell = gdspy.Cell('ms_bondpad')
  bondpad_cell.add(bondpad)

  path = gdspy.Path(ms_width)
  path.segment(sf_len/2, '+x', layer=3)
  long_stretch = gdspy.Cell('ms_long_end')
  long_stretch.add(path)

  return [gdspy.CellReference(bondpad_cell), gdspy.CellReference(long_stretch)]

def make_feedline_end(feed_len, sf_len, ms_width, pos):
  final_feed = gdspy.Cell('final_feed' + pos)
  final_feed.add(gdspy.CellReference('ms_bondpad'))

  path = gdspy.Path(ms_width)
  path.segment(feed_len, '-y', layer=3)
  short_stretch = gdspy.Cell('ms_short_end' + pos)
  short_stretch.add(path)
  final_feed.add(gdspy.CellReference(short_stretch,\
      [bondpad_size/2, -bondpad_size]))

  curr = [bondpad_size/2, -bondpad_size - feed_len]
  final_feed.add(gdspy.CellReference('feed_corner', [curr[0] + 2.5*ms_width,\
    curr[1] - 1.5*ms_width]))
  curr = [curr[0] + 2.5*ms_width, curr[1] - 1.5*ms_width]

  final_feed.add(gdspy.CellReference('ms_long_end',\
      [curr[0] - ms_width, curr[1]]))
  curr = [curr[0] - ms_width + sf_len/2, curr[1]]

  final_feed.add(gdspy.CellReference('feed_corner', [curr[0] - ms_width,\
    curr[1]], rotation=180))
  curr = [curr[0] - ms_width, curr[1]]
  final_feed.add(gdspy.CellReference(short_stretch,\
      [curr[0] + 2.5 * ms_width, curr[1] - 1.5 * ms_width]))
  curr = [curr[0] + 2.5 * ms_width, curr[1] - 1.5 * ms_width]
  final_feed.add(gdspy.CellReference('feed_corner', [curr[0] - 2.5*ms_width,\
    curr[1] - feed_len - 1.5 * ms_width], rotation=180, x_reflection=True))
  return final_feed



def make_main_feedline(ms_width, sf_start, sf_len, nrows, ncols, indspacing):
  row_start = [sf_start[0] + 1.5 * ms_width, sf_start[1]]
  row_feed_arr = add_row_feeds(ms_width, row_start, sf_len -3*ms_width,\
      nrows, ncols, indspacing)
  conn_len = indspacing[1]- 3 * ms_width
  feed = gdspy.Path(ms_width)
  feed.segment(conn_len, '+y', layer=3)
  feed_conn = gdspy.Cell('short_feeds')
  feed_conn.add(feed)

  spacing = [indspacing[0], 2 * indspacing[1]]
  col_start = [sf_start[0], sf_start[1] + 1.5 * ms_width]
  num_lrows = nrows //2
  num_rrows = num_lrows
  if nrows % 2 == 0:
    num_rrows -= 1

  left_connectors = get_left_connector_feeds(feed_conn , col_start,\
      num_lrows, spacing)
  right_connectors = get_right_connector_feeds(feed_conn, col_start, conn_len,\
      sf_len, num_rrows, spacing)
  corners = add_corner_connector_feeds(ms_width, row_start, sf_len,\
      nrows, ncols, indspacing, conn_len, num_lrows, num_rrows)

  row_bbox = row_feed_arr.get_bounding_box()
  feeds_x_extent = row_bbox[1, 0] - row_bbox[0, 0]
  feeds_y_extent = row_bbox[1, 1] - row_bbox[0, 1] - ms_width

  feeds_y_margin = 300 # Distance from bondpad to edge of the chip
  available_space_top = wafer_len - row_bbox[1,1] - feeds_y_margin -0.5 * ms_width
  feed_len_top = (available_space_top - 3 * (1.5 * ms_width) - bondpad_size)/2
  final_feed_top = make_feedline_end(feed_len_top, sf_len, ms_width, '_top')
  # Calculate the position of the final stretch at the top
  ffeed_bbox = final_feed_top.get_bounding_box()
  dx_feed = ffeed_bbox[1, 0] - ffeed_bbox[0, 0] - 0.5 * ms_width
  dy_feed = ffeed_bbox[1, 1] - ffeed_bbox[0, 1] - 0.5 * ms_width
  top_feed = gdspy.CellReference(final_feed_top)

  available_space_bot = row_bbox[0,1] - feeds_y_margin  + 0.5 * ms_width
  feed_len_bot = (available_space_bot - 3 * (1.5 * ms_width) - bondpad_size)/2
  final_feed_bot = make_feedline_end(feed_len_bot, sf_len, ms_width, '_bot')
  # Calculate the position of the final stretch at the bot
  ffeed_bbox_bot = final_feed_bot.get_bounding_box()
  dx_feed_bot = ffeed_bbox_bot[1, 0] - ffeed_bbox_bot[0, 0] - 0.5 * ms_width
  dy_feed_bot = ffeed_bbox_bot[1, 1] - ffeed_bbox_bot[0, 1] - 0.5 * ms_width
  bot_feed = gdspy.CellReference(final_feed_bot)

  # If nrows is even, we connect the last piece of the feedline at the top right
  # end of the rows. If odd, it is connected at the top left end and the last
  # piece needs to be rotated into place.
  if nrows % 2 == 0:
    x_end = row_bbox[1, 0]
    y_end = row_bbox[1, 1] - 0.5 * ms_width
    origin = [x_end - dx_feed + 1.5* ms_width, y_end + dy_feed]
    bot_origin = [origin[0], row_bbox[0, 1] - dy_feed_bot + 0.5 * ms_width ]
  else:
    x_end = row_bbox[0, 0]
    y_end = row_bbox[1, 1] - 0.5 * ms_width
    origin = [x_end + dx_feed - 1.5 * ms_width, y_end + dy_feed]
    overlap = 2 * dx_feed - feeds_x_extent - bondpad_size + 0.5 * ms_width
    bot_origin = [origin[0] - overlap, row_bbox[0, 1] - dy_feed_bot + 0.5 * ms_width ]
    top_feed.x_reflection = True
    top_feed.rotation = 180
  top_feed.origin = origin
  bot_feed.x_reflection = True
  bot_feed.origin = bot_origin
  feedline = gdspy.Cell('mainfeedline')
  feedline.add(top_feed)
  feedline.add(bot_feed)
  feedline.add([row_feed_arr, left_connectors, right_connectors, corners])
  return gdspy.CellReference(feedline)

def make_ground_plane(ms_width, sf_start, sf_len, ms_start, nrows, ncols,\
    yspacing, cap_size, ind_size, y_edge_margin):
  # Making the ground trunk lines running on the sides of the chip
  gnd_width = 500
  gnd_feed = gdspy.Path(gnd_width)
  gnd_feed.segment(wafer_len - 2000, '-y', layer=3)
  gnd_mainfeed = gdspy.Cell('gnd_trunk_line')
  gnd_mainfeed.add(gnd_feed)
  gnd_start_left = (sf_start[0] // 2, ms_start[1])
  gnd_start_right = (wafer_width - gnd_start_left[0], ms_start[1])

  # Make the row tines that serve each row of capacitors
  gnd_sfeed = gdspy.Path(gnd_width)
  gnd_sfeed.segment(sf_len, '+x', layer=3)
  gnd_sidefeed = gdspy.Cell('gnd_row_line')
  gnd_sidefeed.add(gnd_sfeed)
  num_rrows = nrows //2
  num_lrows = nrows - num_rrows
  gnd_lsidefeedarr = gdspy.CellArray(gnd_sidefeed, columns=1,\
      rows=num_lrows ,spacing=[0, 2*yspacing])
  gnd_rsidefeedarr = gdspy.CellArray(gnd_sidefeed, columns=1,\
      rows=num_rrows , spacing=[0, 2*yspacing])
  gnd_sf_offset_y = (cap_size - ind_size)/2 + coup_cap_finger_length + gnd_width//2
  gnd_sf_start = (gnd_start_left[0], y_edge_margin - gnd_sf_offset_y)
  gnd_lsidefeedarr.origin = gnd_sf_start
  gnd_rsidefeedarr.origin = gnd_sf_start
  gnd_rsidefeedarr.translate(2 * gnd_start_left[0], yspacing)

  # Assemble the full ground plane
  gnd_plane = gdspy.Cell('Ground_Plane')
  gnd_plane.add(gnd_lsidefeedarr)
  gnd_plane.add(gnd_rsidefeedarr)
  gnd_plane.add(gdspy.CellReference(gnd_mainfeed, gnd_start_left))
  gnd_plane.add(gdspy.CellReference(gnd_mainfeed, gnd_start_right))

  return gdspy.CellReference(gnd_plane)

def get_inductor():
  fn = '../resobolo_files/new_inductor.gds'
  gdsii = gdspy.GdsLibrary()
  gdsii.read_gds(fn)
  ind = gdsii.extract('Inductor_Winding')
  polys = ind.get_polygons()
  (xmin, ymin), (xmax, ymax) = ind.get_bounding_box()
  #ind_view = gdspy.CellReference(ind, [-xmin, -ymin])
  inductor = gdspy.Cell('inductor')
  for poly in polys:
    polygon = gdspy.Polygon(poly, layer=0)
    polygon = polygon.translate(-xmin, -ymin)
    inductor.add(polygon)
  inductor.flatten()
  return inductor

def get_coupcap(cap_array, model_coupcap):
  coupcap_tofeed = gdspy.CellReference(model_coupcap, rotation = -90)
  coupcap_tognd = gdspy.CellReference(model_coupcap, rotation = -90)
  (xmin, ymin), (xmax, ymax) = coupcap_tofeed.get_bounding_box()
  (xmin_e, ymin_e), (xmax_e, ymax_e) = cap_array[0].get_bounding_box()
  (xmin_s, ymin_s), (xmax_s, ymax_s) = cap_array[-1].get_bounding_box()
  dx = xmax - xmin
  midpoint = (xmin_s + xmax_e)/2
  tofeed_origin = [midpoint- dx/2 , ymax_e + coup_cap_finger_length]
  tognd_origin = [midpoint- dx/2 , ymin_e]
  coupcap_tofeed.origin = tofeed_origin
  coupcap_tognd.origin = tognd_origin

  return [coupcap_tofeed, coupcap_tognd]

def get_array_size(arr):
  return arr.columns * arr.rows

def generate_caps_from_freqs(fstart, df, N_res, L):
  fs = (fstart + df * np.arange(N_res))
  fracs = 1/(fs/fs[0])**2
  Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
  caps = [IDC(1.0) for i in range(len(Cs))]

  # Let's make the coupling capacitors as well
  Qc = Qi
  Ccs = 2 * np.sqrt((Cs * u.F)/((Z0/2) * Qc * 2*np.pi*fs*u.MHz)).to(u.F).value
  coupcaps = [IDC(1.0) for i in range(len(Cs))]

  N_fingers = [0]*N_res
  coup_N_fingers = [0]*N_res
  for i in range(N_res):
    caps[i].finger_length = finger_length
    nfingers, capfrac = getnumfingers(caps[i], Cs[i])
    caps[i].nfinger = nfingers
    caps[i].capfrac = capfrac
    caps[i].C = caps[i].capacitance()
    N_fingers[i] = nfingers

    coupcaps[i].finger_length = coup_cap_finger_length - 2
    nfingers, capfrac = getnumfingers(coupcaps[i], Ccs[i])
    coup_N_fingers[i] = nfingers

  return fs, caps, Cs, N_fingers, coup_N_fingers



def main():
  # Wafer organization all dimensions in microns
  nrows = 18
  ncols = 18
  N_res = nrows * ncols
  width_spacing = np.around((wafer_width - \
      2*edge_margin )/ncols, 0)
  len_spacing = np.around((wafer_len - 2*edge_margin)/nrows, 0)
  chip_outline = makechipoutline(wafer_width, wafer_len)
  wafer = gdspy.Cell('6_inch_wafer')

  # Make the array of inductors first since it is simple enough
  #ind = InductorMeander()
  #indcell = ind.draw()
  indcell = get_inductor()
  (xmin, ymin), (xmax, ymax) = indcell.get_bounding_box()
  global indboxheight
  indboxheight = ymax - ymin
  global indboxwidth
  indboxwidth = xmax - xmin
  indspacing = [width_spacing, len_spacing]
  indarray = gdspy.CellArray(indcell,ncols, nrows, indspacing)
  (xmin, ymin), (xmax, ymax) = indarray.get_bounding_box()
  ind_start = (xmin, ymax)
  ind_dx, ind_dy = (xmax - xmin, ymax - ymin)
  indarray.translate((wafer_width - ind_dx)//2 + island_halfwidth,\
      (wafer_len - ind_dy)//2)
  indbbox = indarray.get_bounding_box()
  ind_start = (indbbox[0,0], indbbox[1,1])
  ind_dx, ind_dy = get_size(indbbox)
  x_edge_margin = indbbox[0, 0]
  y_edge_margin = indbbox[0, 1]
  wafer.add(gdspy.CellReference(chip_outline, (0,0)))

  # Add the connectors from the inductors to the capacitors
  wafer.add(indarray)
  conn_start = indbbox[0]
  top_connectors = add_connector_array(ncols, nrows, indspacing, conn_start)
  connector = top_connectors.ref_cell
  (conn_dx, conn_dy) = get_size(connector.get_bounding_box())
  bot_connectors = gdspy.copy(top_connectors, 0, -indboxheight)
  bot_connectors.x_reflection=True
  #bot_connectors.rotation = 180
  bot_connectors.translate(0, ind_dy - finger_length -4 )
  wafer.add(top_connectors)
  wafer.add(bot_connectors)

  # Generate all the capacitors
  L = 10 #nH
  fstart = 300 # MHz
  df = 2 #MHz
  fs, caps, Cs, N_fingers, coup_N_fingers = generate_caps_from_freqs(fstart, df, N_res, L)
  #unq, unq_counts = np.unique(N_fingers, return_counts=True)
  #fig, ax = plt.subplots(figsize=(10, 10))
  #ax.hist(N_fingers, bins=unq, histtype='step')
  #plt.show()

  #cap_fingers_dict = dict(zip(unq, unq_counts))
  #print (cap_fingers_dict)
  # Template capacitors for constructing all the capacitors
  bin_pow = 1.5
  power_max = int(np.log(np.max(N_fingers))/np.log(bin_pow))
  num_fingers_max =  int(bin_pow ** power_max)
  model_cap_cells, model_cap_nfingers = make_captank_models(num_fingers_max,
      bin_pow)
  print ("We have model capacitors of sizes ", model_cap_nfingers)
  #model_cap = makerestankmodel(np.min(N_fingers))
  #model_cap_cell = model_cap.draw(less_one=True)

  coup_nfinger_step = 5
  model_coupcaps = []
  coup_power_max = int(np.log(np.max(coup_N_fingers))/np.log(2))
  coup_num_fingers_max =  2 ** coup_power_max
  #model_coupcap_nfingers = 2**np.arange(coup_power_max + 1)
  model_coupcap_nfingers = np.arange(coup_nfinger_step,\
      np.max(coup_N_fingers) + 1, coup_nfinger_step)
  all_modelcoupcaps = []
  for nfingers in model_coupcap_nfingers:
    coupcap = IDC(1.0)
    coupcap.nfinger = nfingers
    coupcap.finger_length = coup_cap_finger_length - 2
    coupcap.C = coupcap.capacitance()
    coupcap.cellname = 'Model_CoupRes_Tank_%d' %nfingers
    all_modelcoupcaps.append(coupcap)
    model_coupcaps += [coupcap.draw(less_one=True)]

  print ("We have model coupling capacitors of sizes ", model_coupcap_nfingers)
  all_capacitors = []
  indices = []
  for i in (range(N_res)):
    cap_array = make_capacitor(caps[i], model_cap_cells, model_cap_nfingers)
    origin= get_cap_position(cap_array[0], N_res - (i+1), N_res, nrows, ncols,\
        width_spacing, len_spacing, ind_start, conn_dx - 1.5*2)
    cap_array[0].origin = origin
    update_origins(cap_array)
    # Make and place the coupling capacitor
    coup_nfingers = coup_N_fingers[i]
    #Find which coupcap to use for this cap array
    #index = int(np.log2(coup_nfingers))
    index = rounded_five(coup_nfingers) //coup_nfinger_step  - 1
    indices.append(index)
    coup_caps = get_coupcap(cap_array, model_coupcaps[index])
    all_capacitors.append(cap_array)
    wafer.add(coup_caps)
    wafer.add(cap_array)

  unq, unq_counts = np.unique(indices, return_counts=True)
  #print (unq, unq_counts)
  num_caps = list(map(lambda x: np.sum(list(map\
      (lambda y: get_array_size(y),x))), all_capacitors))
  #print (num_caps)
  print ("Avg number of capacitors shot per full capacitor", np.mean(num_caps))
  print ("Peak number of shots to make a full capacitor", np.max(num_caps))
  print ("The total number of capacitors to be shot is ", np.sum(num_caps))

  # I want to obtain the coupling efficiency for each capacitor based on the
  # current coupling scheme.
  Qcs = []
  coup_fingers = []
  for i in range(N_res):
    C = caps[i].C
    Cc = all_modelcoupcaps[indices[i]].C
    coup_fingers.append(all_modelcoupcaps[indices[i]].nfinger)
    Qc = ((C*u.F)/((Z0/2)*2*np.pi*fs[i]*u.MHz*(Cc/2 * u.F)**2)).to(1).value
    Qcs.append(Qc)

  Qcs = np.array(Qcs)
  rho_c = Qcs/Qi
  chi_c = 4 * rho_c/(1 + rho_c)**2
  #print (Qcs)

  #fig, ax = plt.subplots(figsize=(10, 10))
  #ax.plot(chi_c, 'bd' )
  #ax.plot(coup_fingers, 'b', label='realized')
  #ax.plot(coup_N_fingers, 'r', label='expected')
  #ax.legend(loc='best')
  #plt.show()
  #wafer.add(model_cap_ref)


  #for i in (range(N_res)):
  #  coupcap_array = make_capacitor(coupcaps[i], model_cap_cells, model_cap_nfingers)
  #  origin= get_cap_position(cap_array[0], N_res - (i+1), N_res, nrows, ncols,\
  #      width_spacing, len_spacing, ind_start)
  #  cap_array[0].origin = origin
  #  update_origins(cap_array)
  #  all_cap_arrays.append(cap_array)
  #  wafer.add(cap_array)


  # Feedlines
  ms_start = (wafer_width //2 , wafer_len - 1000)
  ms_width = 40

  # First let's calculate the position of the side tines connecting to the main
  # microstrip
  cap_size = finger_length + 6
  ind_size = indboxheight
  max_nfingers = np.max(N_fingers)
  sf_offset_y = (cap_size - ind_size)/2 + ind_size + 100 + ms_width //2
  sf_offset_x = max_nfingers * 4 + 6
  sf_start = (ind_start[0] - sf_offset_x, y_edge_margin + sf_offset_y)
  sf_len = wafer_width - 2*sf_start[0]
  make_feedline_cells(sf_len, ms_width)
  main_feeds = make_main_feedline(ms_width, sf_start, sf_len, nrows, ncols,\
      indspacing)
  wafer.add(main_feeds)

  gnd_plane = make_ground_plane(ms_width, sf_start, sf_len, ms_start,\
      nrows, ncols, indspacing[1], cap_size, ind_size, y_edge_margin)
  wafer.add(gnd_plane)


  gdspy.write_gds('darkres.gds',unit=1e-6,precision=1e-9)

  # Let's make a map of the number of unit caps needed along the array
  #grid = np.arange(N_res)[::-1].reshape((nrows, ncols))
  #z = np.array(N_fingers)[::-1].reshape((nrows, ncols)) // model_cap.nfinger
  #img = plt.imshow(z, cmap='jet', aspect='equal', interpolation=None)
  ##img.set_cmap('jet')
  #plt.colorbar()
  #plt.show()

if __name__=='__main__':
  main()


