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
finger_length = 600
indboxwidth = 200.
indboxheight = 280.
coup_cap_finger_length = 100
bondpad_size = 140


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

# Lambda functions
rounded_ten = lambda x: int(x/10)*10
rounded_even = lambda x: int(x/2)*2
get_size = lambda x: (x[1,0] - x[0,0], x[1,1]-x[0,1])
get_coords = lambda x: (x[0,0], x[1,1])

def makerestankmodel(nfingers):
  model_cap = IDC(1.0)
  model_cap.nfinger = rounded_ten(nfingers) if nfingers >= 10 else \
      rounded_even(nfingers)
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
def make_captank_models(num):
  if num == 1:
    cap = makerestankmodel(num)
    return [cap.draw(less_one=True)], [num]

  if num == 0:
    return make_captank_models(1)

  cap = makerestankmodel(num)
  cell = cap.draw(less_one=True)

  caps, nfingers = make_captank_models(num // 2)
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
    ind_start):
  dx, dy = get_size(cap.get_bounding_box())
  row = index  // ncols
  col = index - (row * ncols)
  #index += 1
  #row = (N_res - index) // nrows
  #col = (N_res - index) - (row * ncols)
  space = 6 #small gap between cap and ind
  # First let's locate the inductor associated with the capacitor
  x_ind = ind_start[0] + col* (w_spacing)
  y_ind = ind_start[1] - row * (l_spacing)
  x = x_ind - space - dx
  y = y_ind + (dy - indboxheight  + 5)/2 - 2

  return (x,y)

def make_connector(len_conn):
  connector = gdspy.Path(2, (0,0))
  connector.segment(6, '+x', layer=1)
  connector.turn(2, 'r', layer=1)
  connector.segment(len_conn, '-y', layer=1)
  conn = gdspy.Cell('connectors')
  conn.add(connector)
  return conn

def add_connector_array(ncols, nrows, indspacing, start):
  len_conn = (finger_length + 2 - indboxheight + 5 )/2 -1
  connector = make_connector(len_conn)
  origin = (start[0] - 7, start[1] + indboxheight -5 + len_conn + 2)
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
  print (available_space_top)
  feed_len_top = (available_space_top - 3 * (1.5 * ms_width) - bondpad_size)/2
  final_feed_top = make_feedline_end(feed_len_top, sf_len, ms_width, '_top')
  # Calculate the position of the final stretch at the top
  ffeed_bbox = final_feed_top.get_bounding_box()
  dx_feed = ffeed_bbox[1, 0] - ffeed_bbox[0, 0] - 0.5 * ms_width
  dy_feed = ffeed_bbox[1, 1] - ffeed_bbox[0, 1] - 0.5 * ms_width
  top_feed = gdspy.CellReference(final_feed_top)

  available_space_bot = row_bbox[0,1] - feeds_y_margin  + 0.5 * ms_width
  print (available_space_bot)
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
  (xmin, ymin), (xmax, ymax) = ind.get_bounding_box()
  ind_view = gdspy.CellReference(ind, -xmin, -ymin)
  inductor = gdspy.Cell('inductor')
  inductor.add(ind_view)
  inductor.flatten()
  return inductor

def main():
  # Wafer organization all dimensions in microns
  nrows = 18
  ncols = 18
  N_res = nrows * ncols
  width_spacing = np.around((wafer_width - 2*edge_margin)/ncols, 0)
  len_spacing = np.around((wafer_len - 2*edge_margin)/nrows, 0)
  chip_outline = makechipoutline(wafer_width, wafer_len)
  wafer = gdspy.Cell('6_inch_wafer')

  # Make the array of inductors first since it is simple enough
  # ind = InductorMeander()
  # indcell = ind.draw()
  indcell = get_inductor()
  indspacing = [width_spacing, len_spacing]
  indarray = gdspy.CellArray(indcell,ncols, nrows, indspacing)
  (xmin, ymin), (xmax, ymax) = indarray.get_bounding_box()
  ind_start = (xmin, ymax)
  ind_dx, ind_dy = (xmax - xmin, ymax - ymin) 
  indarray.translate((wafer_width - ind_dx)//2, (wafer_len - ind_dy)//2)
  indbbox = indarray.get_bounding_box()
  ind_start = (indbbox[0,0], indbbox[1,1])
  ind_dx, ind_dy = get_size(indbbox)
  x_edge_margin = indbbox[0, 0]
  y_edge_margin = indbbox[0, 1]
  wafer.add(gdspy.CellReference(chip_outline, (0,0)))


  # Let's think about the frequency spacing
  df = 2 #MHz
  fstart = 300 #MHz
  fs = (fstart + df * np.arange(N_res))
  fracs = 1/(fs/fs[0])**2
  L = 10 #nH
  Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
  caps = [IDC(1.0) for i in range(len(Cs))]

  N_fingers = [0]*N_res
  for i in range(N_res):
    caps[i].finger_length = finger_length
    nfingers, capfrac = getnumfingers(caps[i], Cs[i])
    caps[i].nfinger = nfingers
    caps[i].capfrac = capfrac
    N_fingers[i] = nfingers

  # Template capacitors for constructing all the capacitors
  model_cap_cells, model_cap_nfingers = make_captank_models(np.min(N_fingers))
  print ("We have model capacitors of sizes ", model_cap_nfingers)
  #model_cap = makerestankmodel(np.min(N_fingers))
  #model_cap_cell = model_cap.draw(less_one=True)

  all_cap_arrays = []
  for i in (range(N_res)):
    cap_array = make_capacitor(caps[i], model_cap_cells, model_cap_nfingers)
    origin= get_cap_position(cap_array[0], N_res - (i+1), N_res, nrows, ncols,\
        width_spacing, len_spacing, ind_start)
    cap_array[0].origin = origin
    update_origins(cap_array)
    all_cap_arrays.append(cap_array)
    wafer.add(cap_array)

  all_cap_arrays = np.array(all_cap_arrays)
  #wafer.add(model_cap_ref)
  wafer.add(indarray)
  conn_start = indbbox[0]
  top_connectors = add_connector_array(ncols, nrows, indspacing, conn_start)

  bot_connectors = gdspy.copy(top_connectors, 0, -indboxheight)
  bot_connectors.x_reflection=True
  #bot_connectors.rotation = 180
  bot_connectors.translate(0, ind_dy - finger_length +1 )
  wafer.add(top_connectors)
  wafer.add(bot_connectors)

  # Let's make the coupling capacitors
  Qi = 40000
  Qc = Qi
  Z0 = 50 * u.Ohm
  Ccs = 2 * np.sqrt((Cs * u.F)/((Z0/2) * Qc * 2*np.pi*fs*u.MHz)).to(u.F).value
  coupcaps = [IDC(1.0) for i in range(len(Cs))]

  coup_N_fingers = [0]*N_res
  for i in range(N_res):
    coupcaps[i].finger_length = coup_cap_finger_length
    nfingers, capfrac = getnumfingers(coupcaps[i], Ccs[i])
    coupcaps[i].nfinger = nfingers
    coupcaps[i].capfrac = capfrac
    coup_N_fingers[i] = nfingers

  #print (coup_N_fingers)
  #all_coupcap_arrays = []
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
  ind_size = indboxheight - 5
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


