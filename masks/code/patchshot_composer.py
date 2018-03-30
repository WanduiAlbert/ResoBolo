import gdspy
import numpy as np
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
    "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9}
inv_layers = {value:key for key, value in def_layers.items()}

# scales all the values in an array by the scale
def scalearr(arr, scale):
  return tuple(map(lambda x: x/scale, arr))

scale = 1000

class Shot():

  def __init__(self, cell, cell_shift, cell_size, isArray=False, **kwargs):
    try:
      assert(type(cell) == gdspy.Cell or type(cell) == gdspy.CellReference)
    except AssertionError:
      print ("Please provide a valid cell to construct the shot")
    self.cell = cell
    self.cellname = cell.name
    layers = list(cell.get_layers())
    if len(layers) > 1 :
      raise RuntimeError ("Cannot construct a shot from a cell with multiple layers")
    self.layer = layers[0]
    self.cell_size = cell_size
    self.cell_shift = cell_shift
    self.isArray = isArray
    #self.blade_coords = {'xl':, 'xr':, 'yt':, 'yb':}
    if isArray:
      self.ncols = kwargs['num_cols']
      self.nrows = kwargs['num_rows']
      self.center = kwargs['center']
      self.xspacing = kwargs['xspacing']
      self.yspacing = kwargs['yspacing']

  def update_mask_location(self, maskcellref):
    maskcellname = maskcellref.ref_cell.name
    assert(maskcellname.startswith(self.cellname))

    self.maskcell = maskcellref.ref_cell
    self.maskcellsize = scalearr(get_size(maskcellref), scale)
    self.mask_shift = scalearr(maskcellref.origin, scale)

  def __repr__(self):
    if not self.isArray:
      return "Shot({0}, {1}, {2}, {3})".format(self.cellname,self.cell_shift,\
          self.mask_shift, self.maskcellsize)
    return "Shot({0}, {1}, {2}, {3}, {4}, {5})".format(self.cellname,\
         self.center, self.mask_shift, self.maskcellsize, (self.ncols, self.nrows),\
        (self.xspacing, self.yspacing))

  def __lt__(self, other):
    if self.layer == other.layer:
      return self.cellname < other.cellname
    return self.layer < other.layer

  def get_layer(self):
    return inv_layers[self.layer]

  def get_name(self):
    return self.cellname

  def get_mask_shift(self):
    return self.mask_shift

  def get_cell_size(self):
    return self.maskcellsize

  def get_cell_shift(self):
    return self.cell_shift

  def get_array_size(self):
    if not self.isArray: return None, None
    return self.ncols, self.nrows

  def get_array_center(self):
    if not self.isArray: return None, None
    return self.center

  def get_array_stepsize(self):
    if not self.isArray: return None, None
    return self.xspacing, self.yspacing

def translate(shot, shift):
  dx, dy = shot.cell_shift
  sdx, sdy = shift
  shot.cell_shift = (dx + sdx, dy + sdy)

def get_size(cell):
  (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
  dx = (xmax - xmin)
  dy = (ymax - ymin)
  return dx, dy

def makeshot(element):
  isArray = False
  if type(element) == gdspy.CellArray:
    isArray = True
  elif type(element) != gdspy.CellReference:
    return []

  cell = element.ref_cell
  cellname = cell.name
  deps = cell.get_dependencies()

  # If deps is an empty set then immediately construct the shot from the
  # element. Otherwise loop through the elements of the current element and
  # obtain shots from each of them.
  if not deps:
    cell_shift = scalearr(element.origin, scale)
    cell_size = scalearr(get_size(element), scale)
  else:
    shotlist = []
    for el in cell.elements:
      shots = makeshot(el)
      # For each shot that has been made from each element, update its origin
      # relative to the origin for  the element
      for shot in shots:
        translate(shot, scalearr(el.origin, scale))
      shotlist.extend(shots)
    # Update again because all the element are positioned relative to the whole
    # global overlay
    for shot in shotlist:
      translate(shot, scalearr(element.origin, scale))
      if isArray:
        shot.isArray = True
        shot.center = shot.cell_shift
        shot.ncols = element.columns
        shot.nrows = element.rows
        shot.xspacing, shot.yspacing = scalearr(element.spacing, scale)
    return shotlist

  if isArray:
    cell_size = get_size(cell)
    cell_shift = scalearr(element.origin, scale)
    arr_origin = scalearr(element.origin, scale)
    ncols = element.columns
    nrows = element.rows
    xspacing, yspacing = scalearr(element.spacing, scale)
    args = {'num_cols':ncols, 'num_rows':nrows,\
        'center':arr_origin, 'xspacing':xspacing, 'yspacing':yspacing}
    shot = Shot(cell, cell_shift, cell_size, isArray=True, **args)
  else:
    shot = Shot(cell, cell_shift, cell_size, isArray=isArray)

  return [shot]















