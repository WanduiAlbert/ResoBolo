import gdspy
import numpy as np

class Shot():

  def __init__(self, cell, cell_shift, is_arrayed=False, **kwargs):
    try:
      assert(type(cell) == gdspy.Cell)
    except TypeError:
      print ("Please provide a valid cell to construct the shot")
    self.cell = cell
    self.cellname = cell.cellname
    layers = list(cell.get_layers())
    if len(layers) > 1 :
      raise RuntimeError ("Cannot construct a shot from a cell with multiple layers")
    self.layer = layers[0]
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    dx = (xmax - xmin)
    dy = (ymax - ymin)
    self.cell_size = (dx, dy)
    self.cell_shift = cell_shift
    self.is_arrayed = is_arrayed
    if is_arrayed:
      self.ncols = kwargs['num_cols']
      self.nrows = kwargs['num_rows']
      self.center = kwargs['center']
      self.stepsize = kwargs['stepsize']
