import gdspy
import numpy as np
import pdb
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12}
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
            return "Shot({0}, {1}, {2}, {3})".format(self.cellname,\
                    self.cell_shift,self.mask_shift, self.maskcellsize)
            return "Shot({0}, {1}, {2}, {3}, {4}, {5})".format(self.cellname,\
                    self.center, self.mask_shift, self.maskcellsize,\
                    (self.ncols, self.nrows),(self.xspacing, self.yspacing))

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
    shot.cell_shift = (dx - sdx, dy - sdy)

def get_size(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    dx = (xmax - xmin)
    dy = (ymax - ymin)
    return dx, dy

def get_center(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    x0 = (xmax + xmin)//2
    y0 = (ymax + ymin)//2
    return x0, y0

def makeshot(element, parent=None, hierarchy=0):
    isArray = False
    if type(element) == gdspy.CellArray:
        isArray = True
    elif type(element) != gdspy.CellReference:
        return []

    cell = element.ref_cell
    cellname = cell.name
    # If a cell has dependencies, then the elements gives a list of all the cell
    # references in this cell. If the cell has no dependencies, then subelements
    # just gives the polygon set that makes up the cell. I don't want the
    # polygonset stuff.
    subelements = []
    if cell.get_dependencies():
        subelements = cell.elements

    if isArray:
        arr_center = scalearr(get_center(element), scale)
        ncols = element.columns
        nrows = element.rows
        xspacing, yspacing = scalearr(element.spacing, scale)
        args = {'num_cols':ncols, 'num_rows':nrows,\
            'center':arr_center, 'xspacing':xspacing, 'yspacing':yspacing}

    # If deps is an empty set then immediately construct the shot from the
    # element. Otherwise loop through the dependencies of the current element and
    # obtain shots from each of them.
    if not subelements:
        cell_shift = scalearr(element.origin, scale)
        cell_size = scalearr(get_size(element), scale)
    else:
        shotlist = []
        for el in subelements:
            shots = makeshot(el, parent=element, hierarchy=hierarchy + 1)
            # For each shot that has been made from each element, update its origin
            # relative to the origin for  the element
            if hierarchy >= 2:
                translate(shot, scalearr(parent.origin, scale))
            shotlist.extend(shots)
        # If the current shot is part of a larger array, I want to keep track of
        # that information
        if isArray:
            for shot in shotlist:
                shot.isArray = True
                shot.center = arr_center
                shot.ncols = ncols
                shot.nrows = nrows
                shot.xspacing, shot.yspacing = xspacing, yspacing
        return shotlist

    if isArray:
        args = {'num_cols':ncols, 'num_rows':nrows,\
            'center':arr_center, 'xspacing':xspacing, 'yspacing':yspacing}
        shot = Shot(cell, cell_shift, cell_size, isArray=True, **args)
    else:
        args = {'num_cols':1, 'num_rows':1,\
            'center':cell_shift, 'xspacing':0, 'yspacing':0}
        shot = Shot(cell, (0, 0), cell_size, isArray=True, **args)

    if hierarchy >= 2:
        translate(shot, scalearr(parent.origin, scale))
    return [shot]















