#!/usr/bin/env python3

import numpy as np
import gdspy
import sys
import matplotlib.pyplot as plt
import astropy.units as u
import patches
from smallscaledarkresonators import get_size, roundto
number_of_points = 32

mUnits = 1e-6
inv_margin = 200
coup_cap_finger_length = 100
bondpad_size = 140

gnd_box_margin = 200

def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12}

main_lib = gdspy.GdsLibrary('main')
mask_width = 22000
mask_length = 26000

def cellIsPresent(possible_cell_names):
    allcellnames = set(main_lib.cell_dict.keys())
    return any(name in allcellnames for name in possible_cell_names)

def inverter(cell, rotation=0):
    possible_cell_names = set([cell.name + '_r', cell.name + '_R', 'R_' + cell.name])
    if cellIsPresent(possible_cell_names): return
    cell_name = cell.name + '_r'
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
    main_lib.add(inv_cell)

def invert_cell(cell, rotation=0):
    layers = cell.get_layers()

    if len(layers) == 1:
        inverter(cell, rotation)

    for cell in cell.get_dependencies():
        invert_cell(cell, rotation)

def get_cell_list(element, cell_list):
    if not element.get_dependencies():
        cell = element.ref_cell
        main_lib.add(cell)
        cell_list.append(element)
    elements = element.elements
    for el in elements:

        cell_list.extend(get_cell_list(el))
    return cell_list

def get_mask(fn, mask_name):
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)

    return gdsii.extract(mask_name)

if __name__ == "__main__":
    fn = 'corwin_wafer_edited.gds'
    mask = get_mask(fn, 'reticleGND')
    cell_list = []
    to_invert = get_cell_list(mask, cell_list)
    for cell in to_invert:
        invert_cell(cell)

    main_lib.write_gds('corwin_wafer_inverted.gds',unit=1e-6,precision=1e-9)   
