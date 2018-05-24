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
    possible_cell_names = set([cell.name + '_r'])
    if cellIsPresent(possible_cell_names): return
    cell_name = cell.name + '_r'
    inv_cell = gdspy.Cell(cell_name)
    cell = cell.flatten()
    cell_ref = gdspy.CellReference(cell, rotation=rotation)
    (xmin, ymin), (xmax, ymax) = cell_ref.get_bounding_box()
    xmin = roundto(xmin - inv_margin, 50)
    ymin = roundto(ymin - inv_margin, 50)
    xmax = roundto(xmax + inv_margin, 50)
    ymax = roundto(ymax + inv_margin, 50)
    layer = cell.get_layers().pop()
    polys = cell_ref.get_polygons(depth=1)
    polyset = gdspy.PolygonSet(polys, layer=layer)
    bbox = gdspy.Rectangle([xmin, ymin], [xmax, ymax], layer=layer)
    new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
    inv_cell.add(new_polyset)
    main_lib.add(inv_cell)

def invert_cell(cell, rotation=0):
    deps = cell.get_dependencies()    
    if not deps:
        inverter(cell, rotation)
    else:
        for dep in deps:
            invert_cell(deps, rotation)

def get_node_list(root_element, parent_origin=[0,0]):
    root_cell = root_element.ref_cell
    root_origin = root_element.origin
    # Find the absolute origin
    # abs_origin = root_origin
    abs_origin = patches.cellNode.get_new_origin(parent_origin,\
    root_origin)

    haschildren = bool(root_cell.get_dependencies())
    if not haschildren:
        node = patches.cellNode(root_cell, abs_origin)
        node.update_origin(parent_origin)
        return [node]

    children = root_cell.elements

    child_nodes = []
    for child in children:
        child_nodes += get_node_list(child, abs_origin)
    return child_nodes

def get_cell_list(cell, cell_list):
    deps = cell.get_dependencies()
    if not deps:
        main_lib.add(cell)
        cell_list.append(cell)
        return
    
    for dep in deps:
        get_cell_list(dep, cell_list)
    return

def get_mask(fn, mask_name):
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)

    return gdsii.extract(mask_name)

def generate_inverted_mask(node_list):
    inv_list = []
    for node in node_list:
        node_cellname = node.cell.name
        inv_node_cell = main_lib.cell_dict[node_cellname + '_r']
        inv_cell = gdspy.CellReference(inv_node_cell, origin=node.origin)
        inv_list.append(inv_cell)
    return inv_list

if __name__ == "__main__":
    fn = 'corwin_wafer_edited.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    mask = gdsii.extract('reticleGND')
    global_overlay = gdsii.extract('wafer')
    mask2 = gdsii.extract('reticleMS')
    pixel = gdsii.extract('pixel')
    det_bias = gdsii.extract('MS_det_bias')
    node_list = get_node_list(gdspy.CellReference(mask))
    inv_cell_list = []
    get_cell_list(mask, inv_cell_list)
    # print (len(cell_list))
    for cell in inv_cell_list:
        invert_cell(cell)
    newmask = gdspy.Cell('reticleGNDinv')
    newmask.add(generate_inverted_mask(node_list))
    main_lib.add(newmask)
    main_lib.add(mask)
    mask2_cell_list = []
    get_cell_list(mask2, mask2_cell_list)
    main_lib.add(mask2_cell_list)
    main_lib.add(mask2)
    glob_cell_list = []
    get_cell_list(global_overlay, glob_cell_list)
    main_lib.add(pixel)
    main_lib.add(det_bias)
    main_lib.add(glob_cell_list)
    main_lib.add(global_overlay)

    main_lib.write_gds('corwin_wafer_edited.gds',unit=1e-6,precision=1e-9)   
