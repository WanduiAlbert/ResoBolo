#!/usr/bin/env python3

import gdspy
import patches
import numpy as np
import smallscaledarkresonators as ssd
import pdb

inv_margin = 200
mask_dir = "/home/wanduialbert/Desktop/research_stuff/resobolo/masks/mask_files/"

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
		"120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12, 'Wafer Outline':22}

def generate_inverted_cell_in_overlay(cellref, mask_components):
	cellname = cellref.ref_cell.name
	invcellname = cellname + "_inv"
	if invcellname in mask_components:
		invcellref = gdspy.CellReference(invcellname)
	else:
		try:
			invcell = main_lib.cell_dict[invcellname]
			invcellref = gdspy.CellReference(invcell)
		except KeyError:
			invcell = gdspy.Cell(invcellname)
			if cellref.ref_cell.get_dependencies():
				subcells = cellref.ref_cell.elements
			else:
				subcells = [cellref]
			subcells = [x for x in subcells if type(x) in
					patches.allowed_element_types]
			for a_cell in subcells:
				a_cellref = generate_inverted_cell_in_overlay(a_cell, mask_components)
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

def generate_inverted_overlay(wafer, mask_components):
	allcells = main_lib.cell_dict
	masklist = [x for x in allcells if x.endswith('inv')]
	mask_components = set(masklist)
	invwafer = gdspy.Cell('Wafer_Layout_Inverted')
	gcomponents = wafer.elements
	#pdb.set_trace()
	for component in gcomponents:
		if type(component) not in patches.allowed_element_types: continue
		invcomponent = generate_inverted_cell_in_overlay(component, mask_components)

		if type(component) == gdspy.CellArray:
			# continue
			invcomponent = gdspy.CellArray(invcomponent.ref_cell, columns=component.columns, rows=component.rows,\
			spacing=component.spacing)
		invcomponent.origin = component.origin
		invwafer.add(invcomponent)


	return invwafer

def cellIsPresent(cellname):
	allcellnames = main_lib.cell_dict.keys()
	return cellname in allcellnames

def inverter(cell, rotation=0):
	cell_name = cell.name + '_inv'
	if cellIsPresent(cell_name): return
	inv_cell = gdspy.Cell(cell_name)
	cell = cell.flatten()
	cell_ref = gdspy.CellReference(cell, rotation=rotation)
	dx, dy = ssd.get_size(cell_ref)
	dx += 2*inv_margin
	dy += 2*inv_margin
	dx = ssd.roundto(dx, 100)
	dy = ssd.roundto(dy, 100)
	layer = cell.get_layers().pop()
	polys = cell_ref.get_polygons(depth=1)
	polyset = gdspy.PolygonSet(polys, layer=layer)
	bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)
	new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
	inv_cell.add(new_polyset)
	#print (inv_cell)
	return inv_cell

def invert_cell(cell, rotation=0):
	layers = cell.get_layers()

	if len(layers) == 1:
		return [inverter(cell)]

	icells = []
	for subcell in cell.get_dependencies():
		icells += invert_cell(subcell)
	return icells

if __name__ == "__main__":
	# File name to be read. Must be in GDS format.
	fn = 'one_pixel.gds'
	main_lib.read_gds(mask_dir + fn)
	cell_dict = main_lib.cell_dict

	# Give the name of the global overlay cell in the file just read.
	globaloverlay = cell_dict['GLOBAL_Overlay_Export']
	top = cell_dict['Wafer_Layout']
	inverted_cells = []
	for cell in top.get_dependencies():
		inverted_cells += invert_cell(cell)

	inverted_cells = [icell for icell in inverted_cells if icell is not None]
	generate_inverted_overlay(top, set())
	#for cell in inverted_cells:
	#	print (cell)
	gdspy.write_gds(fn, unit=1e-6,precision=1e-9)
