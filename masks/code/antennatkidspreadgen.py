#!/usr/bin/env python3

import gdspy
import patches
import numpy as np
import smallscaledarkresonators as ssd
import pdb

inv_margin = 200
mask_dir = "../mask_files/"
top = None
# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def_layers = {"Wafer_footprint":1, "Reticle_Outline":2, "GP":3, "CIF_LSNSUB":4, "CIF_ILD":5,
		"Cap_Nb_120nm":6, "Capacitor_Etch":7, "Al":8, "Au_Heater":9, "RES_PRO":10,
		"MS":11, "Au_Thermal_Sink":12, 'CIF_LSN1':13, 'XeF2':14,
        "Filler_btwn_cells":15}
#def_layers = {"Thin Gold":1, "RES_PRO":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
#		"120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10,
#		"Au_Thermal_Sink":11, "GP":12, 'Capacitor_Etch':13, 'Wafer Outline':22}

"""
Useful for generating the overlay composed strictly of all the cells having been
inverted. The use for this is to check for improper overlaps (or lack thereof)
between different cells in the overlay.
"""
def generate_inverted_cell_in_overlay(cellref, mask_components):
	cellname = cellref.ref_cell.name
	invcellname = cellname + "_inv"
	#if cellname == "50umX15mm_Hline_r": pdb.set_trace()
	if invcellname in mask_components:
		invcellref = gdspy.CellReference(invcellname)
	elif cellname in mask_components:
		invcellref = gdspy.CellReference(cellname)
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

"""
Creates the inverted overlay by calling generate_inverted_cell_in_overlay
"""
def generate_inverted_overlay(wafer, mask_components):
	allcells = main_lib.cell_dict
	#pdb.set_trace()
	#masklist = [x for x in allcells if x.endswith('inv')]
	invwafer = gdspy.Cell('Wafer_Layout_new_with_permiter_cells_at_center_Inverted')
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
	print (cell.name)
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
	# Try and collapse the polygons together
	new_polyset = gdspy.fast_boolean(new_polyset, None, 'or', layer=layer)
	inv_cell.add(new_polyset)
	#print (inv_cell)
	return inv_cell

def invert_cell(cell, rotation=0):
	if not cell.get_dependencies():
		return [inverter(cell)]

	icells = []
	for subcell in cell.get_dependencies():
		icells += invert_cell(subcell)
	return icells

def get_mask_lens():
	diameter = 31112
	radius = diameter/2
	lens = gdspy.Round([0,0], radius, number_of_points=2000, layer=def_layers['STEPPER'])
	return gdspy.fast_boolean(lens, None, 'or', layer=def_layers['STEPPER'])

def makechipoutline(width, length, outlinename):
	specs = {'layer':10}
	outline = gdspy.Path(2, (0,0))
	#outline.layers = [0]
	outline.segment(width, '+x', layer=def_layers['PRO1'])
	outline.segment(length, '+y', layer=def_layers['PRO1'])
	outline.segment(width, '-x', layer=def_layers['PRO1'])
	outline.segment(length, '-y', layer=def_layers['PRO1'])
	outline.translate(-width/2, -length/2)
	outlinecell = gdspy.Cell(outlinename)
	outlinecell.add(outline)
	return outlinecell

def get_inverted_cells():
	cells = main_lib.cell_dict
	cellnames = cells.keys()

	return [cells[name] for name in cellnames if name.endswith('_inv') ]

def make_inverted_cells():
	allcells = main_lib.cell_dict
	cellnames = allcells.keys()
	global top
	for cell in top.get_dependencies():
		if not cell.name.endswith('_inv'):
			invert_cell(cell)

def generate_mask():
	print ("Generating the mask....\n\n")
	make_inverted_cells()
	print ("\n\nMask Generation Completed.\n\n")
	return None

def validate_mask():
	print ("Validating the mask....\n\n")
	all_cells = main_lib.cell_dict
	mask = all_cells['reticle_nolens']
	maskhasOverlaps = ssd.check_cell_for_overlaps(mask)
	if maskhasOverlaps:
		print ("FIX ME: Some cells on the mask were found to have overlaps.")
	else:
		print ("No overlaps found on the mask file")
	#filler = ssd.fill_empty_space(mask, ssd.mask_width, ssd.mask_length)
	# A bit of a hack but works
	#filler.layers = [def_layers['Wafer Outline'] for _ in filler.layers]
	#mask.add(filler)

	#canon_lens = get_mask_lens()
	#mask.add(canon_lens)
	#maskoutline = makechipoutline(ssd.mask_width, ssd.mask_length, 'MaskOutline')
	#mask.add(maskoutline)

	print ("\n\nMask Validation Completed.\n\n")
	return mask

if __name__ == "__main__":

	# Generate all the missing inverted cells from the base file
	base_fn = 'Antenna_Coupled_TKIDs_20190221_ROB.gds'
	final_fn = 'Antenna_Coupled_TKIDs_20190221_ROB.gds'
	#main_lib.read_gds(base_fn)
	#cell_dict = main_lib.cell_dict
	#top = cell_dict['Wafer_Layout_new_with_perimeter_cells_at_center']
	#gdspy.write_gds(final_fn, unit=1e-6,precision=1e-9)
	toGenerateMask = False
	makeInvertedOverlay = False
	fillmask = True

	if toGenerateMask:
		main_lib.read_gds(final_fn)
		cell_dict = main_lib.cell_dict
		top = cell_dict['Wafer_Layout_new_with_permiter_cells_at_center']
		mask = generate_mask()
		gdspy.write_gds(final_fn, unit=1e-6,precision=1e-9)
		exit()

	if makeInvertedOverlay:
		main_lib.read_gds(final_fn)
		cell_dict = main_lib.cell_dict

		## Give the name of the top cell in the file just read.
		top = cell_dict['Wafer_Layout_new_with_permiter_cells_at_center']
		masklist = cell_dict['reticle1'].get_dependencies()
		mask_components = {x.name:x for x in masklist}
		inverted_cells = []
		for cell in top.get_dependencies():
			inverted_cells += invert_cell(cell)
		generate_inverted_overlay(top, mask_components)
		#for cell in inverted_cells:
		#	print (cell)
		#gdspy.write_gds(base_fn, unit=1e-6,precision=1e-9)
		gdspy.write_gds(final_fn, unit=1e-6,precision=1e-9)
		if not fillmask: exit()

	"""
	Now work with the inverted cells to finalize the design of the reticle.
	"""
	if fillmask:
		main_lib.read_gds(base_fn)
		# I'll assume that all the necessary inverted cells have been generated.
		# I now want to make the reticle cells
		mask = validate_mask()

		#gdspy.write_gds(final_fn, unit=1e-6,precision=1e-9)
		#exit()
	#main_lib.read_gds(base_fn)

	"""
	Now we can generate the spreadsheet for the antenna coupled device
	"""
	main_lib.read_gds(final_fn)
	cell_dict = main_lib.cell_dict

	# Give the name of the global overlay cell in the file just read.
	globaloverlay = cell_dict['Wafer_Layout_new_with_permiter_cells_at_center']
	# Specify all the mask files to be referenced as a list.
	mask_list = [cell_dict['reticle1']]
	# specify any cells in the global overlay that should be ignored because
	# they are not on the reticles.
	to_ignore = set()
	# A dictionary mapping the layer names to the GDSII layer numbers.
	# Specify the order in which the layers should be arranged in the
	# spreadsheet.
	#layer_order = [3, 12, 9, 6, 13, 4, 5, 2, 8, 1, 11, 10 ]
	layer_order = list(range(1, 16))
    # This command creates a bunch of shots by running through the global
	# overlay and the reticle list. If the cells on the mask are inverted
	# specify so.
	allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
			layer_dict=def_layers, layer_order=layer_order, cellsInverted=True)
	# Create a patch table object from the list of shots generated.
	patchtable = patches.PatchTable(allshots, 'antenna_coupled_TKIDs_20190221.xlsx')
	patchtable.generate_spreadsheet()


