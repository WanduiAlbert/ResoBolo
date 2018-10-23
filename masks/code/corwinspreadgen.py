#!/usr/bin/env python3

import gdspy
import patches
import smallscaledarkresonators as ssd
import numpy as np

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def_layers = {"Al_TES":1, "Outline":2, "PRO1":3, "Ti_TES":4, "GP":5, "PRO2":6,\
	"VIA":7, "RES":8, "MS":9, "RES_pro":10, "heat_cap":11, "pic_frame":12,\
	"LSN":13, "XeF2_STS":14, "New Layer":15, "Lens":20, "Die cutout":22,\
	"Anodization":30}

#def_layers = {"Alignment Layer":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "Aluminum":5,
#		"Resistor":6, "STEPPER":7, "Something":14, "Nb Wire":9,\
#				"XeF2":10, "Island":11, "GP":12, "Wafer":22, "Thin Gold":26, "Anodization":30}
def make_inverted_cell(cellref, mask_components):
	cellname = cellref.ref_cell.name
	# if cellname == "terminal_lines_and_pad":
	#	 pdb.set_trace()
	invcellname = cellname
	if invcellname in mask_components:
		invcellref = gdspy.CellReference(invcellname)
	else:
		try:
			invcell = main_lib.cell_dict[invcellname]
			invcellref = gdspy.CellReference(invcell)
		except KeyError:
			invcell = gdspy.Cell(invcellname)
			subcells = cellref.ref_cell.elements
			for a_cell in subcells:
				a_cellref = make_inverted_cell(a_cell, mask_components)
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

def generate_inverted_overlay(wafer, mask_list):
	mask_components = set()
	for mask in mask_list:
		for x in mask.elements:
			if type(x) not in patches.allowed_element_types: continue
			mask_components.add(x.ref_cell.name)
	#print (mask_components)
	# invcomponent = make_inverted_cell(component, mask_components)
	invwafer = gdspy.Cell('Global_Overlay_Inverted')
	gcomponents = wafer.elements
	for component in gcomponents:
		# print (component)
		ctype = type(component)
		if ctype not in patches.allowed_element_types: continue
		if component.ref_cell.name == "frame" : continue
		invcomponent = make_inverted_cell(component, mask_components)

		if ctype == gdspy.CellArray:
			# continue
			invcomponent = gdspy.CellArray(invcomponent.ref_cell, columns=component.columns, rows=component.rows,\
			spacing=component.spacing)
		invcomponent.origin = component.origin
		invwafer.add(invcomponent)


	return invwafer


# def invert_cell(cell, rotation=0):
#	 layers = cell.get_layers()

#	 if len(layers) == 1:
#		 inverter(cell)

#	 for cell in cell.get_dependencies():
#		 invert_cell(cell, rotation)


if __name__ == "__main__":
	# File name to be read. Must be in GDS format.
	fn = 'diplexer_FINAL_AT_more_edits.gds'
	main_lib.read_gds(fn)
	cell_dict = main_lib.cell_dict

	# Give the name of the global overlay cell in the file just read.
	globaloverlay = cell_dict['Module_features_just_tile']
	# Specify all the mask files to be referenced as a list.
	mask_list = [cell_dict['reticle_1'], cell_dict['reticle_2']]
	# specify any cells in the global overlay that should be ignored because
	# they are not on the reticles.
	to_ignore = set(['frame'])
	# A dictionary mapping the layer names to the GDSII layer numbers.
	def_layers = {"Al_TES":1, "Outline":2, "PRO1":3, "Ti_TES":4, "GP":5, "PRO2":6,\
	"VIA":7, "RES":8, "MS":9, "RES_pro":10, "heat_cap":11, "pic_frame":12,\
	"LSN":13, "XeF2_STS":14, "New Layer":15, "Lens":20, "Die cutout":22,\
	"Anodization":30}
	# Specify the order in which the layers should be arranged in the
	# spreadsheet.
	layer_order = [1,3,4,6,5,8,10,7,9,13,11,12,14,30,2,15,20,22]
	# This command creates a bunch of shots by running through the global
	# overlay and the reticle list. If the cells on the mask are inverted
	# specify so.
	allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
			layer_dict=def_layers, layer_order=layer_order, cellsInverted=False)
	# Create a patch table object from the list of shots generated.
	patchtable = patches.PatchTable(allshots, 'diplexer_FINAL_AT_more_edits.xlsx')
	patchtable.generate_spreadsheet()

