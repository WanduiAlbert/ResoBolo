#!/usr/bin/env python3

import gdspy
import patches
import smallscaledarkresonators as ssd
import numpy as np

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
		"120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12, 'Wafer Outline':22}

def generate_inverted_overlay(wafer, mask_list):
	mask_components = set()
	for mask in mask_list:
		for x in mask.elements:
			if type(x) not in patches.allowed_element_types: continue
			mask_components.add(x.ref_cell.name)
	invwafer = gdspy.Cell('Global_Overlay_4in_Inverted')
	gcomponents = wafer.elements
	for component in gcomponents:
		# print (component)
		ctype = type(component)
		if ctype not in patches.allowed_element_types: continue
		if component.ref_cell.name == "frame" : continue
		try:
			invcomponent = ssd.make_inverted_cell(component, mask_components)
		except ValueError:
			invcomponent = component
		if ctype == gdspy.CellArray:
			invcomponent = gdspy.CellArray(invcomponent.ref_cell,
					columns=component.columns, rows=component.rows,
					spacing=component.spacing)
		invcomponent.origin = component.origin
		invwafer.add(invcomponent)


	return invwafer

if __name__ == "__main__":
	fn = '../mask_files/sscale_darkres_4in.gds'
	main_lib.read_gds(fn)
	cell_dict = main_lib.cell_dict
	# bowtie_r = ssd.invert_cell(cell_dict['R_BowtieSlot'])[0]
	# gdsii.add(bowtie_r)
	# gdsii.write_gds('diplexer_ready_for_production.gds',unit=1e-6,precision=1e-9)

	globaloverlay = cell_dict['Global_Overlay_4in']
	mask_list = [cell_dict['ResoArray_Mask_May2018']]
	to_ignore = set(['WaferOutline'])
	layer_order = [12, 9, 6, 3, 8 ]
	allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,
			def_layers, layer_order)
	patchtable = patches.PatchTable(allshots, 'ResonatorArray_4in.xlsx')
	patchtable.generate_spreadsheet()
	invwafer = generate_inverted_overlay(globaloverlay, mask_list)
	gdspy.write_gds(fn,unit=1e-6,precision=1e-9)
