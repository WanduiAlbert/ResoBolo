#!/usr/bin/env python3

import gdspy
import patches
import smallscaledarkresonators as ssd
import numpy as np
import pdb

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
		"120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12, 'Wafer Outline':22}

def make_inverted_cell(cellref, mask_components):
	cellname = cellref.ref_cell.name
	print (cellname)
	#if cellname == "reso_structure_4in": pdb.set_trace()
	invcellname = ssd.get_inv_cellname(cellname)
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
	#for comp in mask_components:
	#	print (comp)
	#exit()
	invwafer = gdspy.Cell('Global_Overlay_4in_Inverted')
	gcomponents = wafer.elements
	for component in gcomponents:
		# print (component)
		ctype = type(component)
		if ctype not in patches.allowed_element_types: continue
		if component.ref_cell.name == "frame" : continue
		try:
			invcomponent = make_inverted_cell(component, mask_components)
		except ValueError:
			invcomponent = component
		if ctype == gdspy.CellArray:
			invcomponent = gdspy.CellArray(invcomponent.ref_cell,
					columns=component.columns, rows=component.rows,
					spacing=component.spacing)
		invcomponent.origin = component.origin
		invwafer.add(invcomponent)


	return invwafer


possible_endings = ['_r', '_inv']
# Finds a suitable name for a cell from the second gds file that does not
# conflict with a name in the 1st file.
def rename(name, ending, index=2):
	if not index: index = ''
	newname = name
	for an_ending in possible_endings:
		if name.endswith(an_ending):
			newname = name.replace(an_ending, ending)
			break
	if not newname.endswith(ending):
		newname += str(index)
	else:
		newname = newname.strip(ending) + str(index) + ending
	return newname

#def inv_mask_cellname(shot, ending):
#	if shot.cellname.endswith('_r'):
#		name = shot.cellname
#	elif shot.cellname.endswith('inv'):
#		name = shot.cellname
#	elif shot.cellname.endswith('-r'):
#		name = shot.cellname
#	else:
#		name = shot.cellname + ending
#	return name

def get_rename_dict(mainlibcells, maskcells, fn2, ending, index=2):
	tmplib = gdspy.GdsLibrary(name='tmp')
	tmplib.read_gds(fn2)
	tmplibcells = list(tmplib.cell_dict.keys())

	renamer = {}
	# Simply rename all the cells in tmp so they are easily identifiable
	for tmpcell in tmplibcells:
		if tmpcell in mainlibcells or tmpcell.startswith('alignment') or tmpcell.startswith('50um'):
			renamer[tmpcell] = rename(tmpcell, ending, index)
		else:
			renamer[tmpcell] = rename(tmpcell, ending, index=None)
			#if newname not in mainlibcells:
			#	renamer[tmpcell] = newname
			#else:
			#	renamer[tmpcell] = rename(tmpcell, ending, index)


	layermap = {}
	# Only want to check single layer cells because those are represented on the
	# reticle
	maintoconsider = list(filter(
		lambda x: len(main_lib.cell_dict[x].get_layers()) == 1,
		mainlibcells))
	tmptoconsider = list(filter(
		lambda x: len(tmplib.cell_dict[x].get_layers()) == 1,
		tmplibcells))
	tmptoconsidermatched = [x.replace('_inv', ending) for x in tmptoconsider]
	tmpmap = dict(zip(tmptoconsidermatched, tmptoconsider))
	#for tlcell in tmplibcells:
	#	renamer[tlcell] = tlcell.replace('_inv', ending)
	for mlcell in maintoconsider:
		isinv = mlcell.endswith(ending)
		invmlcell = mlcell + ending if not mlcell.endswith(ending) else mlcell
		if mlcell in maskcells or invmlcell in maskcells: continue
		layer_1st, = main_lib.cell_dict[mlcell].get_layers()
		layer_2nd = None
		if mlcell in tmptoconsider:
			print (mlcell)
			renamer[mlcell] = mlcell
			layer_2nd, = tmplib.cell_dict[mlcell].get_layers()
		elif not isinv and invmlcell in tmptoconsidermatched:
			print (mlcell)
			renamer[tmpmap[invmlcell]] = invmlcell
			layer_2nd, = tmplib.cell_dict[tmpmap[invmlcell]].get_layers()
		if layer_2nd and layer_2nd not in layermap:
			layermap[layer_2nd] = layer_1st
		#if mlcell in maskcells or mlcell + ending in maskcells:
		#	renamer[mlcell] = rename(mlcell, ending, index)
		#	#print (mlcell, renamer[mlcell])
		#elif mlcell in tmplibcellsmatched:
		#	renamer[tmpmap[mlcell]] = rename(mlcell, ending, index)
		#	#print (mlcell, tmpmap[mlcell], renamer[tmpmap[mlcell]])
		#else:
		#layermap{tmplib.cell_dict[mlcell].get_layers()[0]:mlcell.get_layers()[0]}

	print (layermap)
	return renamer, layermap

if __name__ == "__main__":
	invcell_ending = '_r'
	fn = '../mask_files/sscale_darkres_4in_with_island_newfeedline_20190912_AW.gds'
	main_lib.read_gds(fn)
	mainlibcells = main_lib.cell_dict.keys()
	cellsonmask = [x.name for x in main_lib.cell_dict['ResoArray_Mask_Sept2019'].get_dependencies(recursive=True)]
	fn2 = 'Antenna_Coupled_TKIDs_20190826_AW.gds'
	renamer, layermap = get_rename_dict(mainlibcells, cellsonmask, fn2, invcell_ending)
	main_lib.read_gds(fn2, rename=renamer, layers=layermap)
	cell_dict = main_lib.cell_dict
	#for k in sorted(cell_dict.keys()):
	#	print (k)
	#exit()
	globaloverlay = cell_dict['Global_Overlay_4in']
	mask_list = [cell_dict['ResoArray_Mask_Sept2019'], cell_dict['reticle1']]
	to_ignore = set(['WaferOutline'])
	layer_order = [12, 9, 4, 6, 5, 3, 8, 10]
	#pdb.set_trace()
	allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,
			def_layers, layer_order, invcell_ending=invcell_ending)
	patchtable = patches.PatchTable(allshots,
			'ResonatorArray_4in_newfeedline_20190924.xlsx')
	patchtable.generate_spreadsheet()
	invwafer = generate_inverted_overlay(globaloverlay, mask_list)
	gdspy.write_gds(fn,unit=1e-6,precision=1e-9)
