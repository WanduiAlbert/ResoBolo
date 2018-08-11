from __future__ import division

import gdspy
import numpy as np
import openpyxl

fn_spreadsheet = 'ResonatorArray.xlsx'
fn_mask = 'sscale_darkres.gds'
fn_out = 'sscale_darkres_simulated.gds'

mask_cell = 'ResoArray_Mask_May2018'

print("Loading mask...")
lib = gdspy.GdsLibrary()
lib.read_gds(fn_mask)
pmask = lib.cell_dict[mask_cell].flatten().elements
c = gdspy.Cell(mask_cell)
c.add(pmask)

print("Loading spreadsheet...")
wb = openpyxl.load_workbook(fn_spreadsheet)
ws = wb[wb.sheetnames[0]]
r = 15

next_layer_num = 1
#layer_nums = {}
layer_nums = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,\
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9,\
        "XeF2":10, "GP":12, 'Wafer Outline':22}

cwafer = gdspy.Cell('wafer')

while True:
	
	layer = ws['A' + str(r)].value
	
	if not layer:
		break
	
	name = "row" + str(r) + '_' + ws['C' + str(r)].value + '_r'
	print("Working on cell " + name)
	
	if layer not in layer_nums:
		layer_nums[layer] = next_layer_num
		next_layer_num += 1
		
	b_xl = 1e3*float(ws['H' + str(r)].value)
	b_xr = 1e3*float(ws['I' + str(r)].value)
	b_yu = 1e3*float(ws['J' + str(r)].value)
	b_yd = 1e3*float(ws['K' + str(r)].value)
	ws_x = 1e3*float(ws['O' + str(r)].value)
	ws_y = 1e3*float(ws['P' + str(r)].value)
	ar_nx = int(ws['T' + str(r)].value)
	ar_ny = int(ws['U' + str(r)].value)
	ar_cx = 1e3*float(ws['V' + str(r)].value)
	ar_cy = 1e3*float(ws['W' + str(r)].value)
	ar_xp = 1e3*float(ws['X' + str(r)].value)
	ar_yp = 1e3*float(ws['Y' + str(r)].value)
	
	pblades = gdspy.Rectangle((b_xl,b_yu),(b_xr,b_yd))
	
	c = gdspy.Cell("blades_" + name)
	#c.add(gdspy.CellReference('mask'))
	c.add(pblades)
	
	pcell = gdspy.fast_boolean(pmask, pblades, 'and', layer=layer_nums[layer])
	assert pcell is not None

	c = gdspy.Cell(name)
	c.add(pcell)
	
	x0 = ws_x + ar_cx - (ar_nx - 1) * ar_xp / 2.0
	y0 = ws_y + ar_cy - (ar_ny - 1) * ar_yp / 2.0
	
	for x in range(ar_nx):
		for y in range(ar_ny):
			cwafer.add(gdspy.CellReference(name, (x0 + x*ar_xp, y0 + y*ar_yp)))
	
	r += 1
	
print("Writing file...")
#gdspy.LayoutViewer()
import pdb
pdb.set_trace()
gdspy.write_gds(fn_out, unit=1.0e-6, precision=1.0e-9)
print("Done!")
