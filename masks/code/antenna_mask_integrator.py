#! /usr/bin/env python3

import numpy as np

import gdspy as gds

mask_dir = "../resobolo_files/"
ant_layers = [1, 5, 6, 7, 8,9,11,26]
ant_layer_names = ['ALUMINUM', '120nm_NbWiring', 'ThinGold', 'DC_Wiring', 'Wiring_BUS', '400nm_Niobium', 'FSN', 'Heat_capacity']
ant_layer_mapping = dict(zip(ant_layer_names, ant_layers))

reso_layers = list(range(1,11))
reso_layer_names = ['ThinGold', 'PRO1', 'ALUMINUM', 'LSNSUB', 'LSN1', '120nm_NbWiring', 'STEPPER', '400nm_NbWiring', 'ILD', 'NEW LAYER']
reso_layer_mapping = dict(zip(reso_layer_names, reso_layers))

# We also want to load in the TKID waffle mask file to use as a basis for the calculations here

if __name__ == "__main__":
    ant_maskfile = mask_dir + "270GHz-layout-Sep2016.gds"
    ant_mask = gds.GdsLibrary()
    ant_mask.read_gds(ant_maskfile)

    # print ("\nTop Level Antenna Mask Cells\n")
    # for cell in ant_mask.top_level():
    #     print (cell.name)

    reso_maskfile = mask_dir + "Waffle_TKID_Mask.gds"
    reso_mask = gds.GdsLibrary()
    reso_mask.read_gds(reso_maskfile)

    # print ("\nTKID Mask Cells\n")
    # for cell in reso_mask.cell_dict:
    #     print (cell)

    waffle_islands = reso_mask.extract('LSN_Islands')
    waffle_elems = waffle_islands.elements
    for el in waffle_elems :
        if el.layer == reso_layer_mapping['NEW LAYER'] :
            waffle_elems.remove(el)
    # 
    print (len(waffle_elems))

    # global_overlay = gds.Cell('GLOBAL_overlay')

    # for dep in wafer_layout.get_dependencies():
    #     if (dep.name == 'GLOBAL_overlay'):
    #         global_overlay = dep
    
    ant_mask.add(waffle_islands)

    gds.LayoutViewer()