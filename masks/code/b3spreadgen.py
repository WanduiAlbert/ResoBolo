#!/usr/bin/env python3

import gdspy
import patches
import numpy as np


def_layers = {"Al TES":1, "PRO1":2, "TES":3, "PRO2":4, "GP":5,\
        "RES":6, "SiO_PRO":7, "DC_wire":8, "microstrip":9,"Layer 1":15,\
        "XeF2":10, "FSN":11, "Expose":21, "bling_temp":26,\
        "Via":27, "Anodization":12,"Backshort":13,"Au pictureframe":14}


if __name__ == "__main__":
    fn = 'resobolo_files/B3_v2_Jul13.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    cell_dict = gdsii.cell_dict

    globaloverlay = cell_dict['Overlay_Full_tile']
    mask_list = [cell_dict['Top_Mask_B3_Jul13']]
    to_ignore = set()
    allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
            layer_dict=def_layers, layer_order=None, cellsInverted=False)
    patchtable = patches.PatchTable(allshots, 'b3.xlsx')
    patchtable.generate_spreadsheet()

