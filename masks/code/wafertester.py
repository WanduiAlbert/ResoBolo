#!/usr/bin/env python3

import gdspy
import patches
import numpy as np

def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12}
layer_order = [12, 9, 6, 3, 8 ]


if __name__ == "__main__":
    fn = 'sscale_darkres.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    cell_dict = gdsii.cell_dict

    globaloverlay = cell_dict['wafer']
    mask_list = [cell_dict['reticleMS'], cell_dict['reticleGND'], cell_dict['reticleGNDinv']]
    to_ignore = set('frame')
    allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
            layer_dict=def_layers, layer_order=None, cellsInverted=False)
    patchtable = patches.PatchTable(allshots, 'corwinwafer.xlsx')
    patchtable.generate_spreadsheet()


