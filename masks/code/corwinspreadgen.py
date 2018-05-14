#!/usr/bin/env python3

import gdspy
import patches
import numpy as np


def_layers = {"Alignment Layer":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "Aluminum":5,
        "Resistor":6, "STEPPER":7, "Something":14, "Nb Wire":9,\
                "XeF2":10, "Island":11, "GP":12, "Wafer":22, "Thin Gold":26}


if __name__ == "__main__":
    fn = 'corwin_wafer_edited.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    cell_dict = gdsii.cell_dict

    globaloverlay = cell_dict['wafer']
    mask_list = [cell_dict['reticleMS'], cell_dict['reticleGND']]
    to_ignore = set('frame')
    allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
            layer_dict=def_layers, layer_order=None, cellsInverted=False)
    patchtable = patches.PatchTable(allshots, 'corwinwafer.xlsx')
    patchtable.generate_spreadsheet()

