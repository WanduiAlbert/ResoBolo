#!/usr/bin/env python3

import gdspy
import patches
import smallscaledarkresonators as ssd
import numpy as np


def_layers = {"Alignment Layer":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "Aluminum":5,
        "Resistor":6, "STEPPER":7, "Something":14, "Nb Wire":9,\
                "XeF2":10, "Island":11, "GP":12, "Wafer":22, "Thin Gold":26, "Anodization":30}
 

# def invert_cell(cell, rotation=0):
#     layers = cell.get_layers()

#     if len(layers) == 1:
#         inverter(cell)

#     for cell in cell.get_dependencies():
#         invert_cell(cell, rotation)
   

if __name__ == "__main__":
    fn = 'corwin_wafer_20180829.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    cell_dict = gdsii.cell_dict
    # bowtie_r = ssd.invert_cell(cell_dict['R_BowtieSlot'])[0]
    # gdsii.add(bowtie_r)
    # gdsii.write_gds('diplexer_ready_for_production.gds',unit=1e-6,precision=1e-9)

    globaloverlay = cell_dict['Module_features_just_tile']
    mask_list = [cell_dict['reticle_1'], cell_dict['reticle_2']]
    to_ignore = set()
    allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
            layer_dict=def_layers, layer_order=None, cellsInverted=False)
    patchtable = patches.PatchTable(allshots, 'corwin_spreadsheet_20180829.xlsx')
    patchtable.generate_spreadsheet()

