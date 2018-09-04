#!/usr/bin/env python3

import gdspy
import patches
import smallscaledarkresonators as ssd
import numpy as np

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main')
gdspy.current_library = main_lib

def_layers = {"Al_TES":1, "Outline":2, "PRO1":3, "Ti_TES":4, "GP":5, "PRO2":6,\
    "ILD_Via":7, "RES":8, "MS":9, "RES_pro":10, "heat_cap":11, "pic_frame":12,\
    "nitride":13, "XeF2_STS":14, "New Layer":15, "Lens":20, "Die cutout":22,\
    "Anodization":30}

#def_layers = {"Alignment Layer":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "Aluminum":5,
#        "Resistor":6, "STEPPER":7, "Something":14, "Nb Wire":9,\
#                "XeF2":10, "Island":11, "GP":12, "Wafer":22, "Thin Gold":26, "Anodization":30}
def make_inverted_cell(cellref, mask_components):
    cellname = cellref.ref_cell.name
    # if cellname == "terminal_lines_and_pad":
    #     pdb.set_trace()
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
#     layers = cell.get_layers()

#     if len(layers) == 1:
#         inverter(cell)

#     for cell in cell.get_dependencies():
#         invert_cell(cell, rotation)


if __name__ == "__main__":
    fn = 'diplexer_FINAL_reformating.gds'
    main_lib.read_gds(fn)
    cell_dict = main_lib.cell_dict
    # bowtie_r = ssd.invert_cell(cell_dict['R_BowtieSlot'])[0]
    # gdsii.add(bowtie_r)
    # gdsii.write_gds('diplexer_ready_for_production.gds',unit=1e-6,precision=1e-9)

    globaloverlay = cell_dict['Module_features_just_tile']
    mask_list = [cell_dict['reticle_1'], cell_dict['reticle_2']]
    to_ignore = set('frame')
    allshots = patches.gen_patches_table(globaloverlay, mask_list, to_ignore,\
            layer_dict=def_layers, layer_order=None, cellsInverted=False)
    patchtable = patches.PatchTable(allshots, 'diplexer_FINAL_reformating.xlsx')
    patchtable.generate_spreadsheet()

    invwafer = generate_inverted_overlay(globaloverlay, mask_list)

    main_lib.write_gds(fn,unit=1e-6,precision=1e-9)
