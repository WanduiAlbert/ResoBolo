#!/usr/bin/env python3

import gdspy
import patches
import numpy as np
import smallscaledarkresonators as ssd
import pdb

inv_margin = 200
mask_dir = "../mask_files/"

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10,
        "Au_Thermal_Sink":11, "GP":12, 'Wafer Outline':22}

def generate_inverted_cell_in_overlay(cellref, mask_components):
    cellname = cellref.ref_cell.name
    invcellname = cellname + "_inv"
    if invcellname in mask_components:
        invcellref = gdspy.CellReference(invcellname)
    else:
        try:
            invcell = main_lib.cell_dict[invcellname]
            invcellref = gdspy.CellReference(invcell)
        except KeyError:
            invcell = gdspy.Cell(invcellname)
            if cellref.ref_cell.get_dependencies():
                subcells = cellref.ref_cell.elements
            else:
                subcells = [cellref]
            subcells = [x for x in subcells if type(x) in
                    patches.allowed_element_types]
            for a_cell in subcells:
                a_cellref = generate_inverted_cell_in_overlay(a_cell, mask_components)
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

def generate_inverted_overlay(wafer, mask_components):
    allcells = main_lib.cell_dict
    masklist = [x for x in allcells if x.endswith('inv')]
    mask_components = set(masklist)
    invwafer = gdspy.Cell('Wafer_Layout_Inverted')
    gcomponents = wafer.elements
    #pdb.set_trace()
    for component in gcomponents:
        if type(component) not in patches.allowed_element_types: continue
        invcomponent = generate_inverted_cell_in_overlay(component, mask_components)

        if type(component) == gdspy.CellArray:
            # continue
            invcomponent = gdspy.CellArray(invcomponent.ref_cell, columns=component.columns, rows=component.rows,\
            spacing=component.spacing)
        invcomponent.origin = component.origin
        invwafer.add(invcomponent)


    return invwafer

def cellIsPresent(cellname):
    allcellnames = main_lib.cell_dict.keys()
    return cellname in allcellnames

def inverter(cell, rotation=0):
    cell_name = cell.name + '_inv'
    if cellIsPresent(cell_name): return
    inv_cell = gdspy.Cell(cell_name)
    cell = cell.flatten()
    cell_ref = gdspy.CellReference(cell, rotation=rotation)
    dx, dy = ssd.get_size(cell_ref)
    dx += 2*inv_margin
    dy += 2*inv_margin
    dx = ssd.roundto(dx, 100)
    dy = ssd.roundto(dy, 100)
    layer = cell.get_layers().pop()
    polys = cell_ref.get_polygons(depth=1)
    polyset = gdspy.PolygonSet(polys, layer=layer)
    bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)
    new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
    inv_cell.add(new_polyset)
    #print (inv_cell)
    return inv_cell

def invert_cell(cell, rotation=0):
    layers = cell.get_layers()

    if len(layers) == 1:
        return [inverter(cell)]

    icells = []
    for subcell in cell.get_dependencies():
        icells += invert_cell(subcell)
    return icells

def get_mask_lens():
    diameter = 31112
    radius = diameter/2
    lens = gdspy.Round([0,0], radius, number_of_points=2000, layer=def_layers['STEPPER'])
    return gdspy.fast_boolean(lens, None, 'or', layer=def_layers['STEPPER'])

def makechipoutline(width, length, outlinename):
    specs = {'layer':10}
    outline = gdspy.Path(2, (0,0))
    #outline.layers = [0]
    outline.segment(width, '+x', layer=def_layers['PRO1'])
    outline.segment(length, '+y', layer=def_layers['PRO1'])
    outline.segment(width, '-x', layer=def_layers['PRO1'])
    outline.segment(length, '-y', layer=def_layers['PRO1'])
    outline.translate(-width/2, -length/2)
    outlinecell = gdspy.Cell(outlinename)
    outlinecell.add(outline)
    return outlinecell

def get_inverted_cells():
    cells = main_lib.cell_dict
    cellnames = cells.keys()

    return [cells[name] for name in cellnames if name.endswith('_inv') ]

def make_inverted_cells():
    allcells = main_lib.cell_dict
    cellnames = allcells.keys()
    # print (cellnames)
    top = allcells['Wafer_Layout']
    for cell in top.get_dependencies():
        if not cell.name.endswith('_inv'):
            invert_cell(cell)

def generate_mask():
    print ("Generating the mask....\n\n")
    all_cells = main_lib.cell_dict
    mask = gdspy.Cell('reticle1')
    maskoutline = makechipoutline(ssd.mask_width, ssd.mask_length, 'MaskOutline')
    make_inverted_cells()

    inv_cell_list = get_inverted_cells() #+ [ixef2, i_island]
    not_yet = set(inv_cell_list)
    total_mask_area = mask_length * mask_width
    total_area_needed = np.sum(list(map(lambda x: get_cell_area(x),\
            inv_cell_list)))


    default_spacing=200 # A suitable spacing between the cells on the reticle
    intercap_spacing=50

    # 1. Start with the largest feedline structures
    ivmain = all_cells['vert_main_with_corners_r']
    ivgndsub = all_cells['vert_gndsub_with_corners_r']
    ihgndtopad = all_cells['gndsub_hor_feedline_to_pad_r']
    ihmaintopad = all_cells['main_hor_feedline_to_pad_r']
    ivgndtopad = all_cells['gndsub_vert_feedline_to_pad_r']
    ivmaintopad = all_cells['main_vert_feedline_to_pad_r']
    ivgndopening = all_cells['GP_edge_opening_vert_r']
    ihgndopening = all_cells['GP_edge_opening_hor_r']
    ivgndfiller = all_cells['GP_edge_filler_vert_r']
    ihgndfiller = all_cells['GP_edge_filler_hor_r']
    ihedge = all_cells['50umX15mm_Hline_r']
    ivedge = all_cells['50umX15mm_Vline_r']

    ivmain_ref = gdspy.CellReference(ivmain)
    ivgndsub_ref = gdspy.CellReference(ivgndsub)
    ihgndtopad_ref = gdspy.CellReference(ihgndtopad)
    ihmaintopad_ref = gdspy.CellReference(ihmaintopad)
    ivgndtopad_ref = gdspy.CellReference(ivgndtopad)
    ivmaintopad_ref = gdspy.CellReference(ivmaintopad)
    ivgndopening_ref = gdspy.CellReference(ivgndopening)
    ihgndopening_ref = gdspy.CellReference(ihgndopening)
    ivgndfiller_ref = gdspy.CellReference(ivgndfiller)
    ihgndfiller_ref = gdspy.CellReference(ihgndfiller)
    ihedge_ref = gdspy.CellReference(ihedge)
    ivedge_ref = gdspy.CellReference(ivedge)

    go_dx, go_dy = get_size(ihgndopening_ref)
    ihgndopening_ref.translate(0, -mask_length/2 + go_dy/2 + default_spacing)
    centerx(ihgndopening_ref, ihgndfiller_ref)
    centerx(ihgndopening_ref, ihgndtopad_ref)
    centerx(ihgndopening_ref, ihmaintopad_ref)
    centerx(ihgndopening_ref, ihedge_ref)
    moveabove(ihgndopening_ref, ihgndfiller_ref, spacing=-default_spacing)

    e_dx, e_dy = get_size(ihedge_ref)
    ihedge_ref.translate(0, mask_length/2 - e_dy/2 - default_spacing)

    centery(ivgndopening_ref, ivmain_ref)
    centery(ivgndopening_ref, ivmaintopad_ref)
    centery(ivgndopening_ref, ivgndtopad_ref)
    centery(ivgndopening_ref, ivedge_ref)
    centery(ivgndopening_ref, ivgndsub_ref)

    mp_dx, mp_dy = get_size(ivgndtopad_ref)
    ivedge_ref.translate(-mask_width/2 + mp_dx/2 + default_spacing, 0)
    moveright(ivedge_ref, ivmain_ref, spacing=-default_spacing)
    moveright(ivmain_ref, ivgndsub_ref, spacing=-default_spacing)
    moveright(ivgndsub_ref, ivgndtopad_ref, spacing=-default_spacing)
    moveright(ivgndtopad_ref, ivmaintopad_ref, spacing=-default_spacing)
    moveabove(ihgndfiller_ref, ivmaintopad_ref, spacing=-2*default_spacing)

    mp_dx, mp_dy = get_size(ivgndfiller_ref)
    ivgndfiller_ref.translate(mask_width/2 - mp_dx/2 - default_spacing, 0)
    # moveleft(ivmaintopad_ref, ivgndfiller_ref, spacing=default_spacing)
    moveleft(ivgndfiller_ref, ivgndopening_ref, spacing=default_spacing)

    # Placement of the common pixel on the mask
    icommoncap = all_cells['Capacitor_common_r']
    i_ind = all_cells['Al_inductor_r']
    icap2ind = all_cells['Cap_to_Ind_lines_r']

    icommoncap_ref = gdspy.CellReference(icommoncap)
    i_ind_ref = gdspy.CellReference(i_ind)
    icap2ind_ref = gdspy.CellReference(icap2ind)

    i_dx, i_dy = get_size(i_ind_ref)
    moveright(icap2ind_ref, i_ind_ref, spacing=-intercap_spacing)
    moveabove(icap2ind_ref, i_ind_ref, spacing=i_dy)
    moveleft(icap2ind_ref, icommoncap_ref, spacing=intercap_spacing)

    mask.add(icap2ind_ref)
    mask.add(icommoncap_ref)
    mask.add(i_ind_ref)


    icommon_pixel = set([icommoncap, i_ind, icap2ind])
    not_yet -= icommon_pixel


    # I want to lay the resonator structure cells at the margins of the feedlines
    # Place the resonator structure features
    iresogndsub = all_cells['reso_GP_sub_r']
    iILDsub = all_cells['reso_ILD_sub_r']
    imainfeed = all_cells['feedline_main_r']
    igndsubfeed = all_cells['feedline_GP_sub_r']
    i2feed = all_cells['cap_to_feed_r']
    i2gnd = all_cells['cap_to_gnd_r']
    ivia = all_cells['Via_to_Ground_r']
    iam = all_cells['alignment_marks_patch_new_r']

    iresogndsub_ref = gdspy.CellReference(iresogndsub)
    iILDsub_ref = gdspy.CellReference(iILDsub)
    imainfeed_ref = gdspy.CellReference(imainfeed)
    igndsubfeed_ref = gdspy.CellReference(igndsubfeed)
    i2feed_ref = gdspy.CellReference(i2feed)
    i2gnd_ref = gdspy.CellReference(i2gnd)
    ivia_ref = gdspy.CellReference(ivia)
    iam_ref = gdspy.CellReference(iam)

    ireso_struct = set([iresogndsub, iILDsub, imainfeed,\
      igndsubfeed, ivia, iam])

    movebelow(i_ind_ref, iam_ref, spacing=intercap_spacing)
    moveright(icap2ind_ref, iam_ref, spacing=-intercap_spacing)
    moveright(i_ind_ref, ivia_ref, spacing=-default_spacing)
    centery(i_ind_ref, ivia_ref)
    moveleft(ivgndopening_ref, iresogndsub_ref, spacing=default_spacing)
    centery(ivgndopening_ref, iresogndsub_ref)
    moveleft(ivgndopening_ref, iILDsub_ref, spacing=default_spacing)
    movebelow(iresogndsub_ref, iILDsub_ref, spacing=default_spacing)

    movebelow(ihedge_ref, imainfeed_ref, spacing=default_spacing/2)
    movebelow(imainfeed_ref, igndsubfeed_ref, spacing=default_spacing)
    movebelow(igndsubfeed_ref, ihmaintopad_ref, spacing=default_spacing)
    movebelow(ihmaintopad_ref, ihgndtopad_ref, spacing=default_spacing/2)



    mask.add(iresogndsub_ref)
    mask.add(iILDsub_ref)
    mask.add(imainfeed_ref)
    mask.add(igndsubfeed_ref)
    mask.add(ivia_ref)
    mask.add(iam_ref)

    not_yet -= ireso_struct

    mask.add(ivmain_ref)
    mask.add(ivgndsub_ref)
    mask.add(ihgndtopad_ref)
    mask.add(ihmaintopad_ref)
    mask.add(ivgndtopad_ref)
    mask.add(ivmaintopad_ref)
    mask.add(ivgndopening_ref)
    mask.add(ihgndopening_ref)
    mask.add(ihgndfiller_ref)
    mask.add(ivgndfiller_ref)
    mask.add(ihedge_ref)
    mask.add(ivedge_ref)

    not_yet -= set([ivmain, ivgndsub, ihgndtopad, ihmaintopad,\
      ivgndtopad, ivmaintopad, ivgndopening, ihgndopening,\
      ihedge, ivedge, ivgndfiller, ihgndfiller])
    # Placement of the bondpads
    igndbp = all_cells['gndfeed_bondpad_r']
    imainbp = all_cells['MSfeed_bondpad_r']

    igndbp_ref = gdspy.CellReference(igndbp)
    imainbp_ref = gdspy.CellReference(imainbp)

    gb_dx, gb_dy = get_size(igndbp_ref)
    moveleft(iresogndsub_ref, igndbp_ref, spacing=default_spacing)
    moveabove(iresogndsub_ref, igndbp_ref, spacing=gb_dx)
    moveleft(iresogndsub_ref, imainbp_ref, spacing=default_spacing)
    movebelow(igndbp_ref, imainbp_ref, spacing=default_spacing)

    i2_dx, i2dy = get_size(i2feed_ref)
    moveleft(igndbp_ref, i2feed_ref, spacing=default_spacing)
    moveabove(igndbp_ref, i2feed_ref, spacing=i2dy)
    movebelow(i2feed_ref, i2gnd_ref, spacing=default_spacing)
    centerx(i2feed_ref, i2gnd_ref)
    #moveleft(igndbp_ref, i2gnd_ref)

    mask.add(imainbp_ref)
    mask.add(igndbp_ref)
    mask.add(i2feed_ref)
    mask.add(i2gnd_ref)

    not_yet -= set([igndbp, imainbp, i2feed, i2gnd])

    # Placement of the corner cells
    # igndsubcorner = all_cells['gndsub_corner_r']
    # ifeedcorner = all_cells['feed_corner_r']

    # igndsubcorner_ref = gdspy.CellReference(igndsubcorner)
    # ifeedcorner_ref = gdspy.CellReference(ifeedcorner)

    # centerx(i2feed_ref, ifeedcorner_ref)
    # movebelow(i2gnd_ref, ifeedcorner_ref)
    # moveleft(ivia_ref, igndsubcorner_ref)
    # centery(ivia_ref, igndsubcorner_ref)

    # mask.add(igndsubcorner_ref)
    # mask.add(ifeedcorner_ref)

    # not_yet -= set([igndsubcorner, ifeedcorner])


    # Placement of the unique capacitors on the mask. I'll fit a set of 30 caps above the common pixel
    # another 30 below and 2 flanking each side of the common pixel
    caps_inv = [cell for cell in inv_cell_list \
            if cell.name.startswith('Capacitor') and cell.name.find('common')<0]
    caps_inv.sort(key=lambda x:x.name)

    num_above = 30
    num_below = 30
    num_rflank = 2
    num_lflank = 2

    mrows = 5
    mcols = 6

    num_start, num_end = 0, num_above
    caps_above = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_lflank
    caps_lflank = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_rflank
    caps_rflank = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_below
    caps_below = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:]))

    u_dx, u_dy = get_size(caps_above[0])

    u_dx += intercap_spacing
    u_dy += intercap_spacing

    ci_dx, ci_dy = get_size(icap2ind_ref)
    # Place the right flanking caps
    moveabove(icap2ind_ref, caps_rflank[0], spacing=-intercap_spacing)
    moveleft(icap2ind_ref, caps_rflank[0], spacing=-ci_dx/2)
    moveabove(icap2ind_ref, caps_rflank[1], spacing=-intercap_spacing)
    moveright(icap2ind_ref, caps_rflank[1], spacing=ci_dx/2)
    mask.add(caps_rflank[0])
    mask.add(caps_rflank[1])

    # Place the left flanking caps
    movebelow(icap2ind_ref, caps_lflank[0], spacing=intercap_spacing)
    moveleft(icap2ind_ref, caps_lflank[0], spacing=-ci_dx/2)
    movebelow(icap2ind_ref, caps_lflank[1], spacing=intercap_spacing)
    moveright(icap2ind_ref, caps_lflank[1], spacing=ci_dx/2)
    mask.add(caps_lflank[0])
    mask.add(caps_lflank[1])

    # Place all the caps above the common pixel
    for index, cap in enumerate(caps_above):
        irow, icol = index // mcols, index % mcols
        y_displacement = u_dy * (mrows - irow)
        x_displacement = u_dx * (icol - (mcols//2))
        moveabove(icap2ind_ref, cap, spacing=-y_displacement)
        cap.translate(x_displacement , intercap_spacing)
        mask.add(cap)

    # Place all the caps below the common pixel
    for index, cap in enumerate(caps_below):
        irow, icol = index // mcols, index % mcols
        y_displacement = u_dy * (mrows - irow)
        x_displacement = u_dx * (icol - (mcols//2))
        movebelow(icap2ind_ref, cap, spacing=y_displacement)
        cap.translate(x_displacement , - intercap_spacing)
        mask.add(cap)


    # iucaps_ref = placeuniquecaps(caps_inv, mask, 10, nrows, ncols)

    # map(lambda x:mask.add(x), caps_above)
    # map(lambda x:mask.add(x), caps_below)
    # map(lambda x:mask.add(x), caps_rflank)
    # map(lambda x:mask.add(x), caps_lflank)

    not_yet -= set(caps_inv)

    # Add the island and xef2 release structure
    ixef2 = all_cells['XeF2_release_r']
    i_island = all_cells['LSN_Island_280umlegs_r']
    isl_dx, isl_dy = get_size(i_island)
    ixef2_ref = gdspy.CellReference(ixef2)
    i_island_ref = gdspy.CellReference(i_island)
    moveleft(icommoncap_ref, i_island_ref, spacing=default_spacing)
    moveabove(icommoncap_ref, i_island_ref, spacing=isl_dy)
    moveleft(icommoncap_ref, ixef2_ref, spacing=default_spacing)
    movebelow(i_island_ref, ixef2_ref, spacing=default_spacing)
    mask.add(i_island_ref)
    mask.add(ixef2_ref)


    filler = ssd.fill_empty_space(mask, mask_width, mask_length)
    mask.add(filler)


    print ("\n\nMask Generation Completed.\n\n")
    return mask

if __name__ == "__main__":
    # File name to be read. Must be in GDS format.
    fn = 'one_pixel.gds'
    main_lib.read_gds(fn)
    #cell_dict = main_lib.cell_dict

    ## Give the name of the global overlay cell in the file just read.
    #globaloverlay = cell_dict['GLOBAL_Overlay_Export']
    #top = cell_dict['Wafer_Layout']
    #inverted_cells = []
    #for cell in top.get_dependencies():
    #    inverted_cells += invert_cell(cell)

    #inverted_cells = [icell for icell in inverted_cells if icell is not None]
    #generate_inverted_overlay(top, set())
    ##for cell in inverted_cells:
    ##    print (cell)
    #gdspy.write_gds(fn, unit=1e-6,precision=1e-9)


    # I'll assume that all the necessary inverted cells have been generated.
    # I now want to make the reticle cells
    mask = generate_mask()
    canon_lens = get_mask_lens()
    mask.add(canon_lens)

    gdspy.write_gds(fn, unit=1e-6,precision=1e-9)
