#!/usr/bin/env python3

import gdspy
import patches
import numpy as np

lib = gdspy.GdsLibrary('main')
fn = '../resobolo_files/alex_mask.gds'
lib.read_gds(fn)
overlay = lib.cell_dict['double_res']

(xmin, ymin), (xmax, ymax) = overlay.get_bounding_box()

bbox = gdspy.Rectangle([xmin, ymin], [xmax, ymax], layer=0)

invert = gdspy.fast_boolean(overlay.elements, bbox,'xor', layer=1)

i_overlay = gdspy.Cell('double_res_inv')
i_overlay.add(invert)
lib.add(i_overlay)
lib.write_gds('../resobolo_files/alex_mask_winv.gds')
