#! /usr/bin/env python3

import numpy as np
import gdspy as gds
from collections import namedtuple

Point = namedtuple("Point", ['x', 'y'])
save_dir = "../resobolo_files/"


def generate_heater(pad_start, size):
    pad_end = Point(pad_start.x + size.x, pad_start.y - size.y)

    pad = gds.Rectangle(pad_start, pad_end) # heater pad on the island

    wire_w = 1 #width of the wires to the heater pad

    vl = 10 # short length of the path vertical segment
    hl = 412 # short length of the main heater wires
    spacing = 2

    spec = {"layer":1, "datatype":0}
    left_start = Point(pad_start.x + wire_w/2, pad_start.y)
    # Initialized the start of the heater wires
    left_path = gds.Path(wire_w, left_start)
    left_path.segment((0-pad_end.y + 3), direction='-y', **spec)
    left_path.segment(wire_w/2, direction = '-x', **spec)
    left_path.segment(pad_end.x/2 - spacing/2 - wire_w/2, direction = '+x', **spec)
    left_path.segment(wire_w/2, direction = '+y', **spec)
    left_path.segment(vl + spacing + wire_w, direction = '-y', **spec)
    left_path.segment(wire_w/2, direction = '-x', **spec)
    left_path.segment(hl + spacing + wire_w, direction = '+x', **spec)
    # left_path2.segment(wire_w, direction='+y', layer=1)

    right_start = Point(pad_end.x - wire_w/2, pad_start.y)
    right_path = gds.Path(wire_w, right_start)
    right_path.segment((0-pad_end.y + 3), direction='-y', **spec)
    right_path.segment(wire_w/2, direction = '+x', **spec)
    right_path.segment(pad_end.x/2 - spacing/2 - wire_w/2, direction = '-x', **spec)
    right_path.segment(wire_w/2, direction = '+y', **spec)
    right_path.segment(vl, direction = '-y', **spec)
    right_path.segment(wire_w/2, direction = '-x', **spec)
    right_path.segment(hl, direction = '+x', **spec)
    # left_path2.segment(wire_w, direction='+y', layer=1)

    heater = gds.Cell('heater')
    heater.add(pad)
    heater.add(left_path)
    heater.add(right_path)
    return heater

def generate_microstrip(start, length, s_width, f_width):
    ms_line = gds.Path(s_width, start)
    ms_line.segment(length, direction='+x', layer=1, datatype=0, final_width=f_width )
    microstrip = gds.Cell('microstrip')
    microstrip.add(ms_line)

    return microstrip

island = gds.GdsLibrary("island")

#make the heater pad and lines
heater_start = Point(0, 0)
heater_size = Point(17.5, 7)
heater = generate_heater(heater_start, heater_size)

#make the microstrip line from the filter to the island. We can now easily customize the length and widths of the line
ms_w = 2 # terminal width of the microstrip line
ms_w0 = 7 #maximum width of the microstrip line
ms_length = 5710
ms_start = Point(heater_start.x + 17.5 + 2.25 + 1.5, heater_start.y - 17 + 0.5 -6.5 -4.25) # relative position of the island leg from the heater origin.
microstrip = generate_microstrip(ms_start, ms_length, ms_w, ms_w0)

# We want to see the files in their relative positions so we add them to a top cell
top = gds.Cell('top')
top.add(heater)
top.add(microstrip)

# island.add(top)
island.add(heater)
island.add(microstrip)

savefile = open(save_dir + 'island_design.gds', 'wb')
island.write_gds(savefile)

gds.LayoutViewer()
