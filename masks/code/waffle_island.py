#! /usr/bin/env python3

import numpy as np
import gdspy as gds
from collections import namedtuple

Point = namedtuple('Point', ['x', 'y'])
settings = {"layer":0, "datatype"}
#Island dimensions
l = 500
w = 150

start = Point(0, 0)
end = Point(l ,-w)
fsn = gds.Rectangle(start, end, layer=0, datatype=0)

island = gds.Cell('island')
island.add(fsn)

#Define the perforation for the island
perf_l, perf_w = 10,10
bound_l, bound_w = 20, 20
perf_start = Point(0, 0)
perf_end = Point(perf_l, -perf_w)
bound_start = Point(0,0)
bound_end = Point(bound_l, -bound_w)

perforation = gds.Rectangle(perf_start, perf_end, layer=1, datatype=0)
boundary = gds.Rectangle(bound_start, bound_end, layer=2, datatype=0)
perforation.translate(perf_l/2, -perf_w/2)

waffle = gds.Cell('waffle')
waffle.add(boundary)
waffle.add(perforation)

spacing = 30
ncols = l/(bound_l + spacing//2)
nrows = w/(bound_l + spacing//2)
origin = Point(3.5, 0)

waffle_array = gds.CellArray(waffle, ncols, nrows, spacing, origin=origin)
island.add(waffle_array)

gds.LayoutViewer()
