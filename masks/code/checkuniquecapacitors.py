#!/usr/bin/env python3

import numpy as np
import gdspy
import sys
import matplotlib.pyplot as plt
import patches
import pdb
import pprint
pp = pprint.PrettyPrinter(indent=4)

# Set the library to write all the cells in
main_lib = gdspy.GdsLibrary('main', unit=1e-6, precision=1e-9)
gdspy.current_library = main_lib

def main():
	fn = '../mask_files/sscale_darkres_4in_with_island_20190510_AW.gds'
	main_lib.read_gds(fn)
	cell_dict = main_lib.cell_dict

	N = 64
	df = 3
	freqs = 300 + 3*np.arange(N)

	cellnames = list(map(lambda x: 'Capacitor_' + str(x) + 'MHz_r', freqs))

	for name1 in cellnames:
		for name2 in cellnames:
			if name1 == name2: continue
			poly1 = cell_dict[name1].get_polygons()
			poly2 = cell_dict[name2].get_polygons()
			overlap = gdspy.fast_boolean(poly1, poly2, 'xor')
			if overlap is None:
				print ('ERROR: cells ' + name1 + ' and ' + name2 + 'are identical')
				exit()
			else:
				print ("No overlap between " + name1 + " and " + name2)

if __name__=='__main__':
	main()


