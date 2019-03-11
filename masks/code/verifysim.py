"""
Module for comparing the gds file that is the output of steppersim.py to the 
original layout. In the original gds file, there must be a layout cell made out
of the versions of the constituent cells that are present in the reticles. For
example, if the original layout cell is called "layout", a second cell called
"layout_inverted" can be generated using the cells from the reticles. This
second cell is the one that ought to be compared with the output from
steppersim.py. The default cell in the output gds file from steppersim.py is
called "wafer". Saves the output file in [inputfilename]_verifysim.gds

The features in the cells are flattened and then merged. The exclusive or of the
polygons in the original layout cell is made with the cells from the simulated
cell. This process may take a while depending on how many polygons are
generated. In the ideal case, there are no mismatches and the output is purely
empty.

Example:
	python3 verifysim.py --orig-file maskfile.gds --sim-file
		maskfile_sim-output.gds --orig-cell layout_inverted [--sim-cell wafer]
		[--output-file maskfile_verifysim.gds]

The output file has 3 cells:
	1. original_plus_mismatch: compares all the polygons from the original cell
		(in layer 1) and the mismatch (layer 2)
	2. simulated_plus_mismatch: compares all the polygons from the simulated cell
		(in layer 3) and the mismatch (layer 2)
	3. original_plus_sim_plus_mismatch: overlays all the polygons from the
		original cell (in layer 1), simulated cell (in layer 2) and the
		mismatch (layer 2)
	4. mismatch:
By default, the output maskfile is maskfile_sim-output.gds where maskfile.gds is
the name of the input file.

"""
import gdspy
import argparse
import os

orig = "../mask_files/sscale_darkres_4in_with_island_20190310.gds"
sim = "../mask_files/sscale_darkres_4in_with_island_20190310_sim-output.gds"

parser = argparse.ArgumentParser(description="""Stepper verify. Given an output
		file from steppersim.py and the original file, checks how well the
		steppersim output matches the original file.""")
parser.add_argument("--orig-file", metavar='orig', type=str,
	help='name of the original file in GDS format')
parser.add_argument("--sim-file", metavar='sim', type=str,\
	help="""name of the file output from steppersim.py.""")
parser.add_argument("--orig-cell", metavar='orig_cell', type=str,\
		help="""name of the inverted layout cell in orig-file.""")
parser.add_argument("--output-file", metavar='fn_out', type=str,\
		help="""name of the output file.""")
parser.add_argument("--sim-cell", metavar='sim_cell', nargs="?", type=str,\
		help="""name of the inverted layout cell in the file from steppersim.py.""")

orig_cell = 'Global_Overlay_4in_Inverted'
sim_cell = 'wafer'

def verify_layout(orig_file, sim_file, orig_cellname, sim_cellname, fn_out):
	original = gdspy.GdsLibrary(orig_file)
	simulated = gdspy.GdsLibrary(sim_file)

	original.read_gds(orig)
	simulated.read_gds(sim)

	global_orig = original.extract(orig_cellname).flatten().get_polygons()
	global_sim = simulated.extract(sim_cellname).flatten().get_polygons()
	print ("Original cell has {:d} polygons".format(len(global_orig)))

	print ("Simulated cell has {:d} polygons".format(len(global_sim)))
	#global_orig = gdspy.fast_boolean(global_orig, None, 'or')
	#global_sim = gdspy.fast_boolean(global_sim, None, 'or')
	print ("Calculating the overlap ...")
	overlap = gdspy.fast_boolean(global_orig, global_sim, 'xor', layer=1)
	print ("Overlap has {:d} polygons".format(len(overlap.polygons)))
	print ("Overlap calculation done. Writing to file")

	mismatch = gdspy.GdsLibrary('mismatch', unit=1e-6, precision=1e-9)
	miss = gdspy.Cell('miss')
	miss.add(overlap)

	org = gdspy.Cell('original_plus_mismatch')
	for poly in global_orig:
		org.add(gdspy.Polygon(poly, layer=2))
	for poly in overlap.polygons:
		org.add(gdspy.Polygon(poly, layer=1))

	simd = gdspy.Cell('simulated_plus_mismatch')
	for poly in global_sim:
		simd.add(gdspy.Polygon(poly, layer=3))
	for poly in overlap.polygons:
		simd.add(gdspy.Polygon(poly, layer=1))

	dc = gdspy.Cell('original_plus_sim_plus_mismatch')
	for poly in global_orig:
		dc.add(gdspy.Polygon(poly, layer=2))
	for poly in global_sim:
		dc.add(gdspy.Polygon(poly, layer=3))
	for poly in overlap.polygons:
		dc.add(gdspy.Polygon(poly, layer=1))
	mismatch.add(miss)
	mismatch.add(org)
	mismatch.add(dc)
	#pdb.set_trace()
	mismatch.write_gds(fn_out, [org, simd, miss, dc])


if __name__=='__main__':
	args = parser.parse_args()
	orig = args.orig_file
	sim = args.sim_file
	orig_cellname = args.orig_cell
	sim_cellname = args.sim_cell if args.sim_cell else "wafer"
	fn_out = args.output_file if args.output_file else os.path.splitext(orig)[0] + '_verifysim.gds'
	verify_layout(orig, sim, orig_cellname, sim_cellname, fn_out)
