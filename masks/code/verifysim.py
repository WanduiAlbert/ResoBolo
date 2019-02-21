import gdspy
import pdb

#orig = "diplexer_FINAL_reformating.gds"
#sim = "diplexer_FINAL_reformating_sim-output.gds"
orig = "../mask_files/sscale_darkres_4in_with_island.gds"
sim = "../mask_files/sscale_darkres_4in_with_island_sim-output.gds"

orig_cell = 'Global_Overlay_4in_Inverted'
sim_cell = 'wafer'

original = gdspy.GdsLibrary('orig')
simulated = gdspy.GdsLibrary('sim')

original.read_gds(orig)
simulated.read_gds(sim)

global_orig = original.extract(orig_cell).flatten().get_polygons()
global_sim = simulated.extract(sim_cell).flatten().get_polygons()
print ("Orig has {:d} polygons".format(len(global_orig)))

print ("Sim has {:d} polygons".format(len(global_sim)))
#global_orig = gdspy.fast_boolean(global_orig, None, 'or')
#global_sim = gdspy.fast_boolean(global_sim, None, 'or')
print ("Calculating the overlap ...")
overlap = gdspy.fast_boolean(global_orig, global_sim, 'xor', layer=1)
print ("Overlap has {:d} polygons".format(len(overlap.polygons)))
print ("Overlap calculation done. Writing to file")

mismatch = gdspy.GdsLibrary('mismatch', unit=1e-6, precision=1e-9)
miss = gdspy.Cell('miss')
miss.add(overlap)

org = gdspy.Cell('original')
for poly in global_orig:
    org.add(gdspy.Polygon(poly, layer=2))
for poly in overlap.polygons:
    org.add(gdspy.Polygon(poly, layer=1))

simd = gdspy.Cell('simulated')
for poly in global_sim:
    simd.add(gdspy.Polygon(poly, layer=3))
for poly in overlap.polygons:
    simd.add(gdspy.Polygon(poly, layer=1))

dc = gdspy.Cell('overlapped')
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
mismatch.write_gds("../mask_files/sscale_darkres_4in_with_island_verifysim.gds", [org, simd, miss, dc])
