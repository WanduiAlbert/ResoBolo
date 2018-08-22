import gdspy


orig = "sscale_darkres.gds"
sim = "sscale_darkres_sim-output.gds"

orig_cell = 'Global_Overlay_Inverted'
sim_cell = 'wafer'

original = gdspy.GdsLibrary('orig')
simulated = gdspy.GdsLibrary('sim')

original.read_gds(orig)
simulated.read_gds(sim)

global_orig = original.extract(orig_cell).flatten().get_polygons()
global_sim = simulated.extract(sim_cell).flatten().get_polygons()

global_orig = gdspy.fast_boolean(global_orig, None, 'or')
global_sim = gdspy.fast_boolean(global_sim, None, 'or')
overlap = gdspy.fast_boolean(global_orig, global_sim, 'xor', layer=1)



mismatch = gdspy.GdsLibrary('mismatch', unit=1e-6, precision=1e-9)
miss = gdspy.Cell('miss')
miss.add(overlap)

org = gdspy.Cell('original')
org.add(global_orig)
org.add(overlap)

simd = gdspy.Cell('simulated')
simd.add(global_sim)
simd.add(overlap)

dc = gdspy.Cell('overlapped')
dc.add(global_orig)
global_sim.layers=[2]*len(global_sim.layers)
dc.add(global_sim)
dc.add(overlap)
mismatch.write_gds('sscale_darkres_mismatch.gds', [org, simd, miss, dc])
