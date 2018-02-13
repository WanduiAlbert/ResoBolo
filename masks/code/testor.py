#! /usr/bin/env python3

import darkresonators as dr
import numpy as np
import gdspy
import astropy.units as u


if __name__=="__main__":
  L = 22 * u.nH
  f = 500 * u.MHz

  Cs = (1/L/(2*np.pi*f)**2).to(u.F).value
  cap = dr.IDC(1.0)
  nfingers, capfrac = dr.getnumfingers(cap, Cs)
  cap.nfinger = nfingers + 3
  cap.capfrac = capfrac

  print (cap.nfinger)
  # Template capacitors for constructing all the capacitors
  model_cap_cells, model_cap_nfingers = dr.make_captank_models(60)

  cap_arrays = dr.make_capacitor(cap, model_cap_cells, model_cap_nfingers)

  testcap = gdspy.Cell('test')
  for arr in cap_arrays:
    testcap.add(arr)



  gdspy.LayoutViewer()
  #  origin= get_cap_position(cap_array, N_res - (i+1), N_res, nrows, ncols,\
  #      width_spacing, len_spacing, ind_start)
  #  cap_array.origin = origin
  #nfingers = dr.getnumfingers(cap, C)

