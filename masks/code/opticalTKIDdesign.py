#! /usr/bin/env python3

import numpy as np
from math import pi
from capacitors import IDC
import gdspy

nH = 1e-9
pF = 1e-12
MHz = 1e6


L = 10*nH
Z0 = 50

if __name__=="__main__":
    df = 2 #MHz
    N = 8
    Qi = 16000 #target Qi at 350mK
    fstart = 400 #MHz
    fr = fstart + np.arange(N) * df
    wr = 2*pi*fr*MHz
    C = 1./(wr**2*L)/pF # in pF
    Cc = np.sqrt((8*C*pF)/(Qi*Z0*wr))/pF
    CC = np.average(Cc) # use the same coupling cap for all the resonators
    Qc = (C*pF)/(0.5*Z0*wr*(CC*pF/2)**2)
    Qr = 1./(1./Qc + 1./Qi)
    chi_c = 4*Qr**2/Qi/Qc
    print ("Expected coupling efficiency is in range [%1.3f:%1.3f] "%(
        min(chi_c), max(chi_c)))

    # I want to pick the largest capacitor as the base IDC and then design blade
    # cells to create the other capacitors
    Cmain = np.max(C)
    cap = IDC(1.0)
    trace_width = 4.
    gap_width = 4.
    finger_length = 1000.
    finger_gap = 4.
    contact_width = 24
    nfinger = 10 # For now
    cap.set_dimensions(trace_width, gap_width,
            finger_length, finger_gap, nfinger, contact_width)

    # Want to determine the number of fingers for this capacitor structure.
    Nfingers, capfrac = cap.getnumfingers(Cmain*pF)
    print ("Expected number of fingers is ", Nfingers)
    print ("Expected number of fractional fingers is ", capfrac)

    cap.set_dimensions(trace_width, gap_width,
            finger_length, finger_gap, Nfingers, contact_width)
    C_actual = cap.capacitance()/pF
    print ("The expected capacitance of the structure is %1.3f pF"%C_actual)
    print ("Width/height: (%1.3f, %1.3f) um"%(cap.width, cap.height))

    cap.layer = 0
    cap.cellname = 'Cap'
    capcell = cap.draw()
    cap_cell = gdspy.Cell('Cap_%1.0fMHz'%fr[0])
    cref = gdspy.CellReference(capcell, rotation=90)
    cap_cell.add(cref)

    # Now we should add the coupling capacitors
    coupcap = IDC(1.0)
    coupcap.set_dimensions(trace_width, gap_width,
            100, finger_gap, 10, contact_width)
    coup_nfingers, _ = coupcap.getnumfingers(CC*pF)
    print ("Number of coupcap fingers ", coup_nfingers)
    coupcap.set_dimensions(trace_width, gap_width,
            100, finger_gap, coup_nfingers, contact_width)
    coupcap.layer = 0
    coupcap.cellname = 'CoupCap'
    coupcapcell = coupcap.draw()
    coupcap_ref = gdspy.CellReference(coupcapcell, rotation=90)
    y_translation = cap.width/2 + coupcap.width/2 - contact_width
    coupcap_ref.translate(0, y_translation)
    cap_cell.add(coupcap_ref)
    coupcap_ref2 = gdspy.CellReference(coupcapcell, rotation=90)
    coupcap_ref2.translate(0, -y_translation)
    cap_cell.add(coupcap_ref2)

    cap_cell.flatten()

    # Now I want to design the blading cells for each of the resonators
    all_blades = np.zeros(N)
    for i in range(N):
        cap_nfingers, _ = cap.getnumfingers(C[i]*pF)
        all_blades[i]  = cap.nfinger - cap_nfingers

    print (all_blades)
    all_blades = all_blades[1:]
    blade_width = all_blades *(trace_width + gap_width) + gap_width/2
    bh = cap.width # Weird that not cap height.

    blade_cells = []
    for i in range(N-1):
        bw = blade_width[i]
        rect = gdspy.Rectangle([-bw/2, bh/2], [bw/2, -bh/2], layer=0)
        cellname = "blade_%1.0fMHz"%fr[i+1]
        blade_cell = gdspy.Cell(cellname)
        blade_cell.add(rect)
        blade_cells.append(blade_cell)

    gdspy.write_gds("opticalTKID_IDC.gds",
            cells=[cap_cell] + blade_cells,unit=1e-6,
            precision=1e-9)




