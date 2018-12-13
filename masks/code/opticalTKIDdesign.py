#! /usr/bin/env python3

import numpy as np
from math import pi
from capacitors import IDC
import gdspy
import matplotlib.pyplot as plt
from scipy.special import iv, kn
from scipy.constants import h, k

K_0 = lambda x: kn(0, x)
I_0 = lambda x: iv(0, x)

nH = 1e-9
pF = 1e-12
MHz = 1e6
cm = 1e-2
g = 1e-3
mJ = 1e-3

makeplots = False
# Aluminum material properties
gamma = 1.35 * mJ#/mol/Kelvin**2
density = 2.7 * g/cm**3
A_r = 26.98 * g#/u.mol
N_0 = ((3 * (gamma * (density/A_r)))/(2*np.pi**2 * k**2))

Z0 = 50
L = 10*nH

# Rounds numbers to the nearest base
def roundto(num, base):
    return ((num // base) + ((num % base) > (base//2))) * base

def get_MB_Qi(T, f):
    Q_int = 3e7
    Tc = 1.4
    alpha = 0.5
    Delta = 1.764 * k * Tc
    eta = (h * f / (2 * k * T))
    S_1 = ((2/pi)*np.sqrt(2*Delta/(pi*k*T))* np.sinh(eta)*K_0(eta))
    n_th = (2*N_0 * np.sqrt(2*pi* k * T* Delta)*np.exp(-Delta/(k*T)))
    Q_qp = ((2 * N_0 * Delta)/(alpha * S_1 * n_th))
    Q_i = 1./(1/Q_qp + 1./Q_int)

    return Q_i

if __name__=="__main__":
    df = 5 #MHz
    N = 10
    fstart = 300 #MHz
    fr = fstart + np.arange(N) * df
    target_T = 380e-3
    Qi = get_MB_Qi(target_T, fr*MHz)
    #print (Qi)
    #print (Qi)
    wr = 2*pi*fr*MHz
    C = 1./(wr**2*L)
    C_branch = 0.11 * pF
    C_load = C_branch * N
    Zin = Z0
    Zout = 1/(1./Z0 + 1j*wr*C_load)
    Cc = np.sqrt((2*C)/(Qi*Z0*wr))
    CC = np.average(Cc) # use the same coupling cap for all the resonators
    #CC = 0.1945 * pF # Using a number from actual calculations
    y = wr * CC * Z0
    Gprime = (wr*CC*Zin*y/(Zin + Zout)) - (1j*wr*CC*Zin**2*y**2)/(Zin + Zout)**2
    dQe = Gprime/(wr*C)
    Qe = 1/dQe
    Qc = 1./np.real(dQe)
    #print (Qc)
    #print (Qc)
    phi_c = np.arctan2(dQe.imag, dQe.real)
    #print (phi_c*180/pi)
    #Qc = (C*pF)/(0.5*Z0*wr*(CC*pF/2)**2)
    Qr = 1./(1./Qc + 1./Qi)
    #chi_c = 4*Qr**2/Qi/Qc/np.cos(phi_c)
    #chi_c = 4*np.abs(Qe)*Qi/(Qi*np.cos(phi_c) + np.abs(Qe))**2
    chi_c = 4*Qc*Qi/(Qi + Qc)**2#/np.cos(phi_c)

    plt.figure()
    plt.scatter(np.arange(N), chi_c)
    plt.xlabel('Resonator Index')
    plt.ylabel('Coupling Efficiency chi_c')
    plt.grid()
    plt.savefig('opticalTKID_coupling_efficiency.png')
    plt.close()
    dQr = 1/Qr

    #exit()
    if makeplots:
        for i in range(N):
            frequency = np.linspace(-1.0, 1.0, 1000) + fr[i]
            x = (frequency - fr[i])/fr[i]
            S21 = 1 - dQe[i]/(dQr[i] + 1j * 2 * x)
            label = 'Qr=%d Qi=%d\nQc=%d phi_c=%1.3fdeg'%(Qr[i], Qi[i], Qc[i],
                    phi_c[i]*180/pi)
            S21dB = 10*np.log10(np.abs(S21))


            plt.figure(i)
            plt.plot(frequency, S21dB, 'b', label=label)
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('|S21|')
            plt.legend(loc='center right')
            plt.grid()
            plt.savefig("reso_%1.3fMHz_S21.png"%fr[i])
            plt.close()

            plt.figure(i + 250)
            plt.plot(S21.real, S21.imag, 'b', label=label)
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.legend(loc='center right')
            plt.grid()
            plt.axis('square')
            plt.savefig("reso_%1.3fMHz_IQ.png"%fr[i])
            plt.close()
    #exit()

    print ("Total Coupling Capacitance Needed ", CC/pF)
    CC *= 2 #Needed because we have 2 coupling capacitors in series
    # I want to pick the largest capacitor as the base IDC and then design blade
    # cells to create the other capacitors
    Cmain = np.max(C)
    print (C/pF)
    cap = IDC(1.0)
    trace_width = 2.
    gap_width = 2.
    finger_length = 1000.
    finger_gap = 2.
    contact_width = 2
    nfinger = 10 # For now
    cap.set_dimensions(trace_width, gap_width,
            finger_length, finger_gap, nfinger, contact_width)

    # Want to determine the number of fingers for this capacitor structure.
    Nfingers, capfrac = cap.getnumfingers(Cmain)
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
    coup_nfingers, _ = coupcap.getnumfingers(CC)
    coup_nfingers += 2
    print ("Number of coupcap fingers ", coup_nfingers)
    coupcap.set_dimensions(trace_width, gap_width,
            100, finger_gap, coup_nfingers, contact_width)
    print ("Target coupling capacitance", CC/pF)
    print ("Actual coupling capacitance", coupcap.capacitance()/pF)
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
    coarse_blades = np.zeros(N)
    fine_blades = np.zeros(N)
    for i in range(N):
        cap_nfingers, cap_frac = cap.getnumfingers(C[i])
        coarse_blades[i]  = cap.nfinger - cap_nfingers
        fine_blades[i]  = roundto(cap_frac * cap.width, 1)

    print ("coarse blade cell sizes: ", coarse_blades)
    print ("fine blade cell sizes: ", fine_blades)
    blade_width = coarse_blades *(trace_width + gap_width) + gap_width/2

    #blade_cells = []
    #for i in range(N-1):
    #    bw = blade_width[i]
    #    rect = gdspy.Rectangle([-bw/2, bh/2], [bw/2, -bh/2], layer=1)
    #    cellname = "blade_%1.0fMHz"%fr[i+1]
    #    blade_cell = gdspy.Cell(cellname)
    #    blade_cell.add(rect)
    #    blade_cells.append(blade_cell)

    cdy = cap.width + 4 # Weird that not cap height.
    cdx = 20
    crect = gdspy.Rectangle([-cdx/2, cdy/2], [cdx/2, -cdy/2], layer=1)
    c_cellname = "coarse_blade"
    main_blade = gdspy.Cell(c_cellname)
    main_blade.add(crect)
    fdx, fdy = 100, 2
    frect = gdspy.Rectangle([-fdx/2, fdy/2], [fdx/2, -fdy/2], layer=1)
    f_cellname = "fine_blade"
    fine_blade = gdspy.Cell(f_cellname)
    fine_blade.add(frect)
    # Want to demonstrate the blading technique that we will use for the fab
    demo_blades = []
    for i in range(N):
        cellname = "Cap_%1.0fMHz_with_blades"%fr[i]
        bladed_cell = gdspy.Cell(cellname)
        cap_ref = gdspy.CellReference(cap_cell)
        bladed_cell.add(cap_ref)
        blade_ref = gdspy.CellReference(main_blade)
        (xmin, ymin), (xmax, ymax) = cap_ref.get_bounding_box()
        (cxmin, cymin), (cxmax, cymax) = blade_ref.get_bounding_box()
        cap_dx = (xmax - xmin)
        cap_dy = finger_length + finger_gap + 2*contact_width
        cx0 = coarse_blades[i] * (trace_width + gap_width) - gap_width/2 - cdx/2
        cy0 = 0
        blade_ref.translate(xmin+ cx0, cy0)
        bladed_cell.add(blade_ref)
        # To add the fine blading, we need to account for the last finger that
        # is bladed off by the main blade. The zeroth finger is always attached
        # to the lower side of the capacitor.
        sign = -1 if coarse_blades[i] % 2 else 1
        fine_blade_ref = gdspy.CellReference(fine_blade)
        fx0 = (coarse_blades[i] + 1) * (trace_width + gap_width) - gap_width/2
        fx0 -= fdx/2
        fy0 = sign * (cap_dy/2 - contact_width - finger_gap - finger_length +
                fine_blades[i] + fdy/2)
        #print (cap_dy, fine_blades[i], fdy)
        print (i)
        print (xmin + cx0, cy0)
        print (xmin + fx0, fy0)
        fine_blade_ref.translate(xmin+fx0, fy0)
        bladed_cell.add(fine_blade_ref)

        demo_blades.append(bladed_cell)

    gdspy.write_gds("opticalTKID_IDC.gds",
            cells=[cap_cell, main_blade, fine_blade] + demo_blades,unit=1e-6,
            precision=1e-9)




