#! /usr/bin/env python3


#########################################################
# This code is meant to model the shift in the Qc's with#
# the addition of a G10 spacer that is about 1.5 mm thick
# The spacer should cause the parasitic capacitance to go
# down by a factor of about 8

import numpy as np
from math import pi
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from capacitors import IDC
from scipy.constants import epsilon_0
import pprint
pp = pprint.PrettyPrinter(indent=4)

nH = 1e-9
pF = 1e-12
MHz = 1e6
um = 1e-6
mm = 1e-3
er_si = 11.8
er_G10 = 5
d_si = 0.5*mm
d_G10 = 1/16*25.4*mm
L = 12*nH
Z0 = 50 # Ohms
d = d_si + d_G10
er_eff = er_si * er_G10 * d / (er_si * d_G10 + er_G10 * d_si)
#print (er_eff)

def get_default_colors():
    prop_cycle = plt.rcParams['axes.prop_cycle']
    return prop_cycle.by_key()['color']





w = g = 2
fl = 1000
fg = 1.5
l = np.array([1342, 1466, 1602, 1754, 1930])
nfingers = l/(w+g)
Cs = np.zeros(5)
for i in range(5):
    cap = IDC(1.0)
    cap.set_dimensions(w,g,fl, fg, nfingers[i])
    Cs[i] = cap.capacitance()/pF



fr_expected = 1/2/pi/(L*Cs*pF)**0.5/MHz
wr_expected = 2*pi*fr_expected*MHz
fr_meas = np.array([277.638, 294.601, 310.120, 323.85, 335.85])[::-1]

# print (Cs)
# print (fr_expected - fr_meas)

# First, I want to make a prediction for the expected Qc before the spacer
# is added. This will be a good way to check my understanding of the Physics 
# involved

l_cc = 162
cc_cap = IDC(1.0)
cc_cap.set_dimensions(w,g,75,fg,l_cc/(w+g))
Cc = cc_cap.capacitance()/pF
# Cc = 0.28036 # Let's try this

l_tank = 2002 #um. Full length of the resonator tank
A_contact = 25*l_tank*um**2
A_couplingcap = (l_cc/(w+g)*75*w)*um**2
A_short_edge = (0.5*l/(w+g)*fl*w)*um**2 + A_contact + A_couplingcap
A_long_edge = (0.5*l_tank/(w+g)*fl*w)*um**2 + A_contact

pp.pprint(A_short_edge/A_long_edge)
Cp1 = (epsilon_0 * er_si * A_short_edge/d_si)
Cp2 = (epsilon_0 * er_si * A_long_edge/d_si)
Zb = 1/(1j*wr_expected*Cp2) + 1/((1j*wr_expected*Cp1) + 1/(Z0/2 + 1/(1j*wr_expected*Cc*pF)))
E_stored = 0.5*Cs*pF
# Analytical expression for Pdiss: P_diss = V^2 Re(Z)/Abs(Z)^2
# P_diss = V^2 * (2 Cc^2 Cp2^2 Z0 wr^2)/(4(Cp1 + Cp2 + Cc)^2 + Cc^2(Cp1+Cp2)^2 Z0^2 wr^2)
# In the limit Cp1 -> 0, Cp2 -> Cc, then P_diss ~ Cc^2 Z0 wr^2 / 8
# Qc = wr * (1/2 * C * V^2)/P_diss
P_diss = 0.5*Zb.real/np.abs(Zb)**2
Qc_expected = wr_expected * E_stored/P_diss
Qc_measured = np.array([72925., 85417., 106397., 126574., 155298.])
for fr,qc, qcm in zip(fr_expected, Qc_expected, Qc_measured):
    print ("{0:3.2f} MHz : expected/measured {1:3.2f}".format(fr, qc/qcm))
# ---------------------------------------------------------------- #
# A wrong way to calculate Qc

# A_total = A_short_edge + A_long_edge
# # print (A_total/mm**2)
# Cp = (epsilon_0 * er_si * (A_total/2)/d_si)/pF
# # print(Cp)
# Ccp = (Cc * Cp)/(Cc + 2*Cp)
# Qc_expected = 8*Cs*pF/(Z0*wr_expected*(Ccp*pF)**2)

# ---------------------------------------------------------------- #



# Now we add a spacer between the chip and the chip holder. Here we will try to
# include a correction for the fringing fields
R_short_edge = (A_short_edge/pi)**0.5
R_long_edge = (A_long_edge/pi)**0.5
Cp1_sp = (epsilon_0 * er_eff * A_short_edge/d) +\
    epsilon_0*er_eff*R_short_edge*(np.log(16*pi*R_short_edge/d)-1)
Cp2_sp = (epsilon_0 * er_eff * A_long_edge/d) +\
    epsilon_0*er_eff*R_long_edge*(np.log(16*pi*R_long_edge/d)-1)
Zb_sp = 1/(1j*wr_expected*Cp2_sp) + 1/((1j*wr_expected*Cp1_sp) + 1/(Z0/2 + 1/(1j*wr_expected*Cc*pF)))
P_diss_sp = Zb_sp.real/np.abs(Zb_sp)**2
Qc_sp = wr_expected * E_stored / P_diss_sp
Qc_meas_sp = np.array([59087., 83590., 115168., 202963., 424276.])
print ("""
Recalculating to account for the addition of the G10 spacer. Parallel plate capacitor approx
    """)
for fr,qc, qcm in zip(fr_expected, Qc_sp, Qc_meas_sp):
    print ("{0:3.2f} MHz : expected/measured {1:3.2f}".format(fr, qc/qcm))
# , measured {2:3.2f}
# print (Qc_sp/Qc_expected)
# print (Qc_meas_sp/Qc_measured)

# Retry with an overestimate of the total parasitic capacitance. At the opposite limit, 
# we only see a point charge infront of an infinite grounded metal plane. Use this as 
# the limit
Cp1_ol = 8*epsilon_0*er_eff*pi*d*np.ones(5)
Cp2_ol = 8*epsilon_0*er_eff*pi*d
Zb_ol = 1/(1j*wr_expected*Cp2_ol) + 1/((1j*wr_expected*Cp1_ol) + 1/(Z0/2 + 1/(1j*wr_expected*Cc*pF)))
P_diss_ol = Zb_ol.real/np.abs(Zb_ol)**2
Qc_ol = wr_expected * E_stored / P_diss_ol

print ("""
Recalculating to account for the addition of the G10 spacer. Point charge approx
    """)
for fr,qc, qcm in zip(fr_expected, Qc_ol, Qc_meas_sp):
    print ("{0:3.2f} MHz : expected/measured {1:3.2f}".format(fr, qc/qcm))

#pp.pprint (Cp1_ol/Cp1_sp)
#pp.pprint (Cp2_ol/Cp2_sp)


# Let's try and plot Qc as a function of Cp = (Cp1 + Cp2)/2 = (1 + x)Cp2/2. The
# exact formula is Qc = C (2 Cp + Cc)^2 (1 + x)^2/(4 Cc^2 Cp^2 Z0 wr) + 4 C Z0
# wr (1 + x)^2

# Pick a logarithmic range of Cps
minrange = int(np.log10(Cp2_sp))-1
maxrange = int(np.log10(Cp2_ol))
Cp = np.logspace(minrange, maxrange, 100)
x = Cp1_sp/Cp2_sp
wr = wr_expected
fig, ax = plt.subplots(figsize=(10,10))
# Get the list of default colors
colors = get_default_colors()
for i in range(len(x)):
    qc = (Cs[i]*pF*(2*Cp + Cc*pF)**2*(1 + x[i])**2)/(4*(Cc*pF)**2*Cp**2*Z0*wr[i])
    qc += 4*Cs[i]*pF*Z0*wr[i]*(1+x[i])**2
    cp2 = 2*Cp/(1 + x[i])
    ax.hlines(Qc_meas_sp[i], cp2[0]/pF, cp2[-1]/pF, linestyles='dashed',\
        colors=colors[i])
    ax.plot(cp2/pF, qc, label='x = {0:1.3f}'.format(x[i]))
    #print (qc[0], Qc_meas_sp[i])
    ax.set_yscale('log')
    ax.set_xscale('log')

ax.set_xlabel('Cp2 [pF]')
ax.set_ylabel('Qc')
ax.grid(which='both')
ax.legend(title='x = Cp1/Cp2', loc='best')
ax.set_ylim(1e4, 1e6)
plt.show()
