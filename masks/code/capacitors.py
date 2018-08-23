import numpy as np
from scipy import special, constants
from math import pi
import pprint
import astropy.units as u
pp = pprint.PrettyPrinter(indent=5)

mUnits = 1e-6
pF = 1.0e12
class IDC():
    def __init__(self,capfrac):
        self.trace_width = 2.
        self.gap_width = 2.
        self.finger_length = 900.
        self.finger_gap = 2.
        self.nfinger = 64
        self.contact_width = self.trace_width
        self.width = self.contact_width + self.finger_gap + \
            self.finger_length + self.contact_width
        self.height = self.nfinger*(self.gap_width+self.trace_width) + \
            self.trace_width
        self.capfrac = capfrac
        self.C = self.capacitance()

    def set_dimensions(self, trace_width, gap_width, finger_length, finger_gap,\
            nfinger):
        self.trace_width = trace_width
        self.gap_width = gap_width
        self.finger_length = finger_length
        self.finger_gap = finger_gap
        self.nfinger = nfinger
        self.contact_width = self.trace_width
        self.width = self.contact_width + self.finger_gap + \
            self.finger_length + self.contact_width
        self.height = self.nfinger*(self.gap_width+self.trace_width) + \
            self.trace_width
        self.C = self.capacitance()

    def capacitance2(self):
        p = mUnits*(self.finger_length - self.finger_gap)
        q = mUnits*self.height
        s = mUnits*self.gap_width
        w = mUnits*self.trace_width
        a = s + w
        e = constants.epsilon_0 * 11.7

        c = 0.0
        for n in range(1,100):
            j = special.j0((2*n-1)*pi*s/(2*a))
            c += (1.0/(2*n-1)) * j*j

        c *= p*q*(4.*e/(pi*a))
        return c

    def capacitance(self):
        e0 = constants.epsilon_0
        e1 = 1.0
        es = 11.7

        w = mUnits*self.trace_width
        g = mUnits*self.gap_width
        h = mUnits*self.finger_length
        N = int(self.nfinger*self.capfrac) - 1
        l = 2. * (w + g)
        eta = 2.*w/l
        kiinf = np.sin(np.pi*eta / 2.)
        kiinfp = np.sqrt(1 - kiinf**2)
        keinf = 2.*np.sqrt(eta) / (1 + eta)
        keinfp = np.sqrt(1 - keinf**2)

        Kkiinf = special.ellipk(kiinf)
        Kkiinfp = special.ellipk(kiinfp)
        Kkeinf = special.ellipk(keinf)
        Kkeinfp = special.ellipk(keinfp)

        Ci = e0 * (1+es) * h * Kkiinf/Kkiinfp
        Ce = e0 * (1+es) * h * Kkeinf/Kkeinfp
        Ct = (N-3) * Ci/2 + 2*(Ci*Ce)/(Ci+Ce)
        return Ct

    def getnumfingers(cap, target_C):
        e0 = constants.epsilon_0
        e1 = 1.0
        es = 11.7

        w = mUnits*cap.trace_width
        g = mUnits*cap.gap_width
        h = mUnits*cap.finger_length
        l = 2. * (w + g)
        eta = 2.*w/l
        kiinf = np.sin(np.pi*eta / 2.)
        kiinfp = np.sqrt(1 - kiinf**2)
        keinf = 2.*np.sqrt(eta) / (1 + eta)
        keinfp = np.sqrt(1 - keinf**2)

        Kkiinf = special.ellipk(kiinf)
        Kkiinfp = special.ellipk(kiinfp)
        Kkeinf = special.ellipk(keinf)
        Kkeinfp = special.ellipk(keinfp)

        Ci = e0 * (1+es) * h * Kkiinf/Kkiinfp
        Ce = e0 * (1+es) * h * Kkeinf/Kkeinfp
        n = (target_C -  2*(Ci*Ce)/(Ci+Ce)) * 2/Ci + 4
        nfingers = int(n)
        capfrac = n - nfingers

        return nfingers, capfrac


if __name__=="__main__":
    tracew = 2
    gapw = 2
    finger_length = 75
    finger_gap = 1.5
    nfingers = [39.5*2, 42.5*2, 45.5*2, 48.5*2, 52.5*2]
    main_nfingers = [336.5,  365.5,  399.5,  437.5,  482.5]
    N = len(nfingers)
    Cs = []
    ccs = []
    Qcs = []
    Z0 = 30*u.Ohm
    L = 10*u.nH
    for i in range(N):
        c = IDC(1.0)
        c.contact_width=25
        c.set_dimensions(tracew, gapw, finger_length, finger_gap, nfingers[i]-1)
        ccs.append(c.capacitance()*pF)
        # Now for the main capacitor
        cmain = IDC(1.0)
        cmain.contact_width=25
        cmain.set_dimensions(tracew, gapw, 1000, finger_gap,\
            main_nfingers[i]-0.5)
        Cs.append(cmain.capacitance()*pF)

    Cs = np.array(Cs)*u.pF
    pp.pprint (Cs)
    ccs = np.array(ccs)*u.pF
    wrs = (1/(L*Cs)**0.5).to(1/u.s)
    frs = (wrs/2/pi).to(u.MHz)
    pp.pprint (frs)
    Qcs = (8*Cs/(wrs*ccs**2*Z0)).to(1)
    pp.pprint (Qcs)
    #pp.pprint (ccs)
    #pp.pprint (Cs)

    # Next we do the same for the larger spacing capacitors

    tracew = 4
    gapw = 4
    finger_length = 77
    finger_gap = 3
    nfingers = [60.5, 64.5, 68.5, 74.5, 80.5]
    main_nfingers = [168.5,  184.5,  200.5,  220.5,  242.5]
    N = len(nfingers)
    uCs = []
    uccs = []
    uQcs = []
    for i in range(N):
        c = IDC(1.0)
        c.contact_width=25
        c.set_dimensions(tracew, gapw, finger_length, finger_gap, nfingers[i]-0.5)
        # print (c.width, c.height)
        uccs.append(c.capacitance()*pF)
        # Now for the main capacitor
        cmain = IDC(1.0)
        cmain.contact_width=25
        cmain.set_dimensions(tracew, gapw, 997, finger_gap,\
            main_nfingers[i]-0.5)
        # print (cmain.width, cmain.height)
        uCs.append(cmain.capacitance()*pF)

    uCs = np.array(uCs)*u.pF
    pp.pprint (uCs)
    uccs = np.array(uccs)*u.pF
    uwrs = (1/(L*uCs)**0.5).to(1/u.s)
    ufrs = (uwrs/2/pi).to(u.MHz)
    pp.pprint (ufrs)
    uQcs = (8*uCs/(uwrs*uccs**2*Z0)).to(1)
    pp.pprint (uQcs)
    #pp.pprint (ccs)
    #pp.pprint (Cs)