#!/usr/bin/env python3

import numpy as np
import gdspy
import sys
import matplotlib.pyplot as plt
from scipy import constants,special
import astropy.units as u
from patchshot_composer import Shot, makeshot
from openpyxl import Workbook
from string import ascii_uppercase
number_of_points = 32

mUnits = 1e-6
wafer_width = 101e3
wafer_len = 101e3
inv_margin = 200
edge_margin = 1000
finger_length = 900
indboxwidth = 325.
indboxheight = 127.
coup_cap_finger_length = 100
bondpad_size = 140
island_halfwidth = 800
Qi = 40000
Z0 = 50 * u.Ohm
gnd_box_margin = 200
main_lib = gdspy.GdsLibrary('main')
gdspy.current_library = main_lib
wb = Workbook()
wbsavefile = "ResonatorArray.xlsx"
def_layers = {"Thin Gold":1, "PRO1":2, "ALUMINUM":3, "LSNSUB":4, "LSN1":5,
        "120nm_NbWiring":6, "STEPPER":7, "400nm_NbWiring":8, "ILD":9, "XeF2":10, "GP":12}
allshots = []
feed_cell_length = 100
feed_cell_width = 20
feed_main_width = 8
feed_ms_fraction = 0.2
feedres_dist = 209
mask_width = 22000
mask_length = 26000

# Need this functionality to work with excel spreadsheets
chars = list(ascii_uppercase)
char2num = {char:(ind+1) for ind, char in enumerate(chars)}
num2char = {(ind+1):char for ind, char in enumerate(chars)}


def moveabove_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_top = fixed_ury
    free_bottom = free_lly

    yshift = (fixed_top - free_bottom)
    return yshift

def moveabove(object_fixed,object_free,spacing=0):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    yshift = moveabove_get_shift(box_fixed,box_free)
    object_free.translate(0,yshift-spacing)

def movebelow_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_bottom = fixed_lly
    free_top = free_ury

    yshift = (fixed_bottom - free_top)
    return yshift

def movebelow(object_fixed,object_free,spacing=0):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    yshift = movebelow_get_shift(box_fixed,box_free)
    object_free.translate(0,yshift-spacing)

def moveleft_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_left = fixed_llx
    free_right = free_urx

    xshift = (fixed_left - free_right)
    return xshift

def moveleft(object_fixed,object_free,spacing=0):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    xshift = moveleft_get_shift(box_fixed,box_free)
    object_free.translate(xshift-spacing,0)

def moveright_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_right = fixed_urx
    free_right = free_llx

    xshift = (fixed_right - free_right)
    return xshift

def moveright(object_fixed,object_free,spacing=0):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    xshift = moveright_get_shift(box_fixed,box_free)
    object_free.translate(xshift-spacing,0)

def centerx_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_mid = 0.5*(fixed_llx + fixed_urx)
    free_mid = 0.5*(free_llx + free_urx)

    xshift = (fixed_mid - free_mid)
    return xshift

def centerx(object_fixed,object_free):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    xshift = centerx_get_shift(box_fixed,box_free)
    object_free.translate(xshift,0)

def centery_get_shift(box_fixed,box_free):
    (fixed_llx,fixed_lly),(fixed_urx,fixed_ury) = box_fixed
    (free_llx,free_lly),(free_urx,free_ury) = box_free

    fixed_mid = 0.5*(fixed_lly + fixed_ury)
    free_mid = 0.5*(free_lly + free_ury)

    yshift = (fixed_mid - free_mid)
    return yshift

def centery(object_fixed,object_free):
    box_fixed = object_fixed.get_bounding_box()
    box_free = object_free.get_bounding_box()
    yshift = centery_get_shift(box_fixed,box_free)
    object_free.translate(0,yshift)

def recenter(obj):
    (xmin, ymin), (xmax, ymax) = obj.get_bounding_box()
    midx = (xmax - xmin)//2
    midy = (ymax - ymin)//2
    x0, y0 = obj.origin
    dx = (x0 - midx)
    dy = (y0 - midy)

    obj.translate(dx, dy)

class InductorMeander():
    def __init__(self):
        self.tracewidth = 1.0
        self.gap = 1.0
        self.boxwidth = indboxwidth
        self.boxheight = indboxheight
        self.x0 = 0
        self.y0 = 0
        self.layer = 0
        self.cellname = 'Inductor'
    def draw(self):
        path = gdspy.Path(self.tracewidth,(self.x0,self.y0))
        specseg = {'layer':self.layer}
        specturn = {'layer':self.layer,'number_of_points':number_of_points}
        radius = 0.5*(self.gap + self.tracewidth)
        turnwidth = 0.5*self.gap + self.tracewidth
        path.segment(turnwidth,**specseg)
        seglen = self.boxwidth - 2*turnwidth
        while path.y - self.y0 + 4*(self.gap + self.tracewidth) < self.boxheight:
            path.segment(seglen,**specseg)
            path.turn(radius,'ll',**specturn)
            path.segment(seglen,**specseg)
            path.turn(radius,'rr',**specturn)
        path.segment(seglen,**specseg)
        path.turn(radius,'ll',**specturn)
        path.segment(seglen,**specseg)
        path.segment(turnwidth,**specseg)
        self.x1,self.y1 = path.x,path.y

        self.cell = gdspy.Cell(self.cellname)
        self.cell.add(path)
        return self.cell

class IDC():
    def __init__(self,capfrac):
        self.trace_width = 2.
        self.gap_width = 2.
        self.finger_length = 900.
        self.finger_gap = 2.
        self.nfinger = 64
        self.contact_width = self.trace_width
        self.width = self.contact_width + self.finger_gap + self.finger_length + self.contact_width
        self.height = self.nfinger*(self.gap_width+self.trace_width) + self.trace_width
        self.capfrac = capfrac
        self.cellname = 'Capacitor'
        self.layer= 6 
        self.C = self.capacitance()

    def draw(self, less_one=False):
        dx = (self.finger_length + self.gap_width )/2
        dy = self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
        self.make_fingers(less_one)
        self.make_contacts()
        self.C = self.capacitance()

        self.left_contact.layer = self.layer
        self.left_contact.translate(-dx, -dy/2)
        self.right_contact.translate(-dx, -dy/2)
        self.right_contact.layer = self.layer
        self.cell = gdspy.Cell(self.cellname)
        for f in self.fingers:
            f.layer = self.layer
            f.translate(-dx, -dy/2)
            self.cell.add(f)
        self.cell.add(self.left_contact)
        self.cell.add(self.right_contact)

        return self.cell

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

    def make_fingers(self, less_one):
        self.fingers = []
        xcur = 0.
        ycur = 0.

        left_fingers = []
        right_fingers = []

        nrf = (self.nfinger + 1) / 2
        nrfx = nrf * self.capfrac
        #print "nrfx: ",nrfx
        nrff = int(nrfx)    # Number of right fingers to leave fully intact
        nrfr = nrfx - nrff    # How much to truncate the last right finger

        minx = xcur + self.finger_gap
        maxx = xcur + self.finger_length
        partialx = nrfr * minx + (1-nrfr) * maxx

        range_val = self.nfinger + 1
        if less_one:
            range_val = self.nfinger

        for i in range(range_val):
            if i % 2 == 0:
                lower_leftx = xcur
                upper_rightx = xcur + self.finger_length
            else:
                if i/2 == nrff:
                    lower_leftx = partialx
                    upper_rightx = xcur + self.finger_length + self.finger_gap
                else:
                    lower_leftx = xcur+self.finger_gap
                    upper_rightx = lower_leftx + self.finger_length

            assert lower_leftx >= 0.
            assert upper_rightx <= self.finger_length + self.finger_gap

            lower_lefty = ycur
            upper_righty = lower_lefty + self.trace_width

            ycur += self.trace_width + self.gap_width
            box = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])
            if i % 2 == 0:
                self.fingers.append(box)
            elif (i/2) <= nrff:
                self.fingers.append(box)

    def make_contacts(self):
        lower_leftx = -self.contact_width
        lower_lefty = 0.
        upper_rightx = lower_leftx + self.contact_width
        upper_righty = lower_lefty + self.nfinger*(self.gap_width + self.trace_width) + self.trace_width
        self.left_contact = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])
        lower_leftx = self.finger_gap + self.finger_length
        upper_rightx = lower_leftx + self.contact_width
        self.right_contact = gdspy.Rectangle([lower_leftx,lower_lefty],[upper_rightx,upper_righty])

def Rectangle(x0,y0,x1,y1):
    r = gdspy.Rectangle((x0,y0),(x1,y1))
    def get_bounding_box():
        pt0,pt1,pt2,pt3 = r.points
        return pt0,pt2
    r.get_bounding_box = get_bounding_box
    return r

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

def makechipoutline(width, length, outlinename):
    specs = {'layer':10}
    outline = gdspy.Path(2, (0,0))
    #outline.layers = [0]
    outline.segment(width, '+x', layer=def_layers['PRO1'])
    outline.segment(length, '+y', layer=def_layers['PRO1'])
    outline.segment(width, '-x', layer=def_layers['PRO1'])
    outline.segment(length, '-y', layer=def_layers['PRO1'])
    outline.translate(-width/2, -length/2)
    outlinecell = gdspy.Cell(outlinename)
    outlinecell.add(outline)
    return outlinecell


def cellIsPresent(cellname):
    allcellnames = main_lib.cell_dict.keys()
    return cellname in allcellnames

def inverter(cell, rotation=0):
    cell_name = cell.name + '_r'
    if cellIsPresent(cell_name): return
    inv_cell = gdspy.Cell(cell_name)
    cell = cell.flatten()
    cell_ref = gdspy.CellReference(cell, rotation=rotation)
    dx, dy = get_size(cell_ref)
    dx += 2*inv_margin
    dy += 2*inv_margin
    dx = roundto(dx, 100)
    dy = roundto(dy, 100)
    layer = cell.get_layers().pop()
    polys = cell_ref.get_polygons(depth=1)
    polyset = gdspy.PolygonSet(polys, layer=layer)
    bbox = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=layer)
    new_polyset = gdspy.fast_boolean(polyset, bbox, 'xor', layer=layer)
    inv_cell.add(new_polyset)

def invert_cell(cell, rotation=0):
    layers = cell.get_layers()

    if len(layers) == 1:
        inverter(cell, rotation)

    for cell in cell.get_dependencies():
        invert_cell(cell, rotation)


# Note that I have defined island_halfwidth to be the distance in the x direction 
# between the interface between the capacitor and lines and between the inductor 
# and lines.
def make_cap_to_ind_lines():
    layer = def_layers["400nm_NbWiring"]
    cell_name = 'connector'
    ms_width = 2
    xlen_short = 5 * ms_width
    ind_overlap = 1.25 * ms_width
    xlen_long = island_halfwidth - xlen_short + ms_width + ind_overlap
    cap_overlap = ms_width * 2
    ylen = (finger_length + 3 * ms_width -ms_width)/2
    ylen_short = 6.25 * ms_width
    cap2ind_conn = gdspy.Cell(cell_name)
    #box = gdspy.Rectangle(layer=layer)
    conn1 = gdspy.Path(ms_width, (0,0))
    conn1.segment(xlen_short + cap_overlap, '+x', layer=layer)
    conn2 = gdspy.Path(ms_width, (conn1.x - ms_width/2, conn1.y + ms_width/2))
    conn2.segment(ylen , '-y', layer=layer)
    conn3 = gdspy.Path(ms_width, (conn2.x - ms_width/2, conn2.y + ms_width/2))
    conn3.segment(xlen_long, '+x', layer=layer)
    conn4 = gdspy.Path(ms_width, (conn3.x - ms_width/2, conn3.y - ms_width/2))
    conn4.segment(ylen_short, '+y', layer=layer)
    conn5 = gdspy.Path(2.5*ms_width, (conn4.x, conn4.y - ms_width/2))
    conn5.segment(51, '+y', layer=layer) # Just magic

    dx = island_halfwidth/2 + cap_overlap
    dy = ylen/2 - ms_width/2
    conn1.translate(-dx, dy)
    conn2.translate(-dx, dy)
    conn3.translate(-dx, dy)
    conn4.translate(-dx, dy)
    conn5.translate(-dx, dy)

    cap2ind_conn.add(conn1)
    cap2ind_conn.add(conn2)
    cap2ind_conn.add(conn3)
    cap2ind_conn.add(conn4)
    cap2ind_conn.add(conn5)

    return cap2ind_conn

def get_alignment_marks():
    fn = '../resobolo_files/alignment_marks.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    am = gdsii.extract('alignment_marks_patch_new')
    am_r = gdsii.extract('alignment_marks_patch_new_r')
    return am


def get_inductor():
    fn = '../resobolo_files/new_inductor.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    ind = gdsii.extract('Al_Inductor_Islands_right')
    polys = ind.get_polygons()
    (xmin, ymin), (xmax, ymax) = ind.get_bounding_box()
    dx = (xmax - xmin)
    offset = dx % 2
    #ind_view = gdspy.CellReference(ind, [-xmin, -ymin])
    inductor = gdspy.Cell('Al_inductor')
    for poly in polys:
        polygon = gdspy.Polygon(poly, layer=3)
        polygon = polygon.translate(-xmin, -ymin)
        polygon = polygon.translate(-(xmax - xmin)//2 + offset, -(ymax - ymin)/2)
        inductor.add(polygon)
    inductor.flatten()
    # print (inductor.get_bounding_box())
    return inductor

def get_size(cell):
    (xmin, ymin), (xmax, ymax) = cell.get_bounding_box()
    return (xmax - xmin), (ymax - ymin)

def get_cap_tank(cap, coupcap):
    Dx, Dy = get_size(cap)
    tofeed = gdspy.CellReference(coupcap)
    dx, dy = get_size(tofeed)
    #print (Dx, Dy, dx, dy)
    overlap = 2
    displacement_y = -(Dy - dy)/2
    displacement_x = (Dx + dx)/2 - overlap
    tofeed.origin = [0, 0]
    tofeed.translate(displacement_x, displacement_y)
    cap.add(tofeed)
    tognd = gdspy.CellReference(coupcap)
    tognd.translate(-displacement_x, displacement_y)
    cap.add(tognd)
    cap.flatten()
    return cap

def rot(n, i, j, ri, rj):
    if (rj == 0):
        if (ri == 1):
            i = (n-1) - i
            j = (n-1) - j
        i, j = j, i
    return (i, j)

def double_range(start, stop):
    while start < stop:
        yield start
        start <<= 1

# Maps from a position in an array d, into (i,j) coordinates of an n by n
# square. Assumes we are using a Hilbert curve to fill out the space
# d in {0, n**2 - 1}
def d2ij(d, n):
    i, j = 0, 0
    t = d
    for s in double_range(1,n):
        ri = 1 & (t//2)
        rj = 1 & (t ^ ri)
        i,j = rot(s, i, j, ri,rj)
        i += s*ri
        j += s*rj
        t //= 4
    return (i, j)

uint = np.frompyfunc(int, 1,1)

def get_feedline(xspacing):
    unit_length = 100
    n =  int (xspacing // unit_length )
    # print (n)
    s_gndsub = gdspy.Cell('unit_GP_sub')

    dx = n * feed_cell_length
    dy = feed_main_width
    mainline = gdspy.Rectangle([-dx/2, dy/2],\
            [dx/2, -dy/2], layer=def_layers['400nm_NbWiring'])
    main_cell = gdspy.Cell('feedline_main')
    main_cell.add(mainline)

    # dy = feed_cell_width
    # ild = gdspy.Rectangle([-dx/2, dy/2],\
    #     [dx/2, -dy/2], layer=def_layers['ILD'])
    # ild_cell = gdspy.Cell('feedline_ILD')
    # ild_cell.add(ild)

    # To make each cell more symmetric, I'll divide the ground sub region into 2
    sub_len = (1 - feed_ms_fraction)*feed_cell_length/2
    sub1 = gdspy.Rectangle([-sub_len/2, feed_cell_width/2],\
            [sub_len/2, -feed_cell_width/2], layer=def_layers['LSNSUB'])
    sub2 = gdspy.Rectangle([-sub_len/2, feed_cell_width/2],\
            [sub_len/2, -feed_cell_width/2], layer=def_layers['LSNSUB'])
    sdx = (sub_len + feed_ms_fraction * feed_cell_length)/2
    sub1.translate(-sdx, 0)
    sub2.translate(sdx, 0)
    s_gndsub.add(sub1)
    s_gndsub.add(sub2)
    gndsub_arr = gdspy.CellArray(s_gndsub, n, 1, [feed_cell_length,0])
    gndsub_arr.translate(-(dx - feed_cell_length)/2, 0)
    gndsub_cell = gdspy.Cell('feedline_GP_sub')
    gndsub_cell.add(gndsub_arr)
    gndsub_cell.flatten()

    feedcell = gdspy.Cell('MainFeedline')
    feedcell.add(gdspy.CellReference(main_cell))
    # feedcell.add(gdspy.CellReference(ild_cell))
    feedcell.add(gdspy.CellReference(gndsub_cell))
    # feedcell.add(gdspy.CellReference(main_cell))
    # feedcell.add(gdspy.CellReference(ild_cell))
    # feedcell.add(gdspy.CellReference(gndsub_cell))


    # feedline = gdspy.Cell('large_feedline_sec')
    # feedline.add(feedarr)
    # # feedline.flatten
    # print (get_size(feedarr))

    return feedcell

# Rounds numbers to the nearest base
def roundto(num, base):
    return ((num // base) + ((num % base) > (base//2))) * base

def get_resonator_ILDsub(common_cap, unique_cap, u_xoffset, gndsub_margin):
    c_dx, c_dy = get_size(common_cap)
    u_dx, u_dy = get_size(unique_cap)
    overlap = 2

    xmargin, ymargin = gndsub_margin
    xmargin -= 50
    ymargin -= 50

    dx = c_dx + u_dx - overlap + 2*xmargin
    dy = u_dy + 2 * ymargin

    # I want the ILD sub region to be a rounded to nearest 100
    dx = roundto(dx, 50)
    dy = roundto(dy, 50)

    # The margins have changed. Adjust them accordingly
    xmargin = (dx + overlap - c_dx - u_dx)//2
    excess = (dx + overlap - c_dx - u_dx) % 2
    ymargin = (dy - u_dy)//2

    ild = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=def_layers['ILD'])
    #gnd.translate(-offset, 0)

    i_xoffset = u_xoffset - (dx/2 - xmargin - u_dx/2)

    ildsub = gdspy.Cell('reso_ILD_sub')
    ildsub.add(ild)

    # print (c_dx, u_dx, get_size(ildsub))

    return ildsub, i_xoffset

# cres - common part of resonator
# ures - unique part of resonator
def get_resonator_GPsub(cres_bbox, ures_size):
    u_dx, u_dy = ures_size
    (xmin, ymin), (xmax, ymax) = cres_bbox
    c_dx, c_dy = (xmax - xmin, ymax - ymin)
    ymargin =  200
    xmargin = 300
    overlap = 2
    dx = c_dx + u_dx + 2 * xmargin - overlap
    dy = u_dy + 2 * ymargin

    # I want the gnd sub region to be a rounded to nearest 100
    dx = roundto(dx, 100)
    dy = roundto(dy, 100)
    # The margins have changed. Adjust them accordingly
    xmargin = (dx + overlap - c_dx - u_dx)//2
    excess = (dx + overlap - c_dx - u_dx) % 2
    ymargin = (dy - u_dy)//2
    #offset = dx/2 - (xmax + xmargin + overlap/2)
    gnd = gdspy.Rectangle([-dx/2, dy/2], [dx/2, -dy/2], layer=def_layers['LSNSUB'])
    #gnd.translate(-offset, 0)

    c_xoffset = dx/2 - (xmargin + excess + c_dx//2 + c_dx % 2)
    delta = c_dx//2 - (overlap + c_xoffset)
    u_xoffset = delta + u_dx//2 + u_dx % 1

    # Need to make a small GND plane sub region for connecting the coupling caps to the feedline
    sdx = 20
    sdy = feedres_dist
    small = gdspy.Rectangle([-sdx/2, sdy/2], [sdx/2, -sdy/2], layer=def_layers['LSNSUB'])
    small_xoffset = delta + sdx/2 - overlap + sdx/2 + feed_cell_length/2 -50
    small.translate(-small_xoffset, (dy + sdy)/2) #position relative to common arr

    gndsub = gdspy.Cell('reso_GP_sub')
    gndsub.add(gnd)
    gndsub.add(small)
    return gndsub, c_xoffset, u_xoffset, small_xoffset, [xmargin, ymargin]

def get_vertical_feed():
    allcells = main_lib.cell_dict
    gndsubcell = allcells['feedline_GP_sub']
    # ildcell = allcells['feedline_ILD']
    maincell = allcells['feedline_main']

    gndsubref = gdspy.CellReference(gndsubcell, rotation=90)
    # ildref = gdspy.CellReference(ildcell, rotation=90)
    mainref = gdspy.CellReference(maincell, rotation=90)

  

    vgndsubcell = gdspy.Cell('vert_feedline_GP_sub')
    vgndsubcell.add(gndsubref)
    vmaincell = gdspy.Cell('vert_feedline_main')
    vmaincell.add(mainref)
    # vmaincell.add(mainref2)

    vgndsub_ref = gdspy.CellReference(vgndsubcell)
    vmain_ref = gdspy.CellReference(vmaincell)

    l = feed_cell_width
    s = feed_main_width
    # Small adjustments to make the proper connections between
    # horizontal and vertical feedlines
    feed = gdspy.Rectangle([-s/2,s/2], [s/2, -s/2],\
            layer=def_layers['400nm_NbWiring'])
    feed_corner = gdspy.Cell('feed_corner')
    feed_corner.add(feed)
    feed_box = gdspy.CellReference(feed_corner)
    moveabove(vmain_ref, feed_box, spacing=s/2)
    # ild = gdspy.Rectangle([-l/2, l/2], [l/2, -l/2],\
    #     layer=def_layers['ILD'])
    # ild_corner = gdspy.Cell('ILD_corner')
    # ild_corner.add(ild)
    # ild_box = gdspy.CellReference(ild_corner)
    # moveabove(vild_ref, ild_box, spacing=l/2)
    gndsub = gdspy.Rectangle([-l/2, l/2], [l/2, -l/2],\
            layer=def_layers['LSNSUB'])
    gndsub_corner = gdspy.Cell('gndsub_corner')
    gndsub_corner.add(gndsub)
    gndsub_box = gdspy.CellReference(gndsub_corner)
    moveabove(vgndsub_ref, gndsub_box, spacing=l/2)
    # Same corrections at the bottom side
    feed_box_l = gdspy.CellReference(feed_corner)
    movebelow(vmain_ref, feed_box_l, spacing=-s/2)
    # ild_box_l = gdspy.CellReference(ild_corner)
    # movebelow(vild_ref, ild_box_l, spacing=-l/2)
    gndsub_box_l = gdspy.CellReference(gndsub_corner)
    movebelow(vgndsub_ref, gndsub_box_l, spacing=-l/2)
    main_w_corners = gdspy.Cell('vert_main_with_corners')
    gpsub_w_corners = gdspy.Cell('vert_gndsub_with_corners')

    main_w_corners.add(vmain_ref)
    main_w_corners.add(feed_box)
    main_w_corners.add(feed_box_l)
    main_w_corners.flatten()

    gpsub_w_corners.add(vgndsub_ref)
    gpsub_w_corners.add(gndsub_box)
    gpsub_w_corners.add(gndsub_box_l)
    gpsub_w_corners.flatten()

    vertfeedcell = gdspy.Cell('MainFeedline_vert')
    vertfeedcell.add(gdspy.CellReference(main_w_corners))
    # vertfeedcell.add(vild_ref)
    vertfeedcell.add(gdspy.CellReference(gpsub_w_corners))
    # vertfeedcell.add(feed_box)
    # # vertfeedcell.add(ild_box)
    # vertfeedcell.add(gndsub_box)
    # vertfeedcell.add(feed_box_l)
    # # vertfeedcell.add(ild_box_l)
    # vertfeedcell.add(gndsub_box_l)

    return vertfeedcell

def get_coupling_conns(g_dy, u_dy):
    ymargin = (g_dy - u_dy)//2
    w = 8
    l2gnd = ymargin + 2
    l2feed = l2gnd + feedres_dist
    l2gnd += 21 # A little extra

    c2gnd = gdspy.Rectangle([-w/2,l2gnd/2], [w/2, -l2gnd/2],\
            layer=def_layers['400nm_NbWiring'])
    c2feed = gdspy.Rectangle([-w/2,l2feed/2], [w/2, -l2feed/2],\
            layer=def_layers['400nm_NbWiring'])
    cap2feed = gdspy.Cell('cap_to_feed')
    cap2feed.add(c2feed)

    cap2gnd = gdspy.Cell('cap_to_gnd')
    cap2gnd.add(c2gnd)

    return cap2feed, cap2gnd

def get_common_resonator(ind, common_cap):
    common_resonator = gdspy.Cell('common_resonator')

    ms_width = 2
    connector = make_cap_to_ind_lines()
    conn_dx, conn_dy = get_size(connector)
    dy = conn_dy/2  + ms_width/2
    top_conn_ref = gdspy.CellReference(connector, [0, dy] )
    bot_conn_ref = gdspy.CellReference(connector, [0, -dy], x_reflection=True)
    resonator_connector = gdspy.Cell('Cap_to_Ind_lines')
    resonator_connector.add(top_conn_ref)
    resonator_connector.add(bot_conn_ref)
    resonator_connector.flatten()
    resonator_connector_ref = gdspy.CellReference(resonator_connector)
    ind_ref = gdspy.CellReference(ind)
    cap_ref = gdspy.CellReference(common_cap, rotation = 90)

    conn_dx, conn_dy = get_size(resonator_connector)
    ind_dx, ind_dy = get_size(ind_ref)
    cap_dx, cap_dy = get_size(cap_ref)

    capconn_overlap = 2 * ms_width
    indconn_overlap = 2 * ms_width
    # Get the total length of the common resonator
    dx = conn_dx + ind_dx + cap_dx - indconn_overlap - capconn_overlap
    # print (dx)
    # Currently the center is at the connector.
    offset = dx//2  - (cap_dx + (conn_dx - indconn_overlap - capconn_overlap)/2 )
    resonator_connector_ref.translate(-offset,0)
    common_resonator.add(resonator_connector_ref)
    # print (resonator_connector_ref.get_bounding_box())
    #ind_origin = [island_halfwidth/2  + indboxwidth/2, 0]
    #ind_ref = gdspy.CellReference(ind, origin = ind_origin)
    moveright(resonator_connector_ref, ind_ref, spacing=indconn_overlap)
    common_resonator.add(ind_ref)
    # print (ind_ref.get_bounding_box())
    moveleft(resonator_connector_ref, cap_ref, spacing=-capconn_overlap)
    # print (cap_ref.origin)
    common_resonator.add(cap_ref)
    # print (common_resonator.get_bounding_box())
    # print (ind_ref.get_bounding_box())
    # print (ind_ref.origin)
    return common_resonator

def get_bondpads():
    size = 400
    feedpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['400nm_NbWiring'])
    feedpad_cell = gdspy.Cell('MSfeed_bondpad')
    feedpad_cell.add(feedpad)

    # ildpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['ILD'])
    # ildpad_cell = gdspy.Cell('ILDfeed_bondpad')
    # ildpad_cell.add(ildpad)

    gndpad = gdspy.Rectangle([-size/2, size/2], [size/2, -size/2], layer=def_layers['120nm_NbWiring'])
    gndpad_cell = gdspy.Cell('GNDfeed_bondpad')
    gndpad_cell.add(gndpad)

    pad_cell = gdspy.Cell('Feed_Bonding_Pad')
    pad_cell.add(gdspy.CellReference(feedpad_cell))
    # pad_cell.add(gdspy.CellReference(ildpad_cell))
    pad_cell.add(gdspy.CellReference(gndpad_cell))
    return pad_cell

def make_inverted_cells():
    allcells = main_lib.cell_dict
    cellnames = allcells.keys()
    top = allcells['Global_Overlay']
    for cell in top.get_dependencies():
        if cell.name == 'WaferOutline': continue
        inv_cellname = cell.name + '_r'
        if inv_cellname not in cellnames:
            if inv_cellname.startswith('Capacitor'):
                invert_cell(cell, rotation=90)
            else:
                invert_cell(cell)



def get_inverted_cells():
    cells = main_lib.cell_dict
    cellnames = cells.keys()

    return [cells[name] for name in cellnames if name.endswith('_r') and name !=\
            'WaferOutline_r']

def get_cell_area(cell):
    dx, dy = get_size(cell)
    return dx * dy

def getunitfeedline():
    cells = main_lib.cell_dict
    dx = feed_cell_length
    dy = feed_main_width
    mainline = gdspy.Rectangle([-dx/2, dy/2],\
            [dx/2, -dy/2], layer=def_layers['400nm_NbWiring'])
    main_cell = gdspy.Cell('unit_main_feed')
    main_cell.add(mainline)


    gndsub_cell = cells['unit_GP_sub']
    # unitcell = gdspy.Cell('unit_feedline')
    # unitcell.add([gdspy.CellReference(gndsub_cell),\
    #     gdspy.CellReference(main_cell)])
    return main_cell, gndsub_cell

def getlinetopad():
    cells = main_lib.cell_dict
    umain, ugndsub = getunitfeedline()
    dx, dy = get_size(umain)
    n = int(roundto(6.5*edge_margin/dx, 10))
    main_hor_section = gdspy.CellArray(umain, n, 1, [dx, dy])
    gndsub_hor_section = gdspy.CellArray(ugndsub, n, 1, [dx, dy])

    m_dx, m_dy = get_size(main_hor_section)
    main_hor_section.translate(-(m_dx - dx)/2, 0)
    gndsub_hor_section.translate(-(m_dx - dx)/2, 0)

    # Make a cell for the horizontal part of the feedline too ground
    main_hor_cell = gdspy.Cell('main_hor_feedline_to_pad')
    gndsub_hor_cell = gdspy.Cell('gndsub_hor_feedline_to_pad')
    main_hor_cell.add(main_hor_section)
    main_hor_cell.flatten()
    gndsub_hor_cell.add(gndsub_hor_section)
    gndsub_hor_cell.flatten()

    main_vert_cell = gdspy.Cell('main_vert_feedline_to_pad')
    gndsub_vert_cell = gdspy.Cell('gndsub_vert_feedline_to_pad')
    mainfeed = cells['feedline_main']
    gndsub = cells['feedline_GP_sub']
    vdx, vdy = get_size(mainfeed)
    main_vert_section = gdspy.CellArray(mainfeed, 2, 1, [vdx , vdy ], rotation=90)
    gndsub_vert_section = gdspy.CellArray(gndsub, 2, 1, [vdx , vdy ], rotation=90)
    feed_corner = cells['feed_corner']
    gndsub_corner = cells['gndsub_corner']
    m_dx, m_dy = get_size(main_vert_section)
    u_dx, u_dy = get_size(mainfeed)
    main_vert_section.translate(0, -(m_dy - u_dx)/2)
    gndsub_vert_section.translate(0, -(m_dy - u_dx)/2)

    feed_corner_ref_top = gdspy.CellReference(feed_corner)
    feed_corner_ref_bot = gdspy.CellReference(feed_corner)
    feedspacing = (feed_cell_width + feed_main_width)/2
    moveabove(main_vert_section, feed_corner_ref_top)
    movebelow(main_vert_section, feed_corner_ref_bot)
    main_vert_cell.add(main_vert_section)
    main_vert_cell.add(feed_corner_ref_top)
    main_vert_cell.add(feed_corner_ref_bot)
    main_vert_cell.flatten()

    gndsub_corner_ref_top = gdspy.CellReference(gndsub_corner)
    gndsub_corner_ref_bot = gdspy.CellReference(gndsub_corner)
    moveabove(gndsub_vert_section, gndsub_corner_ref_top, spacing=feedspacing)
    movebelow(gndsub_vert_section, gndsub_corner_ref_bot, spacing=-feedspacing)
    gndsub_vert_cell.add(gndsub_vert_section)
    gndsub_vert_cell.add(gndsub_corner_ref_top)
    gndsub_vert_cell.add(gndsub_corner_ref_bot)
    gndsub_vert_cell.flatten()

    # Calculations needed to properly center the lines at (0, 0)
    #fcdx, fcdy = get_size(feed_corner)
    #gcdx, gcdy = get_size(gndsub_corner)

    #fhdx, fhdy = get_size(main_hor_section)
    #ghdx, ghdy = get_size(gndsub_hor_section)

    #fvdx, fvdy = get_size(main_vert_section)
    #gvdx, gvdy = get_size(gndsub_vert_section)

    #gndsub_xoffset = (ghdx + gvdx/2)/2
    #gndsub_yoffset = (gvdy + ghdy/2)/2

    #main_xoffset = (fhdx + fvdx/2)/2
    #main_yoffset = (fvdy + fhdy/2)/2

    #(gxmin, _), (_, _) = gndsub_hor_section.get_bounding_box()
    #Dx = ghdx/2 + gxmin
    #gndsub_hor_section.translate(-Dx, 0)
    #(fxmin, _), (_, _) = main_hor_section.get_bounding_box()
    #Dx = fhdx/2 + fxmin
    #main_hor_section.translate(-Dx, 0)

    ##Position the gndsub section properly
    #gndsub_hor_section.translate(-gvdx/4, gndsub_yoffset)
    #main_hor_section.translate(-fvdx/4, main_yoffset)

    gndsubmargin = (feed_cell_width - feed_main_width)/2
    u_main_hor_ref = gdspy.CellReference(main_hor_cell)
    u_main_vert_ref = gdspy.CellReference(main_vert_cell)
    u_gndsub_hor_ref = gdspy.CellReference(gndsub_hor_cell)
    u_gndsub_vert_ref = gdspy.CellReference(gndsub_vert_cell)

    moveright(u_main_hor_ref, u_main_vert_ref, spacing=feed_main_width)
    moveabove(u_main_vert_ref, u_main_hor_ref, spacing=feed_main_width)
    movebelow(u_main_hor_ref, u_gndsub_hor_ref, spacing=gndsubmargin)
    moveright(u_main_vert_ref, u_gndsub_vert_ref, spacing=gndsubmargin)

    print (u_main_hor_ref.get_bounding_box())
    print (u_main_vert_ref.get_bounding_box())
    print (u_gndsub_hor_ref.get_bounding_box())
    print (u_gndsub_vert_ref.get_bounding_box())

    sys.exit()
    l_main_hor_ref = gdspy.CellReference(main_hor_cell)
    l_main_vert_ref = gdspy.CellReference(main_vert_cell)
    l_gndsub_hor_ref = gdspy.CellReference(gndsub_hor_cell)
    l_gndsub_vert_ref = gdspy.CellReference(gndsub_vert_cell)

    moveright(l_main_hor_ref, l_main_vert_ref, spacing=feed_main_width)
    movebelow(l_main_vert_ref, l_main_hor_ref, spacing=-feed_main_width)
    moveright(l_gndsub_hor_ref, l_gndsub_vert_ref, spacing=feed_cell_width)
    movebelow(l_gndsub_vert_ref, l_gndsub_hor_ref, spacing=-feed_cell_width)

    #moveright(gndsub_hor_section, gndsub_corner, spacing=feed_cell_width/2)
    #moveright(main_hor_section, feed_corner, spacing=fcdx/2)
    #centery(gndsub_hor_section, gndsub_corner)
    #centery(main_hor_section, feed_corner)

    #centerx(feed_corner, main_vert_section)
    #centerx(feed_corner, gndsub_vert_section)

    #movebelow(gndsub_corner, main_vert_section, spacing=-feed_cell_width/2)
    #movebelow(gndsub_corner, gndsub_vert_section, spacing=-feed_cell_width/2)

    ufeed_terminus = gdspy.Cell('upper_main_feedline_to_pad')
    ufeed_terminus.add(u_main_hor_ref)
    ufeed_terminus.add(u_main_vert_ref)

    ugndsub_terminus = gdspy.Cell('upper_gndsub_feedline_to_pad')
    ugndsub_terminus.add(u_gndsub_hor_ref)
    ugndsub_terminus.add(u_gndsub_vert_ref)

    uterminus = gdspy.Cell('upper_feedline_to_pad')
    uterminus.add([gdspy.CellReference(ufeed_terminus),\
            gdspy.CellReference(ugndsub_terminus)])

    lfeed_terminus = gdspy.Cell('lower_main_feedline_to_pad')
    lfeed_terminus.add(l_main_hor_ref)
    lfeed_terminus.add(l_main_vert_ref)

    lgndsub_terminus = gdspy.Cell('lower_gndsub_feedline_to_pad')
    lgndsub_terminus.add(l_gndsub_hor_ref)
    lgndsub_terminus.add(l_gndsub_vert_ref)

    lterminus = gdspy.Cell('lower_feedline_to_pad')
    lterminus.add([gdspy.CellReference(lfeed_terminus),\
            gdspy.CellReference(lgndsub_terminus)])
    return uterminus, lterminus

def placeuniquecaps(caps_inv, mask, mcols, nrows, ncols):
    # mrows, mcols are the number of rows and columns on the mask.
    N = nrows * ncols
    # All the capacitors have the same size
    cap_dx, cap_dy = get_size(caps_inv[0])

    # startx = -mask_width/2 + cap_dx/2
    # starty = mask_length/2 - cap_dy/2
    startx = -mask_width/2 + cap_dx/2
    starty = (N // mcols + 1) * cap_dy + cap_dy/2

    caps_inv_ref = [None]*N
    caps_inv_ref[0] = gdspy.CellReference(caps_inv[0],\
            [startx, starty])
    for icap in range(1, N):
        cap_ref = gdspy.CellReference(caps_inv[icap])
        if icap % mcols == 0:
            movebelow(caps_inv_ref[icap - mcols], cap_ref)
            centerx(caps_inv_ref[icap - mcols], cap_ref)
        else:
            moveright(caps_inv_ref[icap - 1], cap_ref)
            centery(caps_inv_ref[icap - 1], cap_ref)
        caps_inv_ref[icap] = cap_ref

    mask.add(caps_inv_ref)
    return caps_inv_ref

def get_via(dx):
    via = gdspy.Cell('Via_to_Ground')

    box = gdspy.Rectangle([-dx/2, dx/2], [dx/2, -dx/2], layer=def_layers["ILD"])
    via.add(box)
    return via

def get_resonator_structures(reso_gnd_sub, reso_ild_sub, feed,  i_xoffset, feed_xoffset, u_dy):
    reso_top = gdspy.Cell('reso_structure')
    g_dx, g_dy = get_size(reso_gnd_sub)
    g_dy -= feedres_dist
    gndsub = gdspy.CellReference(reso_gnd_sub)
    ildsub = gdspy.CellReference(reso_ild_sub)
    ildsub.translate(-i_xoffset, 0)
    
    am = gdspy.CellReference(get_alignment_marks())
    centery(ildsub, am)
    moveleft(gndsub, am, spacing=500)

    cap2feed, cap2gnd = get_coupling_conns(g_dy, u_dy)
    cap2feed_dx, cap2feed_dy = get_size(cap2feed)
    cap2gnd_dx, cap2gnd_dy = get_size(cap2gnd)
    overlap = 2
    cap2feed_yoffset = u_dy/2 + cap2feed_dy/2 - overlap
    cap2gnd_yoffset = -u_dy/2 - cap2gnd_dy/2 + overlap

    cap2feed_ref = gdspy.CellReference(cap2feed)
    cap2feed_ref.translate(-feed_xoffset, cap2feed_yoffset)

    cap2gnd_ref = gdspy.CellReference(cap2gnd)
    cap2gnd_ref.translate(-feed_xoffset, cap2gnd_yoffset)

    via = get_via(cap2gnd_dx*0.5)
    via_ref = gdspy.CellReference(via)
    movebelow(cap2gnd_ref, via_ref, spacing = -0.75*cap2gnd_dx)
    centerx(cap2gnd_ref, via_ref)

    feed_ref = gdspy.CellReference(feed)
    # feed_ref.translate(-feed_xoffset, 0)
    centerx(gndsub, feed_ref)
    f_dx, f_dy = get_size(feed_ref)
    moveabove(gndsub, feed_ref, spacing=f_dy/2)
    reso_top.add(gndsub)
    reso_top.add(ildsub)
    reso_top.add(cap2feed_ref)
    reso_top.add(cap2gnd_ref)
    reso_top.add(via_ref)
    reso_top.add(feed_ref)
    reso_top.add(am)

    return reso_top


def populate_column(sheet, col, startrow, dataset):
    N = len(dataset)
    for index in range(N):
        row = startrow + index
        sheet.cell(column = col, row=row, value = dataset[index])

def get_gnd_cutouts(yspacing):
    dx = 1.5*yspacing
    dy = 500
    margin = 2*inv_margin
    outer = gdspy.Rectangle([-dx/2 - margin/2,dy/2 + margin/2],\
                    [dx/2 + margin/2,-dy/2 - margin/2], layer=def_layers['ILD'])
    inner = gdspy.Rectangle([-dx/2,dy/2], [dx/2,-dy/2],layer=def_layers['ILD'])
    hcutout = gdspy.fast_boolean(outer, inner, 'xor', layer=def_layers['ILD'])
    hfiller = gdspy.Rectangle([-dx/2,dy/2],\
                    [dx/2,-dy/2],layer=def_layers['400nm_NbWiring'])
    hcutout_cell = gdspy.Cell('GP_edge_opening_hor_r')
    hcutout_cell.add(hcutout)
    hfiller_cell = gdspy.Cell('GP_edge_filler_hor')
    hfiller_cell.add(hfiller)
    hgnd_cutout = gdspy.Cell('GP_edge_hor')
    hgnd_cutout.add([gdspy.CellReference(hfiller_cell),\
            gdspy.CellReference(hcutout_cell)])


    outer = gdspy.Rectangle([-dy/2 - margin/2,dx/2 + margin/2],\
                    [dy/2 + margin/2,-dx/2 - margin/2], layer=def_layers['ILD'])
    inner = gdspy.Rectangle([-dy/2,dx/2], [dy/2,-dx/2],layer=def_layers['ILD'])
    vcutout = gdspy.fast_boolean(outer, inner, 'xor', layer=def_layers['ILD'])
    vfiller = gdspy.Rectangle([-dy/2,dx/2],\
                    [dy/2,-dx/2],layer=def_layers['400nm_NbWiring'])
    vcutout_cell = gdspy.Cell('GP_edge_opening_vert_r')
    vcutout_cell.add(vcutout)
    vfiller_cell = gdspy.Cell('GP_edge_filler_vert')
    vfiller_cell.add(vfiller)
    vgnd_cutout = gdspy.Cell('GP_edge_ver')
    vgnd_cutout.add([gdspy.CellReference(vfiller_cell),\
            gdspy.CellReference(vcutout_cell)])

    return hgnd_cutout, vgnd_cutout

def get_island():
    fn = '../resobolo_files/LSN_Island_280umlegs-R.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    isl = gdsii.extract('LSN_Island_280umlegs_R')
    island = gdspy.Cell('LSN_Island_280umlegs')
    isl_ref = gdspy.CellReference(isl)
    island.add(isl_ref)
    island.flatten()
    island.copy('LSN_Island_280umlegs_r')

    return island

def get_XeF2_release():
    fn = '../resobolo_files/XeF2_release.gds'
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(fn)
    xef2 = gdsii.extract('XeF2_release')
    invert_cell(xef2)
    return xef2




def main():
    # Wafer organization all dimensions in microns
    nrows = 8
    ncols = 8
    N_res = nrows * ncols

    wafer = gdspy.Cell('Global_Overlay')
    ind = get_inductor()
    island = get_island()
    xef2 = get_XeF2_release()

    # Generate all the capacitors
    L = 10 #nH
    fstart = 300 # MHz
    df = 3 #MHz
    fs = fstart + df * np.arange(N_res)

    capfrac = 1/(fs/fs[0])**2
    Cs = (1/(L*u.nH)/(2*np.pi*fs*u.MHz)**2).to(u.F).value
    idcs = [IDC(capfrac[i]) for i in range(N_res)]

    # Let's make the coupling capacitors as well
    Qc = Qi
    Ccs = 2 * np.sqrt((Cs * u.F)/((Z0/2) * Qc * 2*np.pi*fs*u.MHz)).to(u.F).value
    coupcaps = [IDC(1.0) for i in range(len(Cs))]

    caps = []
    coup_Nfingers = []
    Nfingers = []
    cfracs = []
    for idc, c, cs, coupcap, freq in zip(idcs, Cs, Ccs, coupcaps, fs):
        idc.finger_length = finger_length
        idc.cellname = "Capacitor_{0:3d}MHz".format(freq)
        nfinger, cfrac = getnumfingers(idc, c)
        Nfingers.append(nfinger)
        cfracs.append(cfrac)

        coupcap.finger_length = coup_cap_finger_length - 2
        nfingers, _ = getnumfingers(coupcap, cs)
        nfingers = int(round(nfingers/5)*5)
        coup_Nfingers.append(nfingers)


    coup_Nfingers = np.array(coup_Nfingers)
    unq_nfingers, unq_indices, unq_inv = np.unique(coup_Nfingers,\
            return_index=True, return_inverse=True)
    unq_coupcaps = []
    for index in unq_indices:
        coupcaps[index].nfinger = coup_Nfingers[index]
        coupcaps[index].cellname =\
                "coupling_capacitor_{0:d}".format(coup_Nfingers[index])
        unq_coupcaps.append(coupcaps[index].draw())

    Nfinger = int(round(np.max(Nfingers)/10)*10)
    common_Nfinger = np.min(Nfingers)
    common_capfrac = capfrac[-1]

    common_idc = IDC(1.0)
    common_idc.nfinger = common_Nfinger
    common_idc.cellname = "Capacitor_common"
    common_cap = common_idc.draw(less_one=True)
    common_cap_ref = gdspy.CellReference(common_cap, rotation=90)
    obtained_Cs = []
    obtained_Ccs = []
    for ires, coup_index in zip(range(N_res), unq_inv):
        idcs[ires].nfinger = Nfinger - common_Nfinger
        idcs[ires].capfrac = (idcs[ires].capfrac - common_capfrac)/(1 - common_capfrac)
        caps.append(get_cap_tank(idcs[ires].draw(), unq_coupcaps[coup_index]))
        obtained_Cs.append(idcs[ires].C)
        obtained_Ccs.append(coupcaps[ires].C)

    # We will need the size of the capacitor when positioning it finally
    cap_ref = gdspy.CellReference(caps[0],rotation=90)
    u_dx, u_dy = get_size(cap_ref)

    # Make the common portion of the resonator
    common_resonator = get_common_resonator(ind, common_cap)

    cres_dx, cres_dy = get_size(common_resonator)
    cres_x_uneven = 1 #xmax - (-xmin) for the resonator

    # Make an array from the common resonator
    xspacing = 11000 #roundto(wafer_width/ncols, 100)
    yspacing = 11000 #roundto(wafer_len/nrows, 100)
    arr_xshift = (xspacing/2)*(ncols-1)
    arr_yshift = (yspacing/2)*(nrows-1)
    common_res_arr = gdspy.CellArray(common_resonator, ncols, nrows, [xspacing, yspacing] )
    common_res_arr.translate(-arr_xshift, -arr_yshift)
    # Make the GND plane subtract region around each resonator
    cres_bbox = common_resonator.get_bounding_box()

    # Generate the region around the resonator where the ground plane is to be
    # removed. Also return the offsets in the x direction from the center of this
    # region where the common resonator and the unique capacitor are to be placed.
    reso_gnd_sub, c_xoffset, u_xoffset, feed_xoffset, gndsub_margin =\
            get_resonator_GPsub(cres_bbox, get_size(cap_ref))
    g_dx, g_dy = get_size(reso_gnd_sub)
    g_dy -= feedres_dist

    reso_ild_sub, i_xoffset = get_resonator_ILDsub(common_cap_ref,\
            cap_ref, u_xoffset, gndsub_margin)

    feed_cell = get_feedline(xspacing)

    reso_structure = get_resonator_structures(reso_gnd_sub, reso_ild_sub,\
            feed_cell, i_xoffset, feed_xoffset, u_dy)
    rs_dx, rs_dy = get_size(reso_structure)
    #Make the feedline
    f_dx, f_dy = get_size(feed_cell)
    fncols = int((xspacing//f_dx)*ncols) - 1
    fnrows = nrows
    fxspacing = xspacing
    fyspacing = yspacing


    reso_struct_arr = gdspy.CellArray(reso_structure, ncols, nrows,\
            [xspacing, yspacing] )
    reso_struct_arr.translate(-arr_xshift, -arr_yshift)
    wafer.add(reso_struct_arr)
    cu_overlap = 2
    common_res_arr.translate(c_xoffset, 0)
    wafer.add(common_res_arr)

    # The position of the bottom left corner of the resonator array which is the
    # origin of the i, j coordinate system
    x_origin = -arr_xshift
    y_origin = -arr_yshift
    # Need to now map the variable capacitor pieces relative to the full array
    for ires in range(N_res):
        i, j = d2ij(ires, nrows)
        xpos = x_origin + i*xspacing - u_xoffset
        ypos = y_origin + j*yspacing
        cap_ref = gdspy.CellReference(caps[ires], [xpos, ypos], rotation=90)
        #cap2feed_ref = gdspy.CellReference(cap2feed)
        #moveabove(cap_ref, cap2feed_ref, spacing=2)
        #cap2gnd_ref = gdspy.CellReference(cap2gnd)
        #movebelow(cap_ref, cap2gnd_ref, spacing=-2)
        wafer.add(cap_ref)
        #wafer.add(cap2feed_ref)
        #wafer.add(cap2gnd_ref)



    vertfeedcell = get_vertical_feed()
    leftconns = gdspy.CellArray(vertfeedcell, 1, nrows//2, [0, 2*fyspacing])
    rightconns = gdspy.CellArray(vertfeedcell, 1, nrows//2-1,[0, 2*fyspacing])

    centery(reso_struct_arr, leftconns)
    centery(reso_struct_arr, rightconns)
    moveleft(reso_struct_arr, leftconns, spacing= -feed_cell_width/2 )
    moveright(reso_struct_arr, rightconns, spacing=feed_cell_width/2 )
    leftconns.translate(0, (rs_dy - f_dy)/2)
    rightconns.translate(0, (rs_dy - f_dy)/2)
    pad = get_bondpads()
    p_dx, p_dy = get_size(pad)
    top_pad_ref = gdspy.CellReference(pad, [0,0])
    bot_pad_ref = gdspy.CellReference(pad, [0,0])


    ulinetopad, llinetopad = getlinetopad()
    upperlinetopad = gdspy.CellReference(ulinetopad)
    lowerlinetopad = gdspy.CellReference(llinetopad)
    uc_dx, uc_dy = get_size(upperlinetopad)
    lc_dx, lc_dy = get_size(lowerlinetopad)

    moveright(reso_struct_arr, upperlinetopad)
    moveright(reso_struct_arr, lowerlinetopad)
    moveabove(reso_struct_arr, upperlinetopad, spacing=uc_dy)
    movebelow(reso_struct_arr, lowerlinetopad,\
            spacing=-(lc_dy + rs_dy - feed_cell_width))
    movebelow(upperlinetopad, top_pad_ref)
    moveright(upperlinetopad, top_pad_ref, spacing=(p_dx + feed_cell_width)/2)
    moveright(lowerlinetopad, bot_pad_ref, spacing=(p_dx + feed_cell_width)/2)
    moveabove(lowerlinetopad, bot_pad_ref)

    wafer.add(leftconns)
    wafer.add(rightconns)
    wafer.add(top_pad_ref)
    wafer.add(bot_pad_ref)
    wafer.add(upperlinetopad)
    wafer.add(lowerlinetopad)

    hcutout, vcutout = get_gnd_cutouts(yspacing)
    cu_dx, cu_dy = get_size(hcutout)
    hor_cutouts = gdspy.CellArray(hcutout, 2, 2,\
            [cu_dx + 2* xspacing, (ncols + 1.08) * yspacing])
    # ver_cutouts.translate(-arr_xshift, -arr_yshift)
    recenter(hor_cutouts)
    (hxmin, hymin), (hxmax, hymax) = hor_cutouts.get_bounding_box()
    x0 = (hxmax + hxmin)/2
    y0 = (hymax + hymin)/2
    hor_cutouts.translate(-x0, -y0)
    wafer.add(hor_cutouts)

    cu_dx, cu_dy = get_size(vcutout)
    ver_cutouts = gdspy.CellArray(vcutout, 2, 2,\
            [(ncols + 1.08) * xspacing, cu_dy + 2* yspacing])

    recenter(ver_cutouts)
    (hxmin, hymin), (hxmax, hymax) = ver_cutouts.get_bounding_box()
    x0 = (hxmax + hxmin)/2
    y0 = (hymax + hymin)/2
    ver_cutouts.translate(-x0, -y0)
    wafer.add(ver_cutouts)


    chip_outline = makechipoutline(wafer_width, wafer_len,'WaferOutline')
    wafer.add(gdspy.CellReference(chip_outline))

    ##########################################################################
    #####       MASK GENERATION. GLOBAL OVERLAY COMPLETED.  ##################
    ##########################################################################


    all_cells = main_lib.cell_dict
    make_inverted_cells()

    inv_cell_list = get_inverted_cells() #+ [ixef2, i_island]

    not_yet = set(inv_cell_list)
    total_mask_area = mask_length * mask_width
    total_area_needed = np.sum(list(map(lambda x: get_cell_area(x),\
            inv_cell_list)))

    mask = gdspy.Cell('Mask')
    maskoutline = makechipoutline(mask_width, mask_length, 'MaskOutline')
    outline_width = 2
    mask.add(gdspy.CellReference(maskoutline))

    # 1. Start with the largest feedline structures
    ivmain = all_cells['vert_main_with_corners_r']
    ivgndsub = all_cells['vert_gndsub_with_corners_r']
    iugndtopad = all_cells['upper_gndsub_feedline_to_pad_r']
    iumaintopad = all_cells['upper_main_feedline_to_pad_r']
    ilgndtopad = all_cells['lower_gndsub_feedline_to_pad_r']
    ilmaintopad = all_cells['lower_main_feedline_to_pad_r']
    ivgndopening = all_cells['GP_edge_opening_vert_r']
    ihgndopening = all_cells['GP_edge_opening_hor_r']
    ivgndfiller = all_cells['GP_edge_filler_vert_r']
    ihgndfiller = all_cells['GP_edge_filler_hor_r']

    ivmain_ref = gdspy.CellReference(ivmain)
    ivgndsub_ref = gdspy.CellReference(ivgndsub)
    iugndtopad_ref = gdspy.CellReference(iugndtopad)
    iumaintopad_ref = gdspy.CellReference(iumaintopad)
    ilgndtopad_ref = gdspy.CellReference(ilgndtopad)
    ilmaintopad_ref = gdspy.CellReference(ilmaintopad)
    ivgndopening_ref = gdspy.CellReference(ivgndopening)
    ihgndopening_ref = gdspy.CellReference(ihgndopening)
    ivgndfiller_ref = gdspy.CellReference(ivgndfiller)
    ihgndfiller_ref = gdspy.CellReference(ihgndfiller)

    go_dx, go_dy = get_size(ihgndopening_ref)
    centerx(maskoutline, ihgndopening_ref)
    centerx(maskoutline, ihgndfiller_ref)
    movebelow(maskoutline, ihgndopening_ref, spacing=-go_dy-outline_width/2)

    moveabove(ihgndopening_ref, ihgndfiller_ref)
    moveabove(ihgndopening_ref, ivgndfiller_ref)
    moveabove(ihgndopening_ref, iumaintopad_ref)
    moveabove(ihgndopening_ref, ilmaintopad_ref)
    moveabove(ihgndopening_ref, iugndtopad_ref)
    moveabove(ihgndopening_ref, ilgndtopad_ref)
    moveabove(ihgndopening_ref, ivgndopening_ref)
    moveabove(ihgndopening_ref, ivmain_ref)
    moveabove(ihgndopening_ref, ivgndsub_ref)

    mp_dx, mp_dy = get_size(iumaintopad)
    moveright(maskoutline, iumaintopad_ref, spacing=mp_dx + outline_width/2)
    moveleft(iumaintopad_ref, iugndtopad_ref)
    moveleft(maskoutline, ilmaintopad_ref, spacing=-mp_dx - outline_width/2)
    moveright(ilmaintopad_ref, ilgndtopad_ref)

    moveleft(iugndtopad_ref, ivgndopening_ref)
    moveright(ilgndtopad_ref, ivmain_ref)
    moveright(ivmain_ref, ivgndsub_ref)
    moveright(ivgndsub_ref, ivgndfiller_ref)


    mask.add(ivmain_ref)
    mask.add(ivgndsub_ref)
    mask.add(iugndtopad_ref)
    mask.add(iumaintopad_ref)
    mask.add(ilgndtopad_ref)
    mask.add(ilmaintopad_ref)
    mask.add(ivgndopening_ref)
    mask.add(ihgndopening_ref)
    mask.add(ihgndfiller_ref)
    mask.add(ivgndfiller_ref)

    not_yet -= set([ivmain, ivgndsub, iugndtopad, iumaintopad,\
      ilgndtopad, ilmaintopad, ivgndopening, ihgndopening])


    # I want to lay the resonator structure cells at the margins of the feedlines
    # Place the resonator structure features
    iresogndsub = all_cells['reso_GP_sub_r']
    iILDsub = all_cells['reso_ILD_sub_r']
    imainfeed = all_cells['feedline_main_r']
    igndsubfeed = all_cells['feedline_GP_sub_r']
    i2feed = all_cells['cap_to_feed_r']
    i2gnd = all_cells['cap_to_gnd_r']
    ivia = all_cells['Via_to_Ground_r']
    iam = all_cells['alignment_marks_patch_new_r']

    iresogndsub_ref = gdspy.CellReference(iresogndsub)
    iILDsub_ref = gdspy.CellReference(iILDsub)
    imainfeed_ref = gdspy.CellReference(imainfeed)
    igndsubfeed_ref = gdspy.CellReference(igndsubfeed)
    i2feed_ref = gdspy.CellReference(i2feed)
    i2gnd_ref = gdspy.CellReference(i2gnd)
    ivia_ref = gdspy.CellReference(ivia)
    iam_ref = gdspy.CellReference(iam)

    ireso_struct = set([iresogndsub, iILDsub, imainfeed,\
      igndsubfeed, i2feed, i2gnd, ivia, iam])

    moveabove(ihgndopening_ref, iresogndsub_ref)
    moveleft(ivgndopening_ref, iresogndsub_ref)
    centery(iresogndsub_ref, iILDsub_ref)
    moveleft(iresogndsub_ref, iILDsub_ref)

    mf_dx, mf_dy = get_size(imainfeed_ref)
    moveabove(iugndtopad_ref, igndsubfeed_ref)#, spacing=mf_dy)
    moveabove(igndsubfeed_ref, imainfeed_ref)
    moveabove(imainfeed_ref, ihgndfiller_ref)

    i2_dx, i2dy = get_size(i2feed_ref)

    moveleft(iILDsub_ref, i2feed_ref)
    moveabove(iILDsub_ref, i2feed_ref, spacing=i2dy)
    centerx(i2feed_ref, i2gnd_ref)
    movebelow(i2feed_ref, i2gnd_ref)
    am_dx, am_dy = get_size(iam_ref)
    moveleft(i2feed_ref, iam_ref)
    moveabove(i2feed_ref, iam_ref, spacing=am_dy)
    moveleft(i2feed_ref, ivia_ref)
    movebelow(iam_ref, ivia_ref)

    mask.add(iresogndsub_ref)
    mask.add(iILDsub_ref)
    mask.add(imainfeed_ref)
    mask.add(igndsubfeed_ref)
    mask.add(i2feed_ref)
    mask.add(i2gnd_ref)
    mask.add(ivia_ref)
    mask.add(iam_ref)

    not_yet -= ireso_struct

    # Placement of the bondpads
    igndbp = all_cells['GNDfeed_bondpad_r']
    imainbp = all_cells['MSfeed_bondpad_r']

    igndbp_ref = gdspy.CellReference(igndbp)
    imainbp_ref = gdspy.CellReference(imainbp)

    moveleft(i2gnd_ref, igndbp_ref)
    movebelow(ivia_ref, igndbp_ref)
    moveleft(igndbp_ref, imainbp_ref)
    centery(igndbp_ref, imainbp_ref)

    mask.add(imainbp_ref)
    mask.add(igndbp_ref)

    not_yet -= set([igndbp, imainbp])

    # Placement of the corner cells
    # igndsubcorner = all_cells['gndsub_corner_r']
    # ifeedcorner = all_cells['feed_corner_r']

    # igndsubcorner_ref = gdspy.CellReference(igndsubcorner)
    # ifeedcorner_ref = gdspy.CellReference(ifeedcorner)

    # centerx(i2feed_ref, ifeedcorner_ref)
    # movebelow(i2gnd_ref, ifeedcorner_ref)
    # moveleft(ivia_ref, igndsubcorner_ref)
    # centery(ivia_ref, igndsubcorner_ref)

    # mask.add(igndsubcorner_ref)
    # mask.add(ifeedcorner_ref)

    # not_yet -= set([igndsubcorner, ifeedcorner])

    # Placement of the common pixel on the mask
    icommoncap = all_cells['Capacitor_common_r']
    i_ind = all_cells['Al_inductor_r']
    icap2ind = all_cells['Cap_to_Ind_lines_r']

    icommoncap_ref = gdspy.CellReference(icommoncap)
    i_ind_ref = gdspy.CellReference(i_ind)
    icap2ind_ref = gdspy.CellReference(icap2ind)

    moveright(icap2ind_ref, i_ind_ref)
    moveleft(icap2ind_ref, icommoncap_ref)

    mask.add(icap2ind_ref)
    mask.add(icommoncap_ref)
    mask.add(i_ind_ref)


    icommon_pixel = set([icommoncap, i_ind, icap2ind])
    not_yet -= icommon_pixel

    # Placement of the unique capacitors on the mask. I'll fit a set of 30 caps above the common pixel
    # another 30 below and 2 flanking each side of the common pixel
    caps_inv = [cell for cell in inv_cell_list \
            if cell.name.startswith('Capacitor') and cell.name.find('common')<0]
    caps_inv.sort(key=lambda x:x.name)

    num_above = 30
    num_below = 30
    num_rflank = 2
    num_lflank = 2

    mrows = 6
    mcols = 5

    num_start, num_end = 0, num_above
    caps_above = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_lflank
    caps_lflank = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_rflank
    caps_rflank = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:num_end]))
    num_start, num_end = num_end, num_end + num_below
    caps_below = list(map(lambda x: gdspy.CellReference(x), caps_inv[num_start:]))

    u_dx, u_dy = get_size(caps_above[0])

    # Place all the caps above the common pixel
    for index, cap in enumerate(caps_above):
        irow, icol = index // mrows, index % mcols
        y_displacement = u_dy * (mrows - irow - 1)
        x_displacement = u_dx * (icol - (mcols//2))
        moveabove(icap2ind_ref, caps_above[index], spacing=-y_displacement)
        caps_above[index].translate(x_displacement, 0)
        mask.add(caps_above[index])

    # Place all the caps below the common pixel
    for index, cap in enumerate(caps_below):
        irow, icol = index // mrows, index % mcols
        y_displacement = u_dy * (mrows - irow - 1)
        x_displacement = u_dx * (icol - (mcols//2))
        movebelow(icap2ind_ref, caps_below[index], spacing=y_displacement)
        caps_below[index].translate(x_displacement, 0)
        mask.add(caps_below[index])

    ci_dx, ci_dy = get_size(icap2ind_ref)
    # Place the right flanking caps
    moveabove(icap2ind_ref, caps_rflank[0])
    moveleft(icap2ind_ref, caps_rflank[0], spacing=-ci_dx/2)
    moveabove(icap2ind_ref, caps_rflank[1])
    moveright(icap2ind_ref, caps_rflank[1], spacing=ci_dx/2)
    mask.add(caps_rflank[0])
    mask.add(caps_rflank[1])

    # Place the left flanking caps
    movebelow(icap2ind_ref, caps_lflank[0])
    moveleft(icap2ind_ref, caps_lflank[0], spacing=-ci_dx/2)
    movebelow(icap2ind_ref, caps_lflank[1])
    moveright(icap2ind_ref, caps_lflank[1], spacing=ci_dx/2)
    mask.add(caps_lflank[0])
    mask.add(caps_lflank[1])

    # iucaps_ref = placeuniquecaps(caps_inv, mask, 10, nrows, ncols)

    # map(lambda x:mask.add(x), caps_above)
    # map(lambda x:mask.add(x), caps_below)
    # map(lambda x:mask.add(x), caps_rflank)
    # map(lambda x:mask.add(x), caps_lflank)

    not_yet -= set(caps_inv)

    # Add the island and xef2 release structure
    ixef2 = all_cells['XeF2_release_r']
    i_island = all_cells['LSN_Island_280umlegs_r']
    isl_dx, isl_dy = get_size(i_island)
    ixef2_ref = gdspy.CellReference(ixef2)
    i_island_ref = gdspy.CellReference(i_island)
    moveleft(icommoncap_ref, i_island_ref)
    moveabove(icommoncap_ref, i_island_ref, spacing=isl_dy)
    moveleft(icommoncap_ref, ixef2_ref)
    movebelow(i_island_ref, ixef2_ref)
    mask.add(i_island_ref)
    mask.add(ixef2_ref)




    main_lib.write_gds('sscale_darkres.gds',unit=1e-6,precision=1e-9)

    # Now I want to generate the full patch shot table from the Global Overlay
    # and the Mask file
    globaloverlay = main_lib.cell_dict['Global_Overlay']
    gcomponents = globaloverlay.elements
    mcomponents = mask.elements
    for component in gcomponents:
        if component.ref_cell.name == "WaferOutline": continue
        allshots.extend(makeshot(component))

    for shot in allshots:
        if shot.cellname.endswith('_r'):
            name = shot.cellname
        else:
            name = shot.cellname + '_r'
        #if name == "WaferOutline_r": continue
        match = list(filter(lambda x: x.ref_cell.name == name, mcomponents))[0]
        shot.update_mask_location(match)

    allshots.sort()

    all_layers = np.array(list(map(lambda x: x.get_layer(), allshots)))
    all_names = np.array(list(map(lambda x: x.get_name(), allshots)))
    all_mask_shifts = np.array(list(map(lambda x: x.get_mask_shift(), allshots)))
    all_cell_sizes = np.array(list(map(lambda x: x.get_cell_size(), allshots)))
    all_cell_shifts = np.array(list(map(lambda x: x.get_cell_shift(), allshots)))
    all_array_sizes = np.array(list(map(lambda x: x.get_array_size(), allshots)))
    all_array_centers = np.array(list(map(lambda x: x.get_array_center(), allshots)))
    all_array_stepsizes = np.array(list(map(lambda x: x.get_array_stepsize(), allshots)))

    # print (all_layers.shape)
    # print (all_names.shape)
    # print (all_mask_shifts.shape)
    # print (all_cell_sizes.shape)
    # print (all_cell_shifts.shape)
    # print (all_array_sizes.shape)
    # print (all_array_centers.shape)
    # print (all_array_stepsizes.shape)

    N = len(all_layers)
    #ids = np.stack([all_layers, all_names], axis=-1)
    #table = np.stack([ids, all_mask_shifts, all_cell_sizes,\
    #    all_cell_shifts, all_array_sizes, all_array_centers,\
    #    all_array_stepsizes], axis=-1)

    xl = all_mask_shifts[:, 0] - all_cell_sizes[:,0]/2
    xr = all_mask_shifts[:, 0] + all_cell_sizes[:,0]/2
    yu = all_mask_shifts[:, 1] + all_cell_sizes[:,1]/2
    yd = all_mask_shifts[:, 1] - all_cell_sizes[:,1]/2

    all_wafer_shifts = all_mask_shifts - all_cell_shifts

    ws = wb.active #get the active worksheet in which to save the data

    ws['A14'] = 'Layer'
    ws['C14'] = 'Cell name'
    ws['D13'] = 'Mask Shift'
    ws.merge_cells('D13:E13')
    ws['D14'] = 'x'
    ws['E14'] = 'y'
    ws['F13'] = 'Cell size'
    ws.merge_cells('F13:G13')
    ws['F14'] = 'x'
    ws['G14'] = 'y'
    ws['M13'] = 'Cell Shift'
    ws.merge_cells('M13:N13')
    ws['M14'] = 'x'
    ws['N14'] = 'y'

    startrow = 15
    endrow = 15 + N

    populate_column(ws, char2num['A'], 15, all_layers)
    populate_column(ws, char2num['C'], 15, all_names)
    populate_column(ws, char2num['D'], 15, all_mask_shifts[:,0])
    populate_column(ws, char2num['E'], 15, all_mask_shifts[:,1])
    populate_column(ws, char2num['F'], 15, all_cell_sizes[:, 0])
    populate_column(ws, char2num['G'], 15, all_cell_sizes[:, 1])
    populate_column(ws, char2num['M'], 15, all_cell_shifts[:, 0])
    populate_column(ws, char2num['N'], 15, all_cell_shifts[:, 1])

    ws.merge_cells('H13:K13')
    ws['H13'] = 'Blade Coordinates'
    ws['H14'] = 'xl'
    ws['I14'] = 'xr'
    ws['J14'] = 'yu'
    ws['K14'] = 'yd'
    populate_column(ws, char2num['H'], 15, xl)
    populate_column(ws, char2num['I'], 15, xr)
    populate_column(ws, char2num['J'], 15, yu)
    populate_column(ws, char2num['K'], 15, yd)

    ws.merge_cells('O13:P13')
    ws['O13'] = 'Wafer Shift'
    ws['O14'] = 'x'
    ws['P14'] = 'y'
    populate_column(ws, char2num['O'], 15, all_wafer_shifts[:,0])
    populate_column(ws, char2num['P'], 15, all_wafer_shifts[:,1])

    ws.merge_cells('T12:Y12')
    ws['T12'] = 'Alternate Array Layout'
    ws.merge_cells('T13:U13')
    ws['T13'] = 'size'
    ws['T14'] = 'C'
    ws['U14'] = 'R'
    ws.merge_cells('V13:W13')
    ws['V13'] = 'center'
    ws['V14'] = 'x'
    ws['W14'] = 'y'
    ws.merge_cells('X13:Y13')
    ws['X13'] = 'step size'
    ws['X14'] = 'x'
    ws['Y14'] = 'y'
    populate_column(ws, char2num['T'], 15, all_array_sizes[:,0])
    populate_column(ws, char2num['U'], 15, all_array_sizes[:,1])
    populate_column(ws, char2num['V'], 15, all_array_centers[:, 0])
    populate_column(ws, char2num['W'], 15, all_array_centers[:, 1])
    populate_column(ws, char2num['X'], 15, all_array_stepsizes[:, 0])
    populate_column(ws, char2num['Y'], 15, all_array_stepsizes[:, 1])










    wb.save(wbsavefile)


if __name__=='__main__':
    main()


