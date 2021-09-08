#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from math import pi
from matplotlib.patches import Rectangle, Circle
from matplotlib.collections import PatchCollection
from scipy import interpolate

um = 1e-3
mm = 1

plotdir = 'ledmapper_fig/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)


inductor_length = 325*um
inductor_width = 127*um

led_length = 1.1*mm
led_width = 0.8*mm

#led_length = 3.0*mm
#led_width = 1.5*mm

#collimator_thickness = 1.5875*mm
collimator_thickness = 2.38125*mm
collimator_radius = 0.2667*mm

wafer2collimator = 2.5*mm + collimator_thickness/2.
collimator2mapper = 2.5*mm + collimator_thickness/2.

Nreso = 128


mapper_xcoords = np.array([
       -14.22, -10.98,   2.58,   5.82,  19.38,  22.62, -31.02, -27.78,
       -22.62, -19.38, -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,
        10.98,  14.22,  19.38,  22.62, -31.02, -27.78, -22.62, -19.38,
       -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,
        19.38,  22.62,  27.78,  31.02, -31.02, -27.78, -22.62, -19.38,
       -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,
        19.38,  22.62,  27.78,  31.02, -31.02, -27.78, -22.62, -19.38,
       -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,
        19.38,  22.62,  27.78,  31.02, -31.02, -27.78, -22.62, -19.38,
       -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,
        19.38,  22.62,  27.78,  31.02, -31.02, -27.78, -22.62, -19.38,
       -14.22, -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,
        19.38,  22.62,  27.78,  31.02, -27.78, -22.62, -19.38, -14.22,
       -10.98,  -5.82,  -2.58,   2.58,   5.82,  10.98,  14.22,  19.38,
        22.62,  27.78,  31.02, -22.62, -19.38,  -5.82,  -2.58,  10.98,
        14.22, -31.02, -31.02, -31.02,  31.02,  31.02,  31.02,  22.62])*mm
mapper_ycoords = np.array([
        32.97,  32.97,  32.97,  32.97,  32.97,  32.97,  24.57,  24.57,
        25.14,  25.14,  24.57,  24.57,  25.14,  25.14,  24.57,  24.57,
        25.14,  25.14,  24.57,  24.57,  16.17,  16.17,  16.74,  16.74,
        16.17,  16.17,  16.74,  16.74,  16.17,  16.17,  16.74,  16.74,
        16.17,  16.17,  16.74,  17.58,   6.93,   7.77,   8.34,   8.34,
         7.77,   7.77,   8.34,   8.34,   7.77,   7.77,   8.34,   8.34,
         7.77,   7.77,   8.34,   8.34,  -0.63,  -0.63,  -0.06,  -0.06,
        -0.63,  -0.63,  -0.06,  -0.06,  -0.63,  -0.63,  -0.06,  -0.06,
        -0.63,  -0.63,  -0.06,  0.78,  -9.87,  -9.03,  -8.46,  -8.46,
        -9.03,  -9.03,  -8.46,  -8.46,  -9.03,  -9.03,  -8.46,  -8.46,
        -9.03,  -9.03,  -8.46,  -8.46, -17.43, -17.43, -16.86, -16.86,
       -17.43, -17.43, -16.86, -16.86, -17.43, -17.43, -16.86, -16.86,
       -17.43, -17.43, -16.86, -16.02, -25.83, -25.26, -25.26, -25.83,
       -25.83, -25.26, -25.26, -25.83, -25.83, -25.26, -25.26, -25.83,
       -25.83, -25.26, -25.26, -33.66, -33.66, -33.66, -33.66, -33.66,
       -33.66,   8.61,  -8.19, -25.26,  15.9,  -0.9, -17.7, -33.66])*mm



collimator_xcoords = np.array([-14.07, -11.13,   2.73,   5.67,  19.53,  22.47, -30.87, -27.93,
       -22.47, -19.53, -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,
        11.13,  14.07,  19.53,  22.47, -30.87, -27.93, -22.47, -19.53,
       -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,
        19.53,  22.47,  27.93,  30.87, -30.87, -27.93, -22.47, -19.53,
       -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,
        19.53,  22.47,  27.93,  30.87, -30.87, -27.93, -22.47, -19.53,
       -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,
        19.53,  22.47,  27.93,  30.87, -30.87, -27.93, -22.47, -19.53,
       -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,
        19.53,  22.47,  27.93,  30.87, -30.87, -27.93, -22.47, -19.53,
       -14.07, -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,
        19.53,  22.47,  27.93,  30.87, -27.93, -22.47, -19.53, -14.07,
       -11.13,  -5.67,  -2.73,   2.73,   5.67,  11.13,  14.07,  19.53,
        22.47,  27.93,  30.87, -22.47, -19.53,  -5.67,  -2.73,  11.13,
        14.07, -30.87, -30.87, -30.87,  30.87,  30.87,  30.87,  22.47])*mm
collimator_ycoords = np.array([ 32.97,  32.97,  32.97,  32.97,  32.97,  32.97,  24.57,  24.57,
        25.14,  25.14,  24.57,  24.57,  25.14,  25.14,  24.57,  24.57,
        25.14,  25.14,  24.57,  24.57,  16.17,  16.17,  16.74,  16.74,
        16.17,  16.17,  16.74,  16.74,  16.17,  16.17,  16.74,  16.74,
        16.17,  16.17,  16.74,  16.74,   7.77,   7.77,   8.34,   8.34,
         7.77,   7.77,   8.34,   8.34,   7.77,   7.77,   8.34,   8.34,
         7.77,   7.77,   8.34,   8.34,  -0.63,  -0.63,  -0.06,  -0.06,
        -0.63,  -0.63,  -0.06,  -0.06,  -0.63,  -0.63,  -0.06,  -0.06,
        -0.63,  -0.63,  -0.06,  -0.06,  -9.03,  -9.03,  -8.46,  -8.46,
        -9.03,  -9.03,  -8.46,  -8.46,  -9.03,  -9.03,  -8.46,  -8.46,
        -9.03,  -9.03,  -8.46,  -8.46, -17.43, -17.43, -16.86, -16.86,
       -17.43, -17.43, -16.86, -16.86, -17.43, -17.43, -16.86, -16.86,
       -17.43, -17.43, -16.86, -16.86, -25.83, -25.26, -25.26, -25.83,
       -25.83, -25.26, -25.26, -25.83, -25.83, -25.26, -25.26, -25.83,
       -25.83, -25.26, -25.26, -33.66, -33.66, -33.66, -33.66, -33.66,
       -33.66,   8.34,  -8.46, -25.26,  16.17,  -0.63, -17.43, -33.66])*mm

inductor_xcoords = collimator_xcoords
inductor_ycoords = collimator_ycoords

def plot_circles(x,y,ax,r,c):
    holes = []
    x = np.asarray(x)
    y = np.asarray(y)
    for i, (xp, yp) in enumerate(zip(x, y)):
        hole = Circle((xp, yp),
                r, fill=True)
        holes.append(hole)

    pc = PatchCollection(holes, facecolor=c, alpha=0.5,
            edgecolor=c)
    ax.add_collection(pc)


def plot_rects(x,y,ax, w, h, c):
    inductors = []
    for i, (xp, yp) in enumerate(zip(x, y)):
        rect = Rectangle((xp - w/2., yp - h/2.),
                w, h, fill=True)
        inductors.append(rect)

    pc = PatchCollection(inductors, facecolor=c, alpha=0.5,
            edgecolor=c)
    ax.add_collection(pc)

fig, ax = plt.subplots(figsize=(10,10))
plot_rects(inductor_xcoords, inductor_ycoords, ax, inductor_length,
        inductor_width, 'blue')
plot_circles(collimator_xcoords, collimator_ycoords, ax, collimator_radius,
        'red')
plot_rects(mapper_xcoords, mapper_ycoords, ax, led_length,
        led_width, 'green')
ax.set_xlim(left=-40.0, right=40.0)
ax.set_ylim(bottom=-40.0, top=40.0)
ax.axis('square')
plt.savefig(plotdir + 'wafer_and_LED_alignment.png')
plt.show()
plt.close('all')
#exit()

Thetas = np.zeros((Nreso, Nreso))
Phis = np.zeros((Nreso, Nreso))
Psis = np.zeros((Nreso, Nreso))

for itheta in range(Nreso):
    d = np.sqrt((collimator_xcoords - mapper_xcoords[itheta])**2 +\
            (collimator_ycoords - mapper_ycoords[itheta])**2)
    Thetas[itheta, :] = np.arctan2(d, collimator2mapper)
    Phis[itheta, :] = np.arctan2(collimator_ycoords - mapper_ycoords[itheta],
            collimator_xcoords - mapper_xcoords[itheta])
    lim = 2*collimator2mapper*collimator_radius/collimator_thickness

    num = collimator_radius*(2*collimator2mapper + collimator_thickness)
    denom = d**2 + collimator2mapper**2 - collimator_radius**2 +\
            collimator2mapper*collimator_thickness + collimator_thickness**2/4.
    mask = d <= collimator_radius
    Psis[itheta, mask] = 0.5*np.arctan2(num, denom[mask])

    num = 2*collimator_radius*collimator2mapper - collimator_thickness*d
    denom = d**2 + collimator2mapper**2 - collimator_radius**2 - collimator_thickness**2/4.
    mask = np.logical_and(d > collimator_radius, d <= lim)
    Psis[itheta, mask] = 0.5*np.arctan2(num[mask], denom[mask])

    Psis[itheta, d > lim] = np.nan


#alphapluses[itheta, :] = np.arccos((collimator2mapper*(collimator2mapper +
#    collimator_thickness/2.) + d**2)/np.sqrt((collimator2mapper**2 +
#        d**2)*((collimator2mapper + collimator_thickness/2.)**2 + d**2)))
#alphaminuses[itheta, :] = np.arccos((collimator2mapper*(collimator2mapper -
#    collimator_thickness/2.) + d**2)/np.sqrt((collimator2mapper**2 +
#        d**2)*((collimator2mapper - collimator_thickness/2.)**2 + d**2)))
#au = np.sqrt((collimator2mapper - collimator_thickness/2.)**2 + (d -
#    collimator_radius)**2)
#bu = np.sqrt((collimator2mapper - collimator_thickness/2.)**2 + (d +
#    collimator_radius)**2)
#al = np.sqrt((collimator2mapper + collimator_thickness/2.)**2 + (d -
#    collimator_radius)**2)
#bl = np.sqrt((collimator2mapper + collimator_thickness/2.)**2 + (d +
#    collimator_radius)**2)
#zetauppers[itheta, :] = 0.5*np.arccos(((collimator2mapper -
#    collimator_thickness/2.)**2 + (d**2 - collimator_radius**2))/(au*bu))
#zetalowers[itheta, :] = 0.5*np.arccos(((collimator2mapper -
#    collimator_thickness/2.)**2 + (d**2 - collimator_radius**2))/(al*bl))
#chi = collimator_radius - collimator_thickness*np.tan(Thetas + alphaminuses - zetalowers)
#nu = collimator_thickness*np.tan(Thetas - alphapluses + zetauppers)
#
#mask_upper = chi < collimator_radius
#mask_lower = nu < 2*collimator_radius
#mask = np.logical_and(mask_upper, mask_lower)
#thetamean = Thetas + 0.5*(zetalowers - zetauppers + alphaminuses - alphapluses)
#shift = (collimator_radius - chi)/2.
#rs = (collimator_radius + chi)/2.*(1 + wafer2collimator/collimator2mapper)

rs = (collimator2mapper + wafer2collimator)*np.sqrt(
        np.sin(2*Psis)/2/(np.cos(2*Thetas) + np.cos(2*Psis)))
mask = np.logical_not(np.isnan(rs))
xps = mapper_xcoords[:, np.newaxis] +\
        (collimator2mapper + wafer2collimator)*np.tan(Thetas)*np.cos(Phis)
yps = mapper_ycoords[:, np.newaxis] +\
        (collimator2mapper + wafer2collimator)*np.tan(Thetas)*np.sin(Phis)


plt.figure(figsize=(10,10))
img = plt.imshow(mask, cmap='Greys')
plt.savefig(plotdir + 'illumination_mask.png')
plt.colorbar()
#plt.show()

plt.figure(figsize=(10,10))
plt.imshow(Thetas)
plt.savefig(plotdir + 'thetas_masked.png')
plt.colorbar(label='Theta')
#plt.show()

plt.figure(figsize=(10,10))
plt.imshow(Phis)
plt.savefig(plotdir + 'phis_masked.png')
plt.colorbar(label='Phi')
#plt.show()

plt.figure(figsize=(10,10))
plt.imshow(Psis)
plt.savefig(plotdir + 'psis_masked.png')
plt.colorbar(label='Psi')
plt.show()
plt.close('all')


def calculate_flux(x0, y0, theta0, phi0, z0, rmax):
    xx, yy = np.r_[-rmax:rmax:2000j] + x0, np.r_[-rmax:rmax:2000j] + y0
    XX, YY = np.meshgrid(xx, yy)
    Rho = np.sqrt((XX - x0)**2 + (YY - y0)**2)
    Phi = np.arctan2(YY - y0, XX - x0)
    #xx = x0 + rho*np.cos(phi)
    #yy = y0 + rho*np.sin(phi)
    #Rho, Phi = np.meshgrid(np.r_[0:rmax:2000j], np.r_[0:2*pi:2000j])
    #XX = x0 + Rho*np.cos(Phi)
    #YY = y0 + Rho*np.sin(Phi)

    flux = 1#./z0**2
    flux /= (1 + (Rho/z0)**2 + np.tan(theta0)**2 +
            2*(Rho/z0)*np.tan(theta0)*np.cos(phi0 - Phi))**4
    flux[Rho > rmax] = 0
    return xx, yy, XX, YY, flux

for i in range(Nreso):
    fig, ax = plt.subplots(figsize=(10,10))
    inductor = []
    #plot_rects(, inductor_ycoords, ax, inductor_length,
    #        inductor_width, 'blue')
    ax.plot(mapper_xcoords[i], mapper_ycoords[i], color='r', marker='+', ms=15,
            ls='None')
    interpolators = []
    xmin, ymin, xmax, ymax = inductor_xcoords[i], inductor_ycoords[i],\
            inductor_xcoords[i], inductor_ycoords[i]
    for j in range(Nreso):
        rect = Rectangle((inductor_xcoords[j] - inductor_length/2.,
            inductor_ycoords[j] - inductor_width/2.), inductor_length,
            inductor_width, fill=False)
        inductor.append(rect)
        if not mask[i, j]: continue
        print ("\t\t", i, j)
        xx, yy, XX, YY, flux = calculate_flux(xps[i, j], yps[i, j], Thetas[i,j], Phis[i,j],
                (collimator2mapper + wafer2collimator), rs[i, j])
        if np.min(xx) <= xmin: xmin = np.min(xx)
        if np.min(yy) <= ymin: ymin = np.min(yy)
        if np.max(xx) >= xmax: xmax = np.max(xx)
        if np.max(yy) >= ymax: ymax = np.max(yy)
        #func = interpolate.interp2d(xx, yy, flux, kind='cubic', fill_value=0.0)
        func = interpolate.RectBivariateSpline(xx, yy, flux)
        interpolators.append(func)
        ZZ = func(xx-0.001, yy-0.001)
    xx, yy = np.r_[xmin:xmax:2000j], np.r_[ymin:ymax:2000j]
    print ("\tparameters: ", i, xmin, xmax, ymin, ymax, flux.shape)
    fullXX, fullYY = np.meshgrid(xx, yy)
    fullflux = np.zeros_like(XX)
    for ifn, interpfn in enumerate(interpolators):
        #if ifn >0: continue
        fullflux += interpfn(xx, yy)
    mesh = ax.pcolormesh(fullXX, fullYY, fullflux, cmap='gray_r', shading='auto',
            alpha=0.5, vmin=0, vmax=1)
    #plot_circles(collimator_xcoords, collimator_ycoords, ax, collimator_radius, 'red')
    pc = PatchCollection(inductor, facecolor='none', alpha=0.5,
            edgecolor='b', linestyle='--', linewidth=2)
    ax.add_collection(pc)
    ax.axis('square')
    plt.colorbar(mesh)
    plt.savefig(plotdir + 'point_illumination_LED%d.png'%(i))
    #plt.show()
    plt.close()









