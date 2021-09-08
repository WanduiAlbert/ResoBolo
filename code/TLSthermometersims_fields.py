#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import pi
from scipy import interpolate
from scipy.constants import h, k, c, epsilon_0, mu_0
from scipy import interpolate, optimize
import scipy.signal as signal
import glob
import sys,os
import pdb
#sys.path.append("/home/wanduialbert/bicep/code/instruments/copper")
import reso_fit
import meshanalysisofresonator as mesh
pi = np.pi

MHz = 1e6
GHz = 1e9
kHz = 1e3
pW = 1e-12
um = 1e-6
nm = 1e-9

nH = 1e-9
pH = 1e-12
pF = 1e-12
MHz = 1e6
Z0 = 50
Y0 = 1./Z0

datadir = '../numerical_sims/'
plotdir = 'TLSthermometerfield_fig/'

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'][1:]

plot_diagnostic = True
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

er_sio2 = 3.9
er_si = 11.9
er_air = 1

def load_data(fn, nx, ny, nlevels=5):
    ilevel = 0
    dataset = np.zeros((nlevels, nx, ny), dtype=np.complex128)
    with open(fn) as fp:
        print ("Loading File....")
        for i in range(12): fp.readline()
        iy = 0
        while True:
            line = fp.readline()
            if line.startswith('VER'):
                for i in range(11): fp.readline()
                ilevel += 1
                print (ilevel, iy)
                iy = 0
                line = fp.readline()
            if line is None: break
            line = line.replace(' ', '')
            line = line.replace('\n', '')
            if line is "": continue
            line = line.split(',')[1:]
            dataset[ilevel, :, iy] = np.array([complex(_) for _ in line])
            iy += 1
            if iy == ny and ilevel + 1 == nlevels: break
        print ("Exiting File....")
    return dataset


if __name__=="__main__":
    filenames = glob.glob(datadir +
            "TKID_Module_Thermometer_ShortedCapacitor_samesidecoupling_*.csv")
    filenames.sort(key=lambda x: x.split("/")[-1].split('.')[0].split('_')[-1])
    #for fn in filenames:
    #    print (fn)
    #exit()
    mesh_dim = 2*um
    f0 = 747.051*MHz
    x = np.arange(1, 2851,2)
    y = np.arange(3999, 0, -2)
    X, Y = np.meshgrid(x,y, indexing='ij')
    labels = ['JX', 'JY', 'charge']
    metal_t = np.array([400, 0, 120, 0, 0])*nm
    dielectric_t = np.array([500, 0.15, 0.15, 500, 1000])*um
    kinetic_inductance = np.array([0.158, 0, 0.138, 0, 0])*pH#/sq
    er = np.array([1, 3.9, 3.9, 11.9, 1])

    JX = load_data(filenames[0], x.size, y.size, nlevels=5)
    JY = load_data(filenames[1], x.size, y.size, nlevels=5)
    charge = load_data(filenames[2], x.size, y.size, nlevels=5)
    charge_density = charge/mesh_dim

    EtX = JX*(1j*2*pi*f0*kinetic_inductance[:, np.newaxis, np.newaxis]**2/mu_0)
    EtY = JX*(1j*2*pi*f0*kinetic_inductance[:, np.newaxis, np.newaxis]**2/mu_0)
    Ez = charge_density
    makerawmaps = True

    if makerawmaps:
        plt.figure(1, figsize=(15,10))
        plt.pcolormesh(X, Y, charge[0].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re Q [C/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"charge_density_real_level0.png")
        plt.savefig(path, dpi=400)

        plt.figure(2, figsize=(15,10))
        plt.pcolormesh(X, Y, charge[0].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im Q [C/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"charge_density_imag_level0.png")
        plt.savefig(path, dpi=400)

        plt.figure(3, figsize=(15,10))
        plt.pcolormesh(X, Y, JX[0].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re JX [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_x_real_level0.png")
        plt.savefig(path, dpi=400)

        plt.figure(4, figsize=(15,10))
        plt.pcolormesh(X, Y, JX[0].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im JX [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_x_imag_level0.png")
        plt.savefig(path, dpi=400)

        plt.figure(5, figsize=(15,10))
        plt.pcolormesh(X, Y, JY[0].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re JY [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_y_real_level0.png")
        plt.savefig(path, dpi=400)

        plt.figure(6, figsize=(15,10))
        plt.pcolormesh(X, Y, JY[0].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im JY [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_y_imag_level0.png")
        plt.savefig(path, dpi=400)
        plt.show()
        plt.close('all')

        plt.figure(7, figsize=(15,10))
        plt.pcolormesh(X, Y, charge[2].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re Q [C/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"charge_density_real_level2.png")
        plt.savefig(path, dpi=400)

        plt.figure(8, figsize=(15,10))
        plt.pcolormesh(X, Y, charge[2].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im Q [C/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"charge_density_imag_level2.png")
        plt.savefig(path, dpi=400)

        plt.figure(9, figsize=(15,10))
        plt.pcolormesh(X, Y, JX[2].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re JX [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_x_real_level2.png")
        plt.savefig(path, dpi=400)

        plt.figure(10, figsize=(15,10))
        plt.pcolormesh(X, Y, JX[2].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im JX [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_x_imag_level2.png")
        plt.savefig(path, dpi=400)

        plt.figure(11, figsize=(15,10))
        plt.pcolormesh(X, Y, JY[2].real, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Re JY [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_y_real_level2.png")
        plt.savefig(path, dpi=400)

        plt.figure(12, figsize=(15,10))
        plt.pcolormesh(X, Y, JY[2].imag, cmap='viridis', shading='auto')
        plt.grid()
        plt.colorbar(label='Im JY [A/m]')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')
        #plt.xlim(left=800, right=2000)
        #plt.ylim(bottom=200, top=3800)
        path = os.path.join(plotdir,"surface_current_y_imag_level2.png")
        plt.savefig(path, dpi=400)
        plt.show()


