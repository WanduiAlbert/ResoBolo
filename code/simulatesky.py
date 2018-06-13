#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy import fftpack
from math import pi
from matplotlib import cm
import scipy.signal as sig
import scipy.interpolate as interp

def randn_ft_of_real(N):
    # Makes an array of N normally distributed numbers with unit variance.
    # Ensures that the fourier transform of the numbers is real.
    m = N-1
    n = N//2-1
    a = (np.random.randn(m,n) + 1j* np.random.randn(m,n))/2**0.5
    b = (np.random.randn(n) + 1j* np.random.randn(n))/2**0.5
    c = (np.random.randn(n) + 1j* np.random.randn(n))/2**0.5
    d = (np.random.randn(1,n) + 1j* np.random.randn(1,n))/2**0.5
    # Useful for generating arrays of complenp.random.randn() random numbers with a unit power
    # spectrum i.e. mean(abs^2) = 1
    rnd = np.zeros((N,N), dtype=np.complex128)
    # Now let's fill out the array appropriately
    rnd[1:, 1:N//2] = a
    rnd[1:, N//2+1:] = np.rot90(np.conj(a), 2)

    rnd[1:N//2, 0] = b
    rnd[N//2+1:, 0] = np.flipud(np.conj(b))

    rnd[1:N//2, N//2] = c
    rnd[N//2+1:, N//2] = np.flipud(np.conj(c))

    rnd[0, 1:N//2] = d
    rnd[0, N//2+1:] = np.fliplr(np.conj(d))

    rnd[0, 0] = np.random.randn()
    rnd[N//2, 0] = np.random.randn()
    rnd[0, N//2] = np.random.randn()
    rnd[N//2, N//2] = np.random.randn()

    return rnd

def calc_axesdata(field_size_deg, npix):
    # ad=calc_ad(field_size_deg,npix)
    #
    # Make a structure containing all the axis data for an image
    # in both image and fourier planes
    #
    # NB: Input Field_size in DEGREES
    #
    # Output dict has the following elements:
    #
    # field_size_deg    (input)
    # npix              (input)
    # field_size        (radians)
    # t_del u_del       (spacings in image and Forier planes)
    # t_val u_val l_val (axis values in radians)
    # u_r               (grid of radial values)
    # t_val_deg         (useful for plots)

    # Keep inputs in output structure too
    ad = dict()
    ad['field_size_deg']=field_size_deg
    ad['npix']=npix

    # Most things in radians
    ad['field_size']=field_size_deg*(pi/180)

    # Calc the spacing in image space
    ad['del_t']=ad['field_size']/npix

    # Calc the spacing in Fourier space - reciprocal of span is resolution
    ad['del_u']=1/ad['field_size']

    # The image space axis values
    ad['t_val']=np.arange(-(npix/2)*ad['del_t'], (npix/2)*ad['del_t'], ad['del_t'])

    # The Fourier space axis values
    ad['u_val']=np.arange(-(npix/2)*ad['del_u'],(npix/2)*ad['del_u'], ad['del_u'])

    # The ell space axis values
    ad['l_val']=ad['u_val']*2*pi

    # Make grid of radial values in fourier plane
    [x,y]=np.meshgrid(ad['u_val'],ad['u_val'])
    ad['u_r']=np.sqrt(x**2+y**2)

    # Image space axis values in degrees useful for plots
    ad['t_val_deg']=ad['t_val']*(180/pi)

    return ad

def gaussian(p, x):
    """
        y = gaussian(p, x)

        Compute values at x on a Gaussian shaped curve with parameters:

        p[0] = peak height
        p[1] = mean
        p[2] = sigma

    """
    return p[0] * np.exp(-((x - p[1])/p[2])**2)


# 1.
# Keep track of the image and fourier plane axes
# a 20x20 chunk of sky with 256x256 pixels
# The fourier space analog of the angle theta is u
# In the small sky approximation l = 2 pi u
npix = 2048
field = 20
ad = calc_axesdata(field, npix)
uval_min = ad['u_val'][0]
uval_max = ad['u_val'][-1]
tmax = ad['field_size_deg']/2
tmin = -tmax

# 2.
# Get some CMB power spectrum from CAMB. First column is ell, then TT, TE, EE,
# BB
dataset = np.loadtxt('camb_spectrum.txt').T
# We will get rid of the monopole and dipole terms
l = dataset[0, 2:]
Cs_l = dataset[1:4, 2:]
prefactor = 2*pi/l/(l+1)
C_l = Cs_l*prefactor


f = randn_ft_of_real(ad['npix'])

image = fftpack.fft2(fftpack.fftshift(f)).real
fig1, ax1 = plt.subplots(figsize=(10,10))
img = ax1.imshow(image, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax1.set_xlabel(r'Angle [Deg]')
ax1.set_ylabel(r'Angle [Deg]')
fig1.colorbar(img)
ax1.set_title('Simulated Sky with Flat Power Spectrum')
plt.savefig('flat_powerspectrum_simulated_sky.png')

# We need to use interpolation to map from the integer ells from CAMB to the
# continuous l from the image we generated
ll = 2*pi*ad['u_r']
lmin = np.min(ll)
lmax = np.max(ll)
print ("Probing l in the range %3.1f to %3.1f" %(lmin, lmax))

TT_intp = np.interp(2*pi*ad['u_r'], l, np.sqrt(C_l[0, :]))
fp = f*TT_intp
Ti = fftpack.fft2(fftpack.fftshift(fp)).real
fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(Ti, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax.set_xlabel(r'Angle [Deg]')
ax.set_ylabel(r'Angle [Deg]')
fig.colorbar(img)
ax.set_title('Simulated TT CMB Sky')
plt.savefig('TT_simulated_sky.png')

# Yeah! Let's look at polarization
# Cosmology predicts E and B modes, whereas we actually measure Q and U on the
# sky. As a result, we need to in Fourier space, rotate the E and B modes into
# the Q and U modes.
#f = randn_ft_of_real(ad['npix'])
#ll = 2*pi*ad['u_r']
Ef = f*np.interp(ll, l, np.sqrt(C_l[2, :]))
Bf = np.zeros_like(Ef) # Assume that there are no b-modes

[u, v] = np.meshgrid(ad['u_val'], ad['u_val'])
chi = -np.arctan2(v, u) + pi/2

Qf = Ef*np.cos(2*chi) + Bf*np.sin(2*chi)
Uf = Ef*np.sin(2*chi) - Bf*np.cos(2*chi)

# Now we can render the Q and U maps
Qi = fftpack.fft2(fftpack.fftshift(Qf)).real
Ui = fftpack.fft2(fftpack.fftshift(Uf)).real

fig2, ax2 = plt.subplots(figsize=(10,10))
img = ax2.imshow(Qi, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax2.set_xlabel(r'Angle [Deg]')
ax2.set_ylabel(r'Angle [Deg]')
fig2.colorbar(img)
ax2.set_title('Simulated Q for E Modes')
plt.savefig('QE_simulated_sky.png')

fig3, ax3 = plt.subplots(figsize=(10,10))
img = ax3.imshow(Ui, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax3.set_xlabel(r'Angle [Deg]')
ax3.set_ylabel(r'Angle [Deg]')
fig3.colorbar(img)
ax3.set_title('Simulated U for E Modes')
plt.savefig('UE_simulated_sky.png')
#plt.show()

# Now let's include a telescope scanning strategy
nx, ny = 200, 33
x = np.linspace(-4,4,nx/2)
x = np.hstack([x, np.flip(x, 0)])
y = np.linspace(-4,4,ny)
x = np.tile(x, [ny,1]).T
y = np.tile(y, [1, nx])
x = x.flatten()
y = y.flatten()

#plt.close()
#fig2, ax2 = plt.subplots(figsize=(10,10))
#img = ax2.imshow(Qi, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
#ax2.set_xlabel(r'Angle [Deg]')
#ax2.set_ylabel(r'Angle [Deg]')
#fig2.colorbar(img)
#ax2.set_title('Simulated Q for E Modes')
#ax2.plot(x, y, 'r-')
#ax2.set_xlim([-5, 5])
#ax2.set_ylim([-5, 5])
#plt.show()

# We need gaussian smoothing to mimic the effect of the beam and prevent
# aliasing
n = ad['npix']//2 + 1
xv = ad['t_val_deg'][n-15:n+15]
[xg, yg] = np.meshgrid(xv, xv)
r = np.sqrt(xg**2 + yg**2)
beam = gaussian([1, 0, 0.5/2.35], r)
beam /= np.sum(beam) # So the beam doesn't change the normalization

# Now we convolve the maps with the beams
Ts = sig.convolve2d(Ti, beam, 'same')
Qs = sig.convolve2d(Qi, beam, 'same')
Us = sig.convolve2d(Ui, beam, 'same')

fig4, ax4 = plt.subplots(figsize=(10,10))
img = ax4.imshow(Ts, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax4.set_xlabel(r'Angle [Deg]')
ax4.set_ylabel(r'Angle [Deg]')
fig4.colorbar(img)
ax4.set_title('Simulated TT with Gaussian Smoothing')
plt.savefig('TT_smoothed_sky.png')

fig5, ax5 = plt.subplots(figsize=(10,10))
img = ax5.imshow(Qi, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax5.set_xlabel(r'Angle [Deg]')
ax5.set_ylabel(r'Angle [Deg]')
fig5.colorbar(img)
ax5.set_title('Simulated Q for EE with Gaussian Smoothing')
plt.savefig('QE_smoothed_sky.png')

fig6, ax6 = plt.subplots(figsize=(10,10))
img = ax6.imshow(Ui, cmap=cm.gray, extent=[tmin,tmax,tmin,tmax])
ax6.set_xlabel(r'Angle [Deg]')
ax6.set_ylabel(r'Angle [Deg]')
fig6.colorbar(img)
ax6.set_title('Simulated U for EE with Gaussian Smoothing')
plt.savefig('UE_smoothed_sky.png')
plt.show()

# We can now make a time stream by interpolating over the maps we have made. It
# turns out that the detector timestream  a(t) = T(r(t)) + Q(r(t)) cos (2 Psi) +
# U(r(t)) sin (2 Psi). Psi is the orientation angle, of the polarized detector.
# This angle is usually measured from the celestial North Pole. We also assume
# that Psi is constant with time
psiA = pi/4

As = Ts + Qs*np.cos(2*psiA) + Us*np.sin(2*psiA)

interpfnA = interp.RectBivariateSpline(ad['t_val_deg'], ad['t_val_deg'], As)
at = interpfnA.ev(x, y)
#add some noise
at += 3*np.random.randn(*at.shape)
# Now make the time stream for the orthogonal detector
psiB = psiA + pi/2
Bs = Ts + Qs*np.cos(2*psiB) + Us*np.sin(2*psiB)
interpfnB = interp.RectBivariateSpline(ad['t_val_deg'], ad['t_val_deg'], Bs)
bt = interpfnB.ev(x, y)
bt += 3*np.random.randn(*bt.shape)

plt.close()
fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(at, label='horizontal')
ax.plot(bt, 'r', label='vertical')
ax.axis('tight')
ax.set_title('Detector Time streams')
#ax.set_xlim([0, 400])
ax.legend(loc='best')
plt.savefig('detector_timestream.png')
plt.show()
plt.close()
