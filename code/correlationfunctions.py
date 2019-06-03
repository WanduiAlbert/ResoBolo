#! /usr/bin/env python3

import numpy as np
import scipy.signal as sig
from math import pi
import matplotlib.pyplot as plt

width = 12.8
height = 12.8
extent = (-width/2, width/2, -height/2, height/2)
pix_size = 0.1
nx, ny = (int(width/pix_size), int(height/pix_size))
patch = np.zeros((nx, ny))


# Want to try and simulate random gaussian fields in 2d
dk_x = 0.1
dk_y = 0.1
nk = 128

kxs = np.arange(nk)*dk_x
kys = np.arange(nk)*dk_y

X, Y = np.meshgrid(np.arange(nx)*pix_size - width/2, np.arange(ny)*pix_size - height/2)

for i in range(nx):
	for j in range(ny):
		# modulate the amplitude somehow
		k = (kxs[i]**2 + kys[j]**2)**0.5
		patch += np.sin(kxs[i]*X + kys[j]*Y + np.random.randn()*pi)


fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(patch, aspect='equal', origin='lower', cmap='Greys', extent=extent)
#ax.plot(Xpos, Ypos, 'k.', ms=20)
ax.set_xlabel('[Deg]')
ax.set_ylabel('[Deg]')
plt.colorbar(img)
plt.savefig('gaussianfield_realspace_constructed.png')

exit()
# Simulating galaxy clusters
cluster_probability = 0.2
cluster_sigma = 0.4
cluster_mean_size = 4

airy_sigma = 0.15
airy_X, airy_Y = np.meshgrid(np.arange(10)-5, np.arange(10)-5)
airy_R =  np.sqrt(airy_X**2 + airy_Y**2)*pix_size
airy_disk = (1/np.sqrt(2*pi*airy_sigma))*np.exp(-airy_R**2/2/airy_sigma**2)

plt.imshow(airy_disk, cmap='Greys')
plt.savefig('airydisk.png')
#exit()
# Mean number of galaxies
n = 50

n_galaxies, = np.random.poisson(n, 1)

Xpos = np.random.random(n_galaxies)*width - width/2
Ypos = np.random.random(n_galaxies)*height - height/2
Xpos[Xpos > width/2] = 0
Xpos[Xpos < -width/2] = 0
Ypos[Ypos > width/2] = 0
Ypos[Ypos < -width/2] = 0

# Need to place each galaxy in the patch based on which pixel it falls into
for i in range(n_galaxies):
	patch[int(Xpos[i]/pix_size), int(Ypos[i]/pix_size)] += 1

	# Also want to add in the clusters around each of the base galaxies we've
	# picked
	inCluster = np.random.rand() < cluster_probability
	if inCluster:
		size, = np.random.poisson(cluster_mean_size, 1)
		cluster_x = np.random.randn(size)*cluster_sigma + Xpos[i]
		cluster_y = np.random.randn(size)*cluster_sigma + Ypos[i]

		for j in range(size):
			patch[int(cluster_x[j]/pix_size), int(cluster_y[j]/pix_size)] += 1

# Zero out the edges of the area
patch[0,:] = 0
patch[:,0] = 0
#patch[-1, :] = 0
#patch[:, -1] = 0
print (np.mean(patch))
print (np.std(patch)**2)

# Convolve the patch with the airy function
patch = sig.convolve2d(patch, airy_disk, mode='same', boundary='symm')

fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(patch, aspect='equal', origin='lower', cmap='Greys', extent=extent)
#ax.plot(Xpos, Ypos, 'k.', ms=20)
ax.set_xlabel('[Deg]')
ax.set_ylabel('[Deg]')
plt.colorbar(img)
plt.savefig('galaxy_random.png')

#autocorr = sig.correlate2d(patch, np.roll(patch, (-10,-10)), mode='same', boundary='symm')
autocorr = sig.correlate2d(patch, patch, mode='same', boundary='symm')


fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(autocorr, aspect='equal', origin='lower', cmap='Greys', extent=extent)
ax.set_xlabel('[Deg]')
ax.set_ylabel('[Deg]')
plt.colorbar(img)
plt.savefig('galaxy_autocorr.png')

#nx, ny = autocorr.shape
#x = np.arange(nx)*pix_size - width/2
#fig, ax = plt.subplots(figsize=(10,10))
#ax.plot(x, np.log(autocorr[:, ny//2]))
#ax.set_xlabel('[Deg]')
#ax.set_ylabel('Log Autocorrelation')
#plt.savefig('autocorr_section.png')


# I want to manually compute the 2 point correlation function between all the
# pixels and then also try and do the 3 and 4 point correlation functions

def two_point_correlation(f):
	nx, ny = f.shape
	correlation = np.zeros((nx, ny))
	mean = np.mean(f)
	mean = 0
	for i in range(nx):
		for j in range(ny):
			for k in range(nx):
				for l in range(ny):
					correlation[i, j] += (f[k, l] - mean)*(f[(k+i)%nx, (l+j)%ny] - mean)
	correlation = np.fft.fftshift(correlation)
	#correlation /= (nx*ny)
	return correlation

print ("Computing the naive 2 point correlation function")
my2point = two_point_correlation(patch)
print (my2point)
print ("All done!!!!\n\n")
fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(my2point, aspect='equal', origin='lower', cmap='Greys', extent=extent)
ax.set_xlabel('[Deg]')
ax.set_ylabel('[Deg]')
plt.colorbar(img)
plt.savefig('my2point.png')

powerspec = np.fft.fftshift(np.fft.fft2(my2point))
kx = np.fft.fftfreq(int(width/pix_size), pix_size)
ky = np.fft.fftfreq(int(height/pix_size), pix_size)
kx_extent = (kx[0], kx[-1], ky[0], ky[-1])
freqgrid = np.meshgrid(kx, ky)

fig, ax = plt.subplots(figsize=(10,10))
img = ax.imshow(np.log(np.abs(powerspec)), aspect='equal', origin='lower', cmap='Greys',
		extent=kx_extent)
ax.set_xlabel('[1/Deg]')
ax.set_ylabel('[1/Deg]')
plt.colorbar(img)
plt.savefig('my2point_powerspec.png')

