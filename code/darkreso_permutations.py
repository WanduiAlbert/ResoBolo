#! /usr/bin/env python3

"""
01/03/2019
------------------------------------------------------------------------------
We did measurements of the first dark resonator array and had a yield of 12/36
detectors. In addition the frequencies have shifted around, maybe collided and
reorganized themselves so it isn't immediately clear what the exact mapping
between resonator seen and resonator on wafer is. This problem can be addressed
by having a way of probing each resonator individually, say using an LED-Mapper
as the guys at NIST did. We don't really have that so I was curious if
brute-force numerical simulation of the possibilities could shed more light. The
probability space has exactly 1,251,677,700 possibilities. Easy enough!

"""

import numpy as np
from math import pi
import matplotlib.pyplot as plt
import time
from itertools import permutations, combinations

def double_range(start, stop):
  while start < stop:
    yield start
    start <<= 1

def rot(n, i, j, ri, rj):
  if (rj == 0):
    if (ri == 1):
      i = (n-1) - i
      j = (n-1) - j
    i, j = j, i
  return (i, j)


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

def round_high(num, base):
	return ((num //base) + ((num % base) > (base//2))) * base

def round_low(num, base):
	return ((num //base) ) * base

# Gives the 2d coordinate of each resonator along an originally 8x8 grid
# threaded by a Hilbert curve with the outer rows deleted to give a 6x6 grid.
def get_reso_positions():
	n = 8
	d = np.arange(n*n)
	positions = np.array(list(map(lambda x: d2ij(x, n), d))) - 1
	mask = np.logical_and(positions[:, 0] >= 0, positions[:, 1] >= 0)
	mask = np.logical_and(mask, positions[:, 0] < 6, positions[:, 1] < 6)
	mask = np.logical_and(mask, positions[:, 1] < 6)

	return positions[mask]
	#print (len(positions[mask]))
	#x = positions[:, 0] - (n-1)/2
	#y = positions[:, 1] - (n-1)/2

if __name__=="__main__":
	freq_meas = np.array([290.6, 303.1, 310.5, 312.9, 328.9, 363.2, 367.5, 382.0
		,386.3, 391.1, 400.4, 448.8, 458.5])
	M = len(freq_meas)
	Delta = 3
	freq_designed = np.array([306, 318, 321, 324, 327, 330, 333, 336, 339, 351, 354,
		369, 372, 381, 384, 387, 390, 393, 396, 399, 402, 405, 408, 417,
		420, 435, 438, 450, 453, 456, 459, 462, 465, 468, 471, 483])
	fmax = freq_designed[-1]
	N = len(freq_designed)
	index_array = np.arange(N)
	reso_pos = get_reso_positions()

	# Approach 1: Purely ascending order in the frequencies along the Hilbert
	# curve
	spacings = round_high(np.cumsum(np.diff(freq_meas)),1)
	high_spacings = round_high(spacings, Delta)
	low_spacings = round_low(spacings, Delta)

	all_subsets = set()
	for i in range(N - M):
		low_samples = np.zeros(M)
		low_samples[0] = freq_designed[i]
		low_samples[1:] = low_samples[0] + low_spacings
		high_samples = np.zeros(M)
		high_samples[0] = freq_designed[i]
		high_samples[1:] = high_samples[0] + high_spacings
		if low_samples[-1] <= fmax:
			all_subsets.add(tuple(low_samples))
		if high_samples[-1] <= fmax:
			all_subsets.add(tuple(high_samples))

	print(all_subsets)
	#exit()

	walker = combinations(index_array, M)
	mean_deviations = []
	std_deviations = []
	start = time.time()
	for i,w in enumerate(walker):
		index = list(w)
		delta = freq_designed[index] - freq_meas
		mean_deviations.append(np.mean(delta))
		std_deviations.append(np.std(delta))
		if i > 2e7: break
	end = time.time()
	print ("Enumerating all arrangements took %1.3f seconds"%(end-start))

	plt.figure(figsize=(10,10))
	plt.hist(mean_deviations, histtype='step')
	plt.xlabel('Mean frequency shift')
	plt.title('Histogram of the mean frequency shift')
	plt.savefig('mean_freqshift.png')

	plt.figure(figsize=(10,10))
	plt.hist(std_deviations, histtype='step')
	plt.xlabel('Stdev of frequency shift')
	plt.title('Histogram of the stdev frequency shift')
	plt.savefig('stdev_freqshift.png')

	plt.show()
