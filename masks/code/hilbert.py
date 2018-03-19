#! /usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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

if __name__=="__main__":
  n = 8
  d = np.arange(n*n)
  positions = np.array(list(map(lambda x: d2ij(x, n), d)))

  print (positions.shape)
  x = positions[:, 0] - (n-1)/2
  y = positions[:, 1] - (n-1)/2

  fig, ax = plt.subplots(figsize=(10,10))
  ax.plot(x, y, 'r-')
  ax.scatter(x, y, s=40, color='b', marker='o')
  ax.set_xlabel('X Pos on Wafer')
  ax.set_ylabel('Y Pos on Wafer')
  plt.savefig('hilbert.png')
  plt.show()