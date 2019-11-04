#! /usr/bin/env python3

# tree_subdivision.py
# -----------------------------------------------------------------------------#
# Subdivides a grid with particles in an efficient recursive manner to ensure
# that each leaf has only one particle in it.
#
#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import namedtuple
import pdb

point = namedtuple('point', ['x', 'y'])

class bbox():
	def __init__(self, center, halfdim):
		self.center = center
		self.halfdim = halfdim

	def contains_point(self, other):
		inxbound = other.x < self.center.x + self.halfdim and other.x > self.center.x - self.halfdim
		inybound = other.y < self.center.y + self.halfdim and other.y > self.center.y - self.halfdim
		return inxbound and inybound

	def check_1doverlap(p0, p1):
		return not (p0.y < p1.x and p0.x > p1.y)

	def intersects(self, other):
		x0  = point(self.center.x - self.halfdim, self.center.x + self.halfdim)
		y0  = point(self.center.y - self.halfdim, self.center.y + self.halfdim)
		x1  = point(other.center.x - other.halfdim, other.center.x + other.halfdim)
		y1  = point(other.center.y - other.halfdim, other.center.y + other.halfdim)
		xoverlaps = bbox.check_1doverlap(x0, x1)
		yoverlaps = bbox.check_1doverlap(y0, y1)

		return xoverlaps and yoverlaps

N = 64 # NxN grid
n = 15 # Number of particles with n < NxN

class node():
	max_capacity=1
	num_child=4

	def __init__(self, boundary):
		self.boundary = boundary
		# Initialize child nodes to none
		self.nw = None
		self.ne = None
		self.sw = None
		self.se = None
		self.points = []

	def insert_point(self, p):
		if not self.boundary.contains_point(p):
			return False

		if len(self.points) < node.max_capacity and self.nw == None:
			self.points.append(p)
			return True

		if self.nw == None:
			self.subdivide()

		if self.nw.insert_point(p): return True
		if self.ne.insert_point(p): return True
		if self.sw.insert_point(p): return True
		if self.se.insert_point(p): return True

		return False

	# splits the node into 4 child nodes of equal area
	def subdivide(self):
		delta = self.boundary.halfdim/2
		x = self.boundary.center.x
		y = self.boundary.center.y
		self.nw = node(boundary=bbox(point(x - delta, y + delta), delta))
		self.ne = node(boundary=bbox(point(x + delta, y + delta), delta))
		self.sw = node(boundary=bbox(point(x - delta, y - delta), delta))
		self.se = node(boundary=bbox(point(x + delta, y - delta), delta))

	# Returns all the points in the node that are within the range
	def queryrange(self, searchrange):
		ptsinrange = []
		#pdb.set_trace()
		# If the search range doesn't even intersect the boundary, return empty
		# list
		if not self.boundary.intersects(searchrange):
			return ptsinrange

		for pt in self.points:
			if searchrange.contains_point(pt):
				ptsinrange.append(pt)
		# Also check if the child nodes have pts that are in range
		if not self.nw:
			return ptsinrange

		ptsinrange += self.nw.queryrange(searchrange)
		ptsinrange += self.ne.queryrange(searchrange)
		ptsinrange += self.sw.queryrange(searchrange)
		ptsinrange += self.se.queryrange(searchrange)

		return ptsinrange


class nodeplotter(matplotlib.axes.Axes):

	def __init__(self, node, fig, rect=(0, 0, 1, 1), facecolor=None, frameon=True,\
			sharex=None, sharey=None, label='', xscale=None, yscale=None, **kwargs):
		self.node = node
		self.extent = self.node.boundary
		super().__init__(fig, rect, facecolor=None, frameon=True,\
			sharex=None, sharey=None, label='', xscale=None, yscale=None,
			**kwargs)
		super().set_xlim((self.extent.center.x - self.extent.halfdim,
			self.extent.center.x + self.extent.halfdim))
		super().set_ylim((self.extent.center.y - self.extent.halfdim,
			self.extent.center.y + self.extent.halfdim))

	def plot(self):
		xcoords = list(map(lambda pt: pt.x, self.node.points))
		ycoords = list(map(lambda pt: pt.y, self.node.points))
		print (xcoords)
		print (ycoords)
		super().plot(xcoords, ycoords, 'ko', markersize=12)

	

if __name__ == "__main__":
	np.random.seed(15273)
	N = 64#gridsize
	npts = 10
	center = point(0,0)
	region = node(bbox(center, N/2))

	# initialize the random points in the region
	xpoints = np.random.random(npts)*N + center.x - N/2.0
	ypoints = np.random.random(npts)*N + center.y - N/2.0

	pts = map(lambda x: point(*x), zip(xpoints, ypoints))
	for pt in pts:
		region.insert_point(pt)
		#print (pt)

	searchrange = bbox(point(10, 0), 20)
	inrange = region.queryrange(searchrange)
	for pt in inrange:
		print (pt)

	fig = plt.figure(figsize=(10,10))
	display = nodeplotter(region, fig)
	display.plot()
	plt.show()
