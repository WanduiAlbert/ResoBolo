#! /usr/bin/env python3

# tree_subdivision.py
# -----------------------------------------------------------------------------#
# Subdivides a grid with particles in an efficient recursive manner to ensure
# that each leaf has only one particle in it.
#
#
import numpy as np
import matplotlib.pyplot as plt

N = 64 # NxN grid
n = 15 # Number of particles with n < NxN

class bbox():


class node():
	max_capacity=1
	num_child=4
	def __init__(self, parent=None):
		self.parent = parent
		s
