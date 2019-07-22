#! /usr/bin/env python3

import glob
import sys, os
import subprocess

files = glob.glob("*.pdf")

for afile in files:
	name = afile[:-4]
	savename = name + ".eps"
	print ("Converting ...", afile)
	subprocess.call(["pdftops", "-eps", afile, savename])

print ("All files have been converted")
