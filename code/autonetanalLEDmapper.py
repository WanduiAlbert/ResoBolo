#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import serial
import io
import time
import pdb
import argparse
import sys, os, datetime, time
import numpy as np
import copper
sys.path.append('/home/lebicep/pyhk')
from pyhkcmd import adr
from ledmapper import LEDMapper
sys.path.append('/home/lebicep/instruments/rs')
import kid

vna = copper.Copper("localhost")

cold_centers = [305.8, 318.3, 337.4, 351.5, 414.6, 428.9, 434.1, 451.4, 473.4, 671.7]
span = 2.0
powers = np.linspace(-30,0,16) - 10
#powers = [-4]
magnet_currents = np.array([0.])
magnet_slew_delay = 1200.
print "magnet_currents: ",magnet_currents
if_bw = 1000
num_pts = 3201

vna.init()
vna.set_num_pts(num_pts)
vna.set_if_bw(if_bw)

mapper = LEDMapper('/dev/cu.usbmodem14101', 19200)
mapper.init()

out_folder = 'raw'
if not os.path.exists(out_folder):
	os.mkdir(out_folder)

def get_datestr():
	now = datetime.datetime.now()
	datestr = now.strftime('%Y%m%d_%H%M%S')
	return datestr

def run_netanal(log,center,span,power):
	datestr = get_datestr()
	temp = adr.get_temp()
		if temp == -2.:
			temp = 0.08

	start_freq = center - span/2.
	stop_freq = center + span/2.
	vna.set_start_freq(start_freq)
	vna.set_stop_freq(stop_freq)
	vna.set_power(power)

	# To get a sense of how long acquiring data actually takes
	starttime = time.time()
	fs,re,im = vna.acquire_data()
	endtime = time.time()
	print "Acquired sample in %1.3f seconds" %(endtime - starttime)

	header = 'start_freq %f\n'%start_freq
	header += 'stop_freq %f\n'%stop_freq
	header += 'if_bw %f\n'%if_bw
	header += 'num_pts %d\n'%num_pts
	header += 'frequency %.1fdBm_re %.1fdBm_im\n'%(power,power)
	data = np.array([fs,re,im]).T
	path = os.path.join(out_folder,datestr+'.dat')
	np.savetxt(path,data,header=header,comments='',fmt='%.7e')

		logentry = 'netanal %s %f %f %f %f'%(datestr,temp,center,span,power)
	print>>log,logentry
	log.flush()

def getargs():
	parser = argparse.ArgumentParser()
	parser.add_argument('--noslew',action='store_true',help="Don't slew magnet")
	parser.add_argument('--skipfirstwait',action='store_true',help="Skip first wait")
	args = parser.parse_args()
	return args

def main():
	args = getargs()
	datestr = get_datestr()
	logfn = 'autonetanal_%s.txt'%datestr
	logfile = open(logfn,'w')
	#pdb.set_trace()
	for row in range(LEDMapper.Nrows):
		for col in range(LEDMapper.Ncols):
			mapper.led_on(row, col)
			time.sleep(1)
			print>>mapper.get_state(), logentry
			for magnet_current in magnet_currents:
				slew_date = get_datestr()
				if not args.noslew:
					adr.slew_magnet(magnet_current)
				logentry = 'slew_magnet %s %f'%(slew_date,magnet_current)
				print>>logfile,logentry
				logfile.flush()
				if not (magnet_current == magnet_currents[0] and args.skipfirstwait):
					time.sleep(magnet_slew_delay)
				for power in powers:
					for cold_center in cold_centers:

						k = kid.DarkKIDAl()
						k.f0 = cold_center * kid.u.megahertz
						k.Qiloss = 100000
						temp = adr.get_temp()
						if temp == -2.:
							temp = 0.08
						k.T = temp *  kid.u.kelvin
						k.alphak = 0.56
						k.L = 10.0e-9 * kid.u.henry
						k.Cc = 0.16e-12 * kid.u.farad
						k.Tc = 1.42 * kid.u.kelvin
						k.fill()
						#center = k.fr.to('megahertz').m
						center = cold_center
						run_netanal(logfile,center,span,power)

			mapper.led_off(row, col)
			print>>mapper.get_state(), logentry

if __name__=='__main__':
	main()

