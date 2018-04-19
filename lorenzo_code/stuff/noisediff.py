#this example show the usage of libUSRP library.
import libUSRP
import numpy as np
import matplotlib.pyplot as pl

#tone used in USRP analog mixer. Integer tuning is NOT used
RF      = 2.65e8

#gain used by USRP TX amplifier
tx_gain = 35
rx_gain = 0

#in case this script and the server are not in the same folder, specify the path:
libUSRP.USRP_filepath = ""

tones_all = libUSRP.USRP_VNA_gettones("Power_Scan04/USRP_VNA_20180417_135114")
tones = np.asarray(tones_all)
N_tones = len(tones)

noise_file = "USRP_Noise_20180417_140540" #.h5

libUSRP.USRP_plotnoise(noise_file,  max_freq=1e3)
