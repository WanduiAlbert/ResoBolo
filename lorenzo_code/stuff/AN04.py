#this example show the usage of libUSRP library.
import libUSRP
import numpy as np
import matplotlib.pyplot as pl
import os

#tone used in USRP analog mixer. Integer tuning is NOT used
RF      = 2.6e8

#number of tones handled at the same time

decimation = 500
rate = 5e7
#in case this script and the server are not in the same folder, specify the path:
libUSRP.USRP_filepath = ""

#five seconds to connect to the USRP server
if libUSRP.USRP_connect(5):
    
    ################################
    #   TEST AND EVENTUAL CALIB    #
    ################################

    #libUSRP.USRP_filetest(RF, 1e8)
    
    ################################
    #       Start live plot        #
    ################################
    
    #tones_all = libUSRP.USRP_VNA_gettones("USRP_VNA_20180410_141602")
    tones = [253.305e6, 254.248e6, 254.514e6, 283.687e6,283.950e6, 285.974e6]#tones_all[1:35]
    #tones = np.insert(tones, len(tones)-1, tones_all[4:5])
    tones = np.asarray(tones) - RF
    N_tones = len(tones)
    timer = 600
    
    
    command = "tx on rx on"
    command += " tx_samples:" + str(int(-1))
    command += " rx_samples:" + str(int(-1))
    command += " rate:"+str(rate)
    command += " tone:"+str(RF)
    command += " rx_gain:"+ str(0)
    command += " tx_gain:"+ str(20*np.log10(N_tones))
    command += " tx_delay:" + str(1)
    command += " rx_delay:" + str(1)
    
    for i in range(N_tones):
        command += " wave type:SINE"
        command += " ampl:"+str(1./float(N_tones))#float(N_tones)
        command += " freq:"+str(tones[i])

    for i in range(N_tones):
        command += " dem type:SINE"
        command += " freq:"+str(tones[i])
        command += " decim:"+str(decimation)
    
    
    libUSRP.USRP_user_typed(command)

    USRP_getnoise(tones, rate, 500, decimation, 6, RFtone = -1, ampl = [1./float(N_tones),])

    #libUSRP.USRP_livedemo("live_demo",rate,decimation,RF , np.asarray(tones),timer, max_samp = None, timeout = None, clear_queue = True,fps_skip = 10 )
    
    raw_input("Pressa a key to disconnect")
    
    #nicely shut down client receiver thread
    libUSRP.USRP_disconnect()

