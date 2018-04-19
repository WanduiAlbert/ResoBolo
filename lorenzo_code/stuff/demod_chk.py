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

#five seconds to connect to the USRP server
if libUSRP.USRP_connect(5):
    
    ################################
    #   TEST AND EVENTUAL CALIB    #
    ################################
    #raw_input("Press a key to test...")
    #let's test if everything is working
    
    #libUSRP.USRP_filetest(RF, 1e7)

    if(libUSRP.USRP_connected):
    
        #tones = [ 10e6 ]
        tones_all = libUSRP.USRP_VNA_gettones("Power_Scan04/USRP_VNA_20180417_135114")
        tones = np.asarray(tones_all)
        N_tones = len(tones)
        #ampli = [0.1 for i in range(N_tones)]
        #tones = np.insert(tones, len(tones)-1, tones_all[4:5])
        tones = tones - RF #- 10e3
        
        print "TONES: "+str(tones)

        rate_u = 1e8
        decimation = 1000
        acc_time = 120
        welch = 10


        if(decimation>0):
            my_rate = rate_u/decimation
        else:
            my_rate = rate_u
            
        filename_noise = libUSRP.USRP_getnoise(tones, rate_u, acc_time, decimation, tx_gain, rx_gain, RFtone = RF)
        libUSRP.USRP_analyze_noise(filename_noise, welch)
	my_name = filename_noise + ("_powerpt_%2f" % (tx_gain -10*np.log10(N_tones)))
        libUSRP.USRP_plotnoise(filename_noise, out_filename = my_name, max_freq = 1e3)#
        print "File name is: "+ filename_noise





    raw_input("Pressa a key to disconnect")
    #nicely shut down client receiver thread
    libUSRP.USRP_disconnect()
   
