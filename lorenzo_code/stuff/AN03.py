#this example show the usage of libUSRP library.
import libUSRP
import numpy as np
import matplotlib.pyplot as pl
import os

#tone used in USRP analog mixer. Integer tuning is NOT used
RF      = 2.65e8

#frequency interval of the VNA scan
start_f = -15e6
last_f  = +25e6

#density of the VNA scan
n_points= 8000

#number of tones handled at the same time
N_tones = 8

#in case this script and the server are not in the same folder, specify the path:
libUSRP.USRP_filepath = ""

#number of gain steps
gain_steps = 1

#maximum gain to use
max_gain = 35.

#create a folder for the scan
foldername = "Power_Scan"

#create list for filenames
scans = []

j = 0
foldername = foldername + str(j)

goodname = False

while goodname == False:
    try:
        os.mkdir(foldername+ str(j))
        os.chdir(foldername+ str(j))
        foldername = foldername+ str(j)
        goodname = True
    except OSError:
        j+=1


#five seconds to connect to the USRP server
if libUSRP.USRP_connect(5):
    
    ################################
    #   TEST AND EVENTUAL CALIB    #
    ################################

    libUSRP.USRP_filetest(RF, 1e7)
    
    ################################
    #          VNA SCAN            #
    ################################
    
    #a = libUSRP.USRP_calibrate(1.3e8, 0, attenuator = 16, USRP_maxpower = libUSRP.USRP_outpower)
    #print "CALIBRATION: "+str(a)
    for i in range(gain_steps):
        
        #change gain used by USRP TX amplifier 
        tx_gain = int(float(max_gain*i)/float(gain_steps))
        
        print "Scanning with gain: "+str(tx_gain)+" dB"
    
        #perform VNA scan
        filename = libUSRP.USRP_vna(35, n_points, start_f, last_f, tone = RF, ppt = 10, analyze = False)
        
        #analyze data and store in an other group analyzed data (see HDF5)
        libUSRP.USRP_VNA_analysis(filename)
        
        #save file informations
        scans.append(filename)
    
    
    #plot results
    libUSRP.USRP_plotpowerscan(scans, out_filename = None, expectedNpeaks = 7)

    #exit the folder
    os.chdir("..")

    raw_input("Pressa a key to disconnect")
    
    #nicely shut down client receiver thread
    libUSRP.USRP_disconnect()

