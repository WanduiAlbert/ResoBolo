###########################################################################
#   THIS LIBRARY IS USEFUL FOR COMMUNIICATING WITH THE USRP SERVER.       #
#   WILL SOON BECOME AN APPLICATION NOTE INCLUDED IN THE DOCUMENTATION    #
###########################################################################
import os
#math libraries
import numpy as np
from scipy import signal#, fftpack
#from scipy.signal import find_peaks_cwt
#from scipy import stats
#from scipy.signal import chirp, sweep_poly
from scipy.optimize import curve_fit
from scipy import optimize
from math import pi
#from scipy import signal, fftpack

#utility libraries
import datetime
import socket
import select
import sys
import time
import threading
from threading import Thread
import signal as ossignal
from functools import partial
import glob
from Queue import Empty
import tempfile
import shutil
from shutil import copyfile
import csv
import gc
from joblib import Parallel, delayed
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import Pipe
from multiprocessing.queues import SimpleQueue
from Queue import *
from struct import *
import array
import h5py
import affinity
import multiprocessing 
import peakutils
from scipy import signal

#plotting libraries
import matplotlib.pyplot as pl
import matplotlib.animation as animation
import matplotlib
#import mpld3
#from mpld3 import plugins
from matplotlib.pyplot import cm 
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import plotly
import colorlover as cl
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import Scatter, Layout
import matplotlib.colors
from matplotlib.colors import colorConverter
from vispy import gloo
from vispy import app
from matplotlib.animation import FuncAnimation
from plotly import tools
from plotly.tools import FigureFactory as FF

#disable core affinity setting
#os.system("taskset -p 0xff %d" % os.getpid())
#affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)

################################
#    GLOBAL VARIABLE DECL      #
################################

#ip address of USRP server
USRP_IP_ADDR = '192.168.5.20'

#soket used for command
USRP_server_address = (USRP_IP_ADDR, 22000)
USRP_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

#address used for data
USRP_server_address_data = (USRP_IP_ADDR, 61360)
USRP_data_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

#queue for fassing data from network
USRP_data_queue = Queue()
#USRP_data_queue = SimpleQueue()

#enable debug msgs
USRP_server_debug   = True

#enable TCP data conection
USRP_data_service = True

#enable TCP data conection state
USRP_data_service_connected = False

#enable TCP command connecion
USRP_recive_service = True

#needed in lat2file fcn for stopping live file writing
USRP_stop_live_filewrite = False

#server is connected
USRP_connected = False

#end of a meas
USRP_endFlag    = False

#start of meas
USRP_startFlag  = False

#end of a burst
USRP_burstFlag  = False

#ready ti interp command
USRP_readyFlag  = False

USRP_filepath = ""

#sinence the empty filepath after one showing
USRP_pathwarning = True

#global calibration constant (Vrms per ADC unit)
USRP_calibration =  5.775091506536823e-06

#line delay
USRP_LineDelay = 0

#the rate associate with the dealy
# in USRP the delay is function of the rate
USRP_LineDelay_Rate = 0

#Maximum rate allowed
USRP_maxrate = 1.e8

#maximum output power of USRP in dBm @ 50 Ohm
USRP_outpower = -6

################################
#            CLASSES           #
################################

class tx_tone:
    '''
    wave_type    = "DC"
    frequency    = 0
    amplitude    = 0
    steps        = 0
    to_frequency = 0
    lapse        = 0
    '''
    def __init__ (self,wave_type_= None,frequency_= None,to_frequency_= None,lapse_= None,steps_= None,amplitude_= None):
        if(wave_type_ != None):
            self.wave_type = wave_type_
        else:
            self.wave_type    = "DC"
        if(frequency_ != None): 
            self.frequency = frequency_
        else:
            self.frequency    = 0
        if(amplitude_ != None):
            self.amplitude = amplitude_
        else:
            self.amplitude    = 0
        if(steps_ != None):
            self.steps = steps_
        else:
            self.steps        = 0
        if(to_frequency_ != None):
            self.to_frequency = to_frequency_
        else:
            self.to_frequency = 0
        if(lapse_ != None):
            self.lapse = lapse_
        else:
            self.lapse        = 0


    
    def set(self,wave_type_,frequency_= None,to_frequency_= None,lapse_= None,steps_= None,amplitude_= None):
        if wave_type_ == "SINE":
            wave_type = wave_type_
            frequency = frequency_
            amplitude = amplitude_
            
        elif wave_type_ == "SWIPE":
            wave_type = wave_type_
            frequency = frequency_
            amplitude = amplitude_
            steps = steps_
            to_frequency = to_frequency_
            lapse = lapse_
            
        elif wave_type_ == "NOISE":
            wave_type = wave_type_
            amplitude = amplitude_
            
        elif wave_type_ == "DC":
            wave_type = wave_type_
            amplitude = amplitude_
            
        elif wave_type_ == "RAMP":
            wave_type = wave_type_
            frequency = frequency_
            amplitude = amplitude_
            
        elif wave_type_ == "SQUARE":
            wave_type = wave_type_
            frequency = frequency_
            amplitude = amplitude_ 
        else:
            print "WARNING: a TX wave type is not well set!"
            
    def to_string(self):
        command = " wave"
        command += " wave_type:"+str(self.wave_type)
        command += " ampl:" + str(self.amplitude)
        
        if self.wave_type != "DC" and self.wave_type != "NOISE":
            command += " freq:" + str(self.frequency)
            
        if self.wave_type == "SWIPE":
            command += " to:" + str(self.to_frequency)
            command += " lapse:" + str(self.lapse)
            command += " steps:" +str( self.steps)
            
        return command
    
class rx_tone:
    '''
    wave_type    = "DC"
    frequency    = 0
    decimation    = 0
    steps        = 0
    to_frequency = 0
    lapse        = 0
    '''
    def __init__ (self,wave_type_= None,frequency_= None,to_frequency_= None,lapse_= None,steps_= None,decimation_= None):
        if(wave_type_ != None):
            self.wave_type = wave_type_
        else:
            self.wave_type    = "DC"
        if(frequency_ != None): 
            self.frequency = frequency_
        else:
            self.frequency    = 0
        if(decimation_ != None):
            self.decimation = decimation_
        else:
            self.decimation    = 0
        if(steps_ != None):
            self.steps = steps_
        else:
            self.steps        = 0
        if(to_frequency_ != None):
            self.to_frequency = to_frequency_
        else:
            self.to_frequency = 0
        if(lapse_ != None):
            self.lapse = lapse_
        else:
            self.lapse        = 0


    
    def set(self, wave_type_,frequency_= None,to_frequency_= None,lapse_= None,steps_= None,decimation_= None):
        if wave_type_ == "SINE":
            wave_type = wave_type_
            frequency = frequency_
            decimation = decimation_
            
        elif wave_type_ == "SWIPE":
            wave_type = wave_type_
            frequency = frequency_
            decimation = decimation_
            steps = steps_
            to_frequency = to_frequency_
            lapse = lapse_
            
        elif wave_type_ == "DC":
            wave_type = wave_type_
            decimation = decimation_
            
        else:
            print "WARNING: a RX wave type is not well set!"

    def to_string(self):
        command = " dem"
        command += " wave_type:"+str(self.wave_type)
        command += " decim:" + str(self.decimation)
        
        if self.wave_type != "DC":
            command += " freq:" + str(self.frequency)
            
        if self.wave_type == "SWIPE":
            command += " to:" + str(self.to_frequency)
            command += " lapse:" + str(self.lapse)
            command += " steps:" + str(self.steps)
            
        return command
            
class parameters:
    '''
    tx = False
    rx = False
    rx_samples = 0
    tx_samples = 0
    rx_delay = 0
    tx_delay = 0
    rx_gain = 0
    tx_gain = 0
    rx_tones = []
    tx_tones = []
    tone = 0
    rate = 0
    '''
    
    def __init__ (self,tx_ = None,rx_ = None,rx_samples_= None,tx_samples_= None, tx_delay_ = None, rx_delay_= None, rate_= None, tx_gain_= None, rx_gain_= None, tone_= None,  rx_tones_= None, tx_tones_= None):
    
        self.rx_tone_num = 0
        self.tx_tone_num = 0
        
        if(tx_ ==None):
            self.tx = False
        else:
            self.tx = tx_
        if(rx_ ==None):    
            self.rx = False
        else:
            self.rx = rx_
        if(rx_samples_ ==None):   
            self.rx_samples = 0
        else:
             self.rx_samples = rx_samples_
        if(tx_samples_ ==None):   
            self.tx_samples = 0
        else:
             self.tx_samples = tx_samples_
        if(rx_delay_ ==None):   
            self.rx_delay = 0
        else:
             self.rx_delay = rx_delay_
        if(tx_delay_ ==None):   
            self.tx_delay = 0
        else:
             self.tx_delay = tx_delay_
        if(rx_gain_ ==None):   
            self.rx_gain = 0
        else:
             self.rx_gain = rx_gain_
        if(tx_gain_ ==None):   
            self.tx_gain = 0
        else:
             self.tx_gain = tx_gain_
        if(rx_tones_ ==None):   
            self.rx_tones = []
        else:
            self.rx_tones = []
            for i in rx_tones_:
                self.rx_tone_num += 1
                self.rx_tones.append(i) 
        if(tx_tones_ ==None):   
            self.tx_tones = []
        else:
            self.tx_tones = []
            for i in tx_tones_:
                self.tx_tone_num += 1
                self.tx_tones.append(i)
        if(rate_ ==None):   
            self.rate = 0
        else:
             self.rate = rate_
        if(tone_ ==None):   
            self.tone = 0
        else:
             self.tone = tone_
             


    
    def add_tx_tone( self,tx_single_tone  ):
        
        if tx_single_tone.amplitude == 0 :
            print "Cannot add tone to parameters array, amplitude is 0!"
        else:
            self.tx = True
            self.rx_tone_num += 1
            self.tx_tones.append(tx_single_tone)
            
    def add_rx_tone(self, rx_single_tone  ):
       
        if rx_single_tone.decimation > 0 :
            self.rx = True
            self.tx_tone_num += 1
            self.rx_tones.append(rx_single_tone)
        else:
            print "Cannot add tone to parameters array, decimation is < 0!"
    
    #def check(self):
        #consistency
        #nyquinst
        
    #def from_string(self, string):
        
    def to_string(self):
        command = ""
        if(self.tx == True):
            command += " tx on"
        if(self.rx == True):
            command += " rx on"
            
        command += " rx_samples:" + str(self.rx_samples)
        command += " tx_samples:" + str(self.tx_samples)
        command += " rx_dealy:"+str(self.rx_delay)
        command += " tx_delay:"+str(self.tx_delay)
        command += " rx_gain:"+str(self.rx_gain)
        command += " tx_gain:"+str(self.tx_gain)
        command += " rate:"+ str(self.rate)
        command += " tone:"+str(self.tone)
        
        for single_tone in self.tx_tones:
             command += single_tone.to_string()

        for single_tone in self.rx_tones:
             command += single_tone.to_string()
        
        return command
             
################################
#    COMMUNICATION FUNCTIONS   #
################################

#push in a queue the rececived samples from USRP server
def USRP_data():
    print 'Starting UDP data service...'
    global USRP_data_socket
    global USRP_data_service
    global USRP_server_debug
    global USRP_data_queue
    global USRP_connected
    global USRP_data_service_connected

    
    USRP_socket.settimeout(1)
    USRP_data_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    
    #asock.bind(USRP_server_address_data)#may give problems.... try writing down the addr manually (windows)
    if(USRP_socket_bind(USRP_data_socket, USRP_server_address_data, 5)):
        
        print 'UDP data service started.'
        USRP_data_service = True
        USRP_data_service_connected = True
        
        #keep track of different number of channels in measures
        old_num_channels = -1
        
        #warn the user for potential data loss (queue full)
        full_warning = True
        
        #will be flexible TODO
        data_type = np.dtype([('re', np.int16), ('im', np.int16)])
        #data_type = np.dtype([('re', np.float32), ('im', np.float32)])
        
        
        header_type = np.dtype([('num_ch', np.int32), ('act_ch', np.int32), ('data_type', np.int32), ('data_length', np.int32)])
        recv = 0
        while(USRP_data_service):
            h_error = False
            num_channels = 1
            actual_channel = 0
            header = []
            #receive 4 different 4 bytes attributes as header
            try:
                header_data = ""
                h_data = 4*4
                while len(header_data) < 4*4 :
                    header_data += USRP_data_socket.recv(min(4*4, 4*4-len(header_data)))
                
            except socket.error as msg:
                if(USRP_server_debug):
                    print("USRP DATA: " + str(msg) + " , " + "Closing connection...")
                if(not USRP_connected):
                    #USRP_data_socket.shutdown(2)
                    USRP_data_socket.close()
                    
                    USRP_data_service = False
                    print "USRP_S DATA conection down"
                    
            try:
                header = np.fromstring(header_data, dtype=header_type, count=1)
            except ValueError as msg:
                print "Server disconnected!, wrong DATA header. disconnecting..."
                print msg
                #USRP_data_service = False
                h_error = True
                #break
                
            if h_error:
                USRP_data_socket.recv(packet_len)
            
            if(USRP_data_service and not h_error):
            
                actual_channel = int(header[0]['act_ch'])
                
                num_channels = int(header[0]['num_ch'])
                #print "Reveived: ch: "+str(actual_channel)+"/"+str(num_channels)
                
                
                # if the number of channel changes, reshape z
                if(old_num_channels!=num_channels):
                    #print "RESHAPING"
                    z = [[] for i in range(num_channels+1)]
                    old_num_channels = num_channels
                    recv = 0
                    
                #print "Receiving: "+str(header[0]['data_length'] * header[0]['data_type'] * 2) + " Bytes"
                buffer_data = ""
                packet_len = (header[0]['data_length']) * header[0]['data_type'] * 2
                #print "data len: " + str(header[0]['data_length']) + "data type: " + str(header[0]['data_type'])
                while packet_len > (len(buffer_data)):
                    try:
                        recv_len = min(10000, (packet_len-len(buffer_data) ) )
                        buffer_data += USRP_data_socket.recv(recv_len)
                    except socket.error as msg:
                        if(USRP_server_debug):
                            print("USRP DATA: " + str(msg) + " , " + "Closing connection...")
                        if(not USRP_connected):
                            #USRP_data_socket.shutdown(2)
                            USRP_data_socket.close()
                            
                            USRP_data_service = False
                            print "USRP_S DATA conection down"
                       
                data = np.fromstring(buffer_data, dtype=data_type, count=header[0]['data_length'])

                
                #will become more flwxible TODO
                z_ch = data.view(np.int16).astype(np.float32).view(np.complex64)
                #z_ch = data.view(np.float32).astype(np.float32).view(np.complex64)
                
                
                
                #print "Buffer size: "+str(len(z_ch))
                #print "actual: " + str(actual_channel-1) + " total "+str(len(z))
                z[actual_channel-1] = z_ch[:]
                
                #when all the channels are collected, put them in a queue
                #!!!I have to know channel order!
                #print "actual: " + str(actual_channel-1) + " len: " + str(len(z[0]))
                if(actual_channel == len(z)):
                    recv+=1
                    #print "buffers gathered: channels:" + str(actual_channel)+"/"+str(len(z)) + " buffn: "+ str(recv)
                    try:
  
                        USRP_data_queue.put((recv,z[:]))

                        #del z
                    except Full:
                        #should never happen as i is infinite...
                        if(full_warning):
                            print "WARNING: data queue is full! losing samples!"
                            full_warning = False
                        continue

    else:
        print "WARNING: USRP client cannot connect to TCP data socket. Maybe the server is streaming samples to an other address?"
        USRP_data_service_connected = False
        
        
#Get data from lan. if clear_queue is true at the end clear the queue so that next meas is ok.
# if used for continuous acquisition is not a good idea to clear_queue
# if length is not given, will acquire data until EORX flag is risen. if it is givn, accquire until ch1 produced AT LEAST length samples
# timeout handle the case in wich the usrp server is down
#Retuns a list of list z[channel][sample]
def USRP_get_data( length = None, clear_queue = True, timeout = None):#USRP_data_queue = USRP_data_queue,
    global USRP_data_service
    global USRP_server_debug
    global USRP_data_queue
    global USRP_connected
    global USRP_endFlag
    global USRP_data_service_connected
    time.sleep(10)
    if(not USRP_connected):
        print "USRP client recv service error: connect USRP_S first!"
        return None
        
    if(length == None):
        length = sys.maxint
        length_opt = True
    else:
        length_opt = False
        
    if(timeout == None):
        timeout = sys.maxint
        
    acc_samp = 0
    
    #if the numer of channels change during a meas. bad.
    num_chs = -1
    
    while( not (USRP_endFlag and length_opt) and acc_samp < length and USRP_data_service):
        #print "CYCLING"
        try:
            recv , z_tmp = USRP_data_queue.get()
            #print "PACKET INDEX: "+str(recv)
            
            #counter internal to Queue class
            USRP_data_queue.task_done()
            
            #may actually differ by channel to channel...   
            acc_samp += len(z_tmp[0])
            #print "Recv: "+str(acc_samp)  + "/" + str(length)
            
            #set the number of channels and init the list
            if len(z_tmp) != num_chs:
                num_chs = len(z_tmp)
                z = [np.empty([0,1])for i in range(num_chs)]
                #print "test "+str(len(z))
            for i in range(num_chs):
                z[i] = np.append(z[i],np.asarray(z_tmp[i]))
                '''
                deb = z[0]
                pl.plot(deb)
                pl.show()
                '''
                
                #print "incrementing: " + str(len(z[0])) +" because i is: "+str(i)
            #pl.plot(z[0])
            #pl.show()
        except Empty:
            time.sleep(0.1)
            wait += 0.1
            
            if(wait>timeout):
                print "timeout for receiving data. Maybe some data went lost."
                acc_samp = sys.maxint
                break
                
            time.sleep(0.1)
            continue
            

    if clear_queue:
        while(not USRP_data_queue.empty()):
            USRP_data_queue.get()
            USRP_data_queue.task_done()
        USRP_data_queue.join()

    print "returning: " + str(len(z[0]))        
    return z
    



#receive messages from the server an rise relative flags.
#This function have to be launched as thread and suppose a non blocking socket.
def USRP_receive(binded_USRP_socket):
    global USRP_server_debug
    global USRP_recive_service
    global USRP_endFlag
    global USRP_startFlag
    global USRP_burstFlag
    global USRP_connected
    global USRP_readyFlag

    
    USRP_recive_service = True
    
    if(not USRP_connected):
        print "USRP client recv service error: connect USRP_S first!"
        return
    
    rx_flag = 0
    tx_flag = 0
    start_rx_flag = 0
    start_tx_flag = 0
    burst_rx_flag = 0
    burst_tx_flag = 0
    while USRP_recive_service:
        data = ""
        while data.find("\n", 0, len(data))==-1 and USRP_recive_service:
            try:
                new_data = binded_USRP_socket.recv(1)
                data += new_data
            except socket.error as msg:
                if(USRP_server_debug):
                    #print msg
                    continue
                if msg.errno == 107:
                    print "Cannot send/recive command from server! disconnecting..."
                    USRP_connected = False
                    USRP_recive_service = False
                
            if(new_data=="" and USRP_recive_service):
                if USRP_recive_service:
                    try:
                        binded_USRP_socket.send("ping\n")
                        control_data = binded_USRP_socket.recv(4)
                        print control_data
                    except socket.error as msg:
                        if(USRP_server_debug):
                            print("Ping response: " + str(msg) + " , " + "Closing connection...")

                        #binded_USRP_socket.shutdown(2)
                        binded_USRP_socket.close()
                        
                        USRP_recive_service = False
                        USRP_connected = False
                        print "USRP_S conection down"
                else:
                    break

        if(USRP_server_debug and USRP_recive_service):
            print "SERVER_DEBUG: " + str(data).rstrip('\n')
        if USRP_recive_service:
            #end of streaming condition
            rx_cond = str(data).rstrip('\n') == str('EORX')
            if rx_cond:
	        	rx_flag = 1
            tx_cond = str(data).rstrip('\n') == str('EOTX')
            if tx_cond:
	        	tx_flag = 1
            slow_cond = rx_flag + tx_flag == 2
            if slow_cond:
                USRP_endFlag = True
                rx_flag=0
                tx_flag=0
             
            #start of streaming condition
            start_rx_cond = str(data).rstrip('\n') == str('BORX')
            if start_rx_cond:
	        	start_rx_flag = 1
            start_tx_cond = str(data).rstrip('\n') == str('BOTX')
            if start_tx_cond:
	        	start_tx_flag = 1
            start_cond = start_rx_flag + start_tx_flag == 2
            if start_cond:
                USRP_startFlag = True
                start_rx_flag=0
                start_tx_flag=0
                
            #burst streaming condition
            burst_rx_cond = str(data).rstrip('\n') == str('EOBRX')
            if burst_rx_cond:
	        	burst_rx_flag = 1
            burst_tx_cond = str(data).rstrip('\n') == str('EOBTX')
            if burst_tx_cond:
	        	burst_tx_flag = 1
            burst_cond = burst_rx_flag + burst_tx_flag == 2
            if burst_cond:
                USRP_burstFlag = True
                burst_rx_flag=0
                burst_tx_flag=0
                
            #server ready condition
            ready_cond = str(data).rstrip('\n') == str('S_READY')
            if(ready_cond):
                USRP_readyFlag = True
                
            
    if USRP_recive_service:        
        binded_USRP_socket.close()        
    USRP_recive_service = False
    USRP_connected = False

#try to bind the socked for a tcp connection with the usrp server
def USRP_socket_bind(USRP_socket, server_address, timeout):
    if timeout < 0:
        print "WARNING: no USRP_S connection established after timeout."
        return False
    else:
        try:
            USRP_socket.connect(server_address)
            return True
        except socket.error as msg:
            print("Socket binding error" + str(msg) + " , " + "Retrying...")
            time.sleep(1)
            timeout = timeout - 1
            return USRP_socket_bind(USRP_socket, server_address, timeout)


#Block until the ready condition has been reached.
def USRP_ready(timeout = sys.maxint):
    global USRP_readyFlag
    global USRP_recive_service
    
    if(not USRP_recive_service):
        print "USRP_S is not connected!"
        return False
        
    if(USRP_readyFlag):
        USRP_readyFlag = False
        return True
    
    waiting_step = 0.001
    total_wait = 0
    while not USRP_readyFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "ERROR: Missing ready signal!"
            USRP_recive_service = False
            exit(-1)
    USRP_readyFlag = False
    return True


#connect to the USRP_S and launch the streaming service
def USRP_connect(timeout = sys.maxint, LAN = True):#, USRP_data_queue = USRP_data_queue
    global USRP_socket
    global USRP_connected
    global USRP_server_address
    global USRP_stop_live_filewrite
    global USRP_data_queue
    USRP_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    
    #USRP_socket.setblocking(0)
    USRP_socket.settimeout(1)
    
    USRP_connected = USRP_socket_bind(USRP_socket, USRP_server_address, timeout)
    if(USRP_connected):
        print "initializing USPR_S connection\n"
        rec_thread = Thread(target = USRP_receive, args = (USRP_socket,))
        rec_thread.start()
        if(LAN):
            rec_data_thread = Thread(target = USRP_data, args = ())
            rec_data_thread.start()
            USRP_stop_live_filewrite = False
            
        USRP_socket.send("CLIENT_KEY\n")
        USRP_ready(10)
        print "Connected to USRP_S."
    return USRP_connected

#connect to the USRP data port only

#disconnect from USRP Server
def USRP_disconnect():
    global USRP_socket
    global USRP_server_debug
    global USRP_recive_service
    global USRP_endFlag
    global USRP_startFlag
    global USRP_burstFlag
    global USRP_connected
    global USRP_data_service
    global USRP_stop_live_filewrite
    if(USRP_endFlag or USRP_startFlag or USRP_burstFlag):
        print "WARNING: disconnecting from USRP_S with oerations suspended. DATA WILL NOT BE SAVED."
    
    if(USRP_connected):
        USRP_recive_service = False
        USRP_data_service = False
        USRP_stop_live_filewrite = True
        time.sleep(1.5)
        try:
            USRP_socket.shutdown(2)
            USRP_socket.close()
            print "USRP_S conection down"
        except socket.error as msg:
            if(USRP_server_debug):
                print("DEBUG: closig error" + str(msg))
            
        if(USRP_connected):
            print "USRP client error: cannot close connection recv service got stuck."
            return
    else:
        print "USRP_S already disconnected"
        
#reset the USRP and the USRP server
def USRP_reset():

    global USRP_connected
    global USRP_recive_service
    if(USRP_connected):
        USRP_disconnect()
        time.sleep(1)
    if(USRP_connected):
        print "USRP client error: cannot close connection recv service got stuck."
        return
    USRP_connect(10)
    
#Reset all the flag befor a new measurement.
def USRP_resetflags():

    global USRP_endFlag
    global USRP_startFlag
    global USRP_burstFlag
    
    USRP_endFlag    = False
    USRP_startFlag  = False
    USRP_burstFlag  = False


#Block until the flag condition has been reached.
def USRP_burstend(timeout= sys.maxint):
    global USRP_burstFlag
    global USRP_recive_service
    
    if(not USRP_recive_service):
        print "USRP_S is not connected!"
        return False
        
    if(USRP_burstFlag):
        USRP_burstFlag = False
        return True
    
    waiting_step = 0.001
    total_wait = 0
    while not USRP_burstFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "WARNING: Missing end of burst signal!"
            return False
            break
    USRP_burstFlag = False
    time.sleep(0.1)
    return True
    
#Block until the flag condition has been reached.
def USRP_streambegin(timeout = sys.maxint):
    global USRP_startFlag
    global USRP_recive_service
    
    if(not USRP_recive_service):
        print "USRP_S is not connected!"
        return False
        
    if(USRP_startFlag):
        USRP_startFlag = False
        return True
    
    waiting_step = 0.001
    total_wait = 0
    while not USRP_startFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "WARNING: Missing begin of stream signal!"
            return False
            break
    USRP_startFlag = False
    return True
    
#Block until the flag condition has been reached.
def USRP_streamend(timeout = sys.maxint):
    global USRP_endFlag
    global USRP_recive_service
    
    if(not USRP_recive_service):
        print "USRP_S is not connected!"
        return False
        
    if(USRP_endFlag):
        USRP_endFlag = False
        return True
    
    waiting_step = 0.001
    total_wait = 0
    while not USRP_endFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "WARNING: Missing begin of stream signal!"
            return False
            break
    USRP_endFlag = False
    time.sleep(0.1)
    return True


#send command and wait for measurement to end. Than shut TX/RX threads in server
#warning does not reset flags!!
def USRP_send_wait(command, timeout = sys.maxint):
    global USRP_connected
    global USRP_socket
    global USRP_server_debug
    global USRP_burstFlag
    global USRP_endFlag

    try:
    
        USRP_socket.send(command)
        
    except socket.error as msg:
    
        if(USRP_server_debug):
            print msg
        if msg.errno == 107:
            print "Cannot send/recive command from server! disconnecting..."
            USRP_connected = False
    

    waiting_step = 0.01
    total_wait = 0
    while not USRP_burstFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "WARNING: Missing end of burst signal!"
            return False
            break
            

    USRP_socket.send("rx off tx off")
    
    waiting_step = 0.01
    total_wait = 0
    while not USRP_endFlag:
        total_wait += waiting_step
        time.sleep(waiting_step)
        if total_wait > timeout:
            print "WARNING: Missing begin of stream signal!"
            return False
            break


#send a command, wait for meas to end but it's not blocking
def USRP_send_nowait(command , timeout = sys.maxint):

    send_wait = Thread(target = USRP_send_wait, args = (command,timeout))
    send_wait.start()

    return

################################
#    UTILITY FUNCTIONS         #
################################

def vrms2dbm(vp):
    return 10.*np.log10(20.*(vp)**2)
    
def dbm2vrms(dbm):
    return np.sqrt((10**(dbm/10.))/20.)
    
def USRP_filelock(filepath):
    locked = None
    file_object = None
    if os.path.exists(filepath):
        try:
            buffer_size = 32
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'r', buffer_size)
            if file_object:
                #print "%s is not locked." % filepath
                locked = False
        except IOError, message:
            print "File is locked (unable to open in read mode). %s." % \
                  message
            locked = True
        finally:
            if file_object:
                file_object.close()
                #print "%s closed." % filepath
    else:
        print "%s not found." % filepath
    return locked
    
def USRP_filewait(timeout, filename):
    global USRP_filepath
    global USRP_pathwarning
    filepath = USRP_filepath + filename
    if(USRP_filepath == "" and USRP_pathwarning):
        print "WARNING: filepath empty. client has to be in the same location where the USRP server saves files."
        USRP_pathwarning = False
    wait_time = 2
    total_waiting = 0
    while not os.path.exists(filepath):
        if(total_waiting>timeout):
            print "WARNING: file is not accessible."
            return False
        print "%s hasn't arrived. Waiting %s seconds." % \
              (filepath, wait_time)
        total_waiting += wait_time
        time.sleep(wait_time)
    while USRP_filelock(filepath):
        if(total_waiting>timeout):
            print "WARNING: file is not accessible."
            return False
        print "%s is currently in use. Waiting %s seconds." % \
              (filepath, wait_time)
        total_waiting += wait_time
        time.sleep(wait_time)
    return True

#returns the filename corresponding to a channel
def USRP_rxfilename(channel):
    global USRP_filepath
    channel_int = int(channel)
    return USRP_filepath+"rx_samples_chan_" + str(channel_int)+".dat"



#returns the whole set of data in the file corresponding to a certain channel
def USRP_loadfile(channel, start_sample = None, lengh = None):
    if(USRP_filewait(2, USRP_rxfilename(channel))):
        if(start_sample == None and lengh == None):
            return np.memmap(USRP_rxfilename(channel),dtype=np.int16).astype(np.float32).view(np.complex64)
        elif(start_sample == None):
            return np.memmap(USRP_rxfilename(channel),dtype=np.int16,shape=(int(lengh)*2)).astype(np.float32).view(np.complex64)
        else:
            return np.memmap(USRP_rxfilename(channel),dtype=np.int16,offset=int(start_sample)*4 ,shape=(int(lengh)*2)).astype(np.float32).view(np.complex64)
    else:
        print "Cannot map the file!"
        return 0

#open h5 file containing mesurement and return a list of list [channel][sample]
#channel specifies a single channel. if specified returns a list [smaples]
def USRP_openH5file(filename, channel= None, start_sample = None, lengh = None):
    
    #the server has to be update to write h5 files... right now it's only binary
    #so this fcn should be used for USRP_lan2file fcn generated
    global USRP_filepath
    
    if(USRP_filewait(2, USRP_filepath + filename +".h5" ) ):

        f = h5py.File(filename+".h5",'r')

        meas = f["raw_data"]
        Ntones = meas.attrs.get("Ntones")
        
        try:
            range(Ntones)
        except TypeError:
            Ntones = 1
            print "WARNING: cannot find attribute Ntones. setting in auto"
        z_ch = [[] for y in range(Ntones)] #not very elegant
        
        print "Opening "+filename+".h5 in read mode"
        
        
        
        for dataset in meas:

            for ch in range(len(meas[dataset])):
                if(lengh != None and start_sample != None):
                    z_ch[ch].extend(meas[dataset][ch][start_sample:start_sample+lengh])
                elif(lengh != None):
                    z_ch[ch].extend(meas[dataset][ch][:lengh])
                elif(start_sample != None):
                    z_ch[ch].extend(meas[dataset][ch][start_sample:])
                else:
                    z_ch[ch].extend(meas[dataset][ch][:])
        
        f.close()
        
        if(channel != None):
            if(channel>Ntones):
                print "the selected file doues not contain that channel!"
                return [(0,0),(0,0)]
            if(lengh != None and start_sample != None):
                return np.asarray(z_ch[channel][start_sample:start_sample+lengh])
            elif(lengh != None):
                return np.asarray(z_ch[channel][:lengh])
            elif(start_sample != None):
                return np.asarray(z_ch[channel][start_sample:])
            else:
                return np.asarray(z_ch[channel][:])
        else:
            for element in z_ch:
                element = np.asarray(element)
            return z_ch

    

#get from the queue all data untill end of streaming, write a file per channel containing data in a folder
#will be upgraded to h5py
#return True if operation has success False otherwise
#the timeout option is for ignoring the end of stream flag otherwise is blocking
#live plot display a realtime plot of the input data. look at console for timing info
def USRP_lan2file(USRP_data_queue, filename, max_samp = None, timeout = None, clear_queue = True ):#USRP_data_queue = USRP_data_queue
    print "Writing LAN samples to file: "+filename+".h5"
    global USRP_filepath
    global USRP_data_service
    global USRP_endFlag
    global USRP_burstFlag
    #global USRP_data_queue
    global USRP_data_service_connected
    global USRP_stop_live_filewrite
    

    
    #determines how many dataset to create in the file by the lenght of the first element in the queue
    #this will be changed with the implementation of the parameters class. params should be attributes of h5file
    firtst_time_rcive = True
    
    if(not (USRP_data_service and USRP_data_service_connected)):
        print "ERROR: Data streaming socket not connected."
        return False
    
    #create h5 file
    h5f = h5py.File(filename+".h5", 'w')
    grp = h5f.create_group("raw_data")
    grp.attrs.create(name = "date", data="current_date_TODO")
    #counts the received samples in channel 0
    acc_samp = 1
    
    #account for waiting time
    wait = 0
    
    #debug counter
    x = 0
    
    #needed for closing the cycle if max_samp is not defined
    end_condition = True
    
    #signal the end of stream condition have been reached once
    flag_raised = False
    
    if timeout == None:
        timeout = sys.maxint
        
    if max_samp == None:
        max_samp = sys.maxint
        
    datachunk = 1000000
    k = 0
    

    
    
    while (end_condition and not USRP_stop_live_filewrite):
        if(not USRP_data_queue.empty()):
            start = time.time()
            recv,z = USRP_data_queue.get()
            #print "PACKET INDEX "+str(recv)
            USRP_data_queue.task_done()
            
            #if first time recive, create datasets
            if(firtst_time_rcive):
                firtst_time_rcive = False
                N_tones = len(z)
                
                buffer_len = len(z[0])
                line = [None for u in range(N_tones)]
                test = [None for u in range(N_tones)]
                
                #z_accum = [[]for t in range(N_tones)]
                z_accum = [np.empty([0,1])for i in range(N_tones)]
                grp.attrs.create(name = "Ntones", data=N_tones)
                j = 0
            


            for i in range(N_tones):
                #z_accum[i].extend(z[i][:])
                z_accum[i] = np.append(z_accum[i],np.asarray(z[i]))

            j +=  len(z[0])
                
            if j > datachunk:
                ds = grp.create_dataset("dataset_"+str(int(k)), data = z_accum )
                ds.attrs.create(name = "Nsamples", data = j)
                #print "dumping " +str(j) +" on file"
                j = 0
                k += 1
                #z_accum = [[]for t in range(N_tones)]

                z_accum = [np.empty([0,1])for i in range(N_tones)]

            acc_samp += len(z[0])
            #print "LAN recv: "+str(acc_samp)+"/"+str(max_samp)

        else:
            wait += 0.05
            time.sleep(0.05)
            if(wait > timeout):
                print "WARNING: Timeout reached."
                return 



        if(USRP_endFlag):
            flag_raised = True
            
        if max_samp == sys.maxint:
            end_condition = not flag_raised or not USRP_data_queue.empty()
            #print "conditions -flag raised: " + str(flag_raised)+" - data empty: "+ str(USRP_data_queue.empty())
        else:
            end_condition = acc_samp < max_samp
            
    #create last dataset
    ds = grp.create_dataset("dataset_"+str(int(k)), data = z_accum )
    ds.attrs.create(name = "Nsamples", data = j)
    print "samples saved on file"
                
    #sUSRP_endFlag = False
    h5f.close()
    
    if clear_queue:
        while(not USRP_data_queue.empty()):
            print "cleaning"
            USRP_data_queue.get()
            USRP_data_queue.task_done()
        USRP_data_queue.join()
    
    return

#launch lan2file as a thread to write samples on file durning acquisition
def USRP_live_lan2file(filename, max_samp = None, timeout = None, wait4finish = False):
    global USRP_data_queue
    #thread = threading.Thread(target=USRP_lan2file(filename, max_samp, timeout,USRP_data_queue))
    #thread.start()
    
    data_go = Thread(target=USRP_lan2file, args = (USRP_data_queue, filename, max_samp, timeout, True))
    data_go.start()
        
    
    if(wait4finish):
        #thread.join()
        data_go.join()
        
    
#cut delay fx in single tone portion of a SWIPE signal
def USRP_cut_delay_fx(data_single_tone, y_percent_thrashold, x_delta):
    ratio = len(data_single_tone)/float(x_delta)
    x_delta = int(x_delta)
    if ratio < 3:
        print "delta X is too big for cutting delay effect!"
        print "sample size is: "+str(len(data_single_tone))+". the size of the delta is "+str(x_delta)
        return (0,len(data_single_tone))
    if ratio - int(ratio) != 0:
        data_single_tone = data_single_tone[x_delta:int((ratio-1)*x_delta)]
        
    original_lenght = len(data_single_tone)
    ratio = int(len(data_single_tone)/float(x_delta))
    expected = np.abs(np.mean(data_single_tone[int(original_lenght/2-x_delta/2):int(original_lenght/2+x_delta/2)]))
    expected += y_percent_thrashold*expected
    good_index = []
    for i in range(0,ratio-1):
        if np.abs(np.mean(data_single_tone[i*x_delta:(i+1)*x_delta])) < expected:
            good_index = np.append(good_index,i*x_delta)
            
    if len(good_index) == 0:
        print "WARNING: delay has not been optimally set. A tone went missing."
        return (0, len(data_single_tone)-1)
        
    start_index = min(good_index)*x_delta
    last_index = max(good_index)*x_delta
    if last_index == start_index:
        print "cannot find good packet to parse with average!"
        print "average: "+str(expected)
        return (0,len(data_single_tone))
    return (int(start_index),int(last_index))
    
################################
#       TEST FUNCTIONS         #
################################

#test if the filesize resulting from a TX/RX loop has expected size
#return true if the test is passed.false otherwise
def USRP_filetest(tone, rate):
    print "Testing USRP_S..."
    global USRP_socket
    global USRP_connected
    
    
    if(not USRP_connected):
        print "USRP_S is not connected."
        return False
        
   
    
    test_command = parameters(tone_ = int(tone), rate_ = int(rate), tx_samples_ = str(rate), rx_samples_ = str(rate))
    test_command.tx = True
    test_command.rx = True
    
    
    command_string = test_command.to_string()
    #this is due to the buffer of the disk....
    #time.sleep(2)
    
    global USRP_data_service_connected
    
    if(USRP_data_service_connected):

        USRP_send_nowait(command_string)
        print "nowait ended"
        USRP_live_lan2file("test", max_samp = rate, wait4finish = True)
        print "lan2file ended"
        filesize = len(USRP_openH5file("test", channel= 0))
        USRP_resetflags()
        if filesize < rate :
            print "TEST: failed. some error occured:"
            print "filesize: "+str(filesize)+"/"+str(rate)
            return False
        else:
            print "TEST: passed"
            return True
        
    else:
        USRP_send_wait(command_string)
        
        filename = USRP_rxfilename(0)
        
        USRP_filewait(2, filename)
        
        filesize = int(os.stat(filename).st_size)
        
	USRP_resetflags()
        if filesize < 4*rate :
            print "TEST: failed. some error occured:"
            print "filesize: "+str(filesize)+"/"+str(4*rate)
            return False
        else:
            print "TEST: passed"
            return True
        
#check if the delay has been setted according to the rate
def USRP_delayset(tone, rate):

    global USRP_LineDelay_Rate
    global USRP_LineDelay
    
    if(USRP_LineDelay_Rate != rate or USRP_LineDelay == 0):
        print "The delay has not been measured with current rate. measuring the delay..."
        USRP_LineDelay_Rate = rate
        USRP_LineDelay = USRP_delay(tone, rate)

#send user typed command
def USRP_user_typed(string = None):
    global USRP_socket
    global USRP_connected
    global USRP_server_debug
    
    if(not USRP_connected):
        print "USRP_S is not connected."
        return False
        
    if string == None:
        string = raw_input("Type command:\n")

    try:
    
        USRP_socket.send(string)
        return True
        
    except socket.error as msg:
    
        if(USRP_server_debug):
            print msg
            
        if msg.errno == 107: ##does not work....?
            print "Cannot send/recive command from server! disconnecting..."
            USRP_connected = False
            
        return False
        
################################
#    PLOTTING FUNCTIONS        #
################################



#live demo of single tone standard readout
def USRP_livedemo(filename, rate, decimation, RFtone, tones, max_time ,max_samp = None, timeout = None, clear_queue = True, fps_skip = None ):

    global USRP_data_queue
    print "Writing LAN samples to file: "+filename+".h5"
    global USRP_filepath
    global USRP_data_service
    global USRP_endFlag
    global USRP_burstFlag
    global USRP_data_service_connected
    global USRP_stop_live_filewrite
    
    #needed for closing the cycle if max_samp is not defined
    end_condition = True
    time_sig = 0
    
    #determines how many dataset to create in the file by the lenght of the first element in the queue
    #this will be changed with the implementation of the parameters class. params should be attributes of h5file
    firtst_time_rcive = True
    
    if(not (USRP_data_service and USRP_data_service_connected)):
        print "ERROR: Data streaming socket not connected."
        return False
    
    #create h5 file
    h5f = h5py.File(filename+".h5", 'w')
    grp = h5f.create_group("raw_data")

    #counts the received samples in channel 0
    acc_samp = 1
    
    #account for waiting time
    wait = 0
    
    #debug counter
    x = 0

    #signal the end of stream condition have been reached once
    flag_raised = False
    
    if timeout == None:
        timeout = sys.maxint
        
    if max_samp == None:
        max_samp = sys.maxint
    
    #to store on file
    datachunk = 50
    k = 0

    #not sure about location of next line
    #signal.pause()

    color_base=iter(cm.rainbow(np.linspace(0,1,10)))
    color_ = []
    for d in range(10):
        color_.append(next(color_base))
    fig = pl.figure(figsize=(12,6) )
    ax1 = fig.add_subplot(1,1,1)
    ax1.cla()
    ax1.set_title("Beam mapping utility")
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel("Delta dBm")
    #pl.ion()
    
    
    my_range = 1e6 #plotting range in samples (pre fps-skip)
    fps_skip = int(1e4) #additional decimation factor(only gaphics)

    usrp_buff = 10e6
    time_factor = usrp_buff/float(rate)
    acc_time = 0
    n_frame = 0
    spp = rate/float(usrp_buff)
    limit_axis = -int(my_range*decimation/1.e6)

    #pass_band = 0.0001
    #stop_band = 0.001
    #c_nom,c_denom = signal.iirdesign(pass_band,stop_band,1,80,ftype="butter")
    re_scale = True
    manual_stop = True
    j = 0
    kk = 0
    while (end_condition and not USRP_stop_live_filewrite and manual_stop):
        if(not USRP_data_queue.empty()):
            start = time.time() 
            recv,z = USRP_data_queue.get()
            USRP_data_queue.task_done()
            
            #if first time recive, create datasets
            if(firtst_time_rcive):
                firtst_time_rcive = False
                N_tones = len(z)
                buffer_len = len(z[0])
                line = [None for u in range(N_tones)]
                test = [None for u in range(N_tones)]
                #z_accum = [[]for t in range(N_tones)]
                j = 0
                jj = 0
                data = [None for i in range(N_tones)]
                off_data = []
                grp.attrs.__setitem__("Ntones", N_tones)
                grp.attrs.create(name = "rate", data=rate)
                grp.attrs.create(name = "decim", data=decimation)
                grp.attrs.create(name = "freq", data=tones)
                grp.attrs.create(name = "rf", data=RFtone)
                grp.attrs.create(name = "tx_gain", data=0)
                grp.attrs.create(name = "meas_type", data="Recorded Streams")
                grp.attrs.create(name = "time_const", data=time_factor)

                x_data = np.asarray(j*time_factor)
                
                z_accum = [np.empty([0,1])for i in range(N_tones)]
                for i in range(N_tones):

                    z[i] = np.abs(z[i])
                    #signal.lfilter(c_nom,c_denom, z[i])
                    tmp = vrms2dbm(USRP_calibration* np.mean(z[i]) )
                    data[i] = tmp
                    z_accum[i] = np.append(z_accum[i],tmp)
                    #z[i] = np.angle(z[i])
                    #data[i] = np.asarray(np.mean(z[i]) )
   
                    cl=color_[i%10]
                    
                    off_data.append(data[i])
                    #off_data.append(0)
                    line[i], = ax1.plot(0, 0, alpha=0.8, color=cl, label =  ("%.1f dBm\n%.3f Mhz" % (data[i],tones[i]/1.e6)))
                
                
                pl.show(block=False)
                print 'initializing normalization mask:\n' + str(off_data)
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig.show()
                xmin, xmax, ymin, ymax = [max(0, (buffer_len+jj)-my_range)*time_factor, max( (buffer_len+jj),my_range)*time_factor , -1000,1000]
                pl.axis([xmin, xmax, ymin, ymax])
                fig.canvas.draw()
        
            else:
                x_data = np.append(x_data, j*time_factor)
                
                for i in range(N_tones):

                    z[i] = np.abs(z[i])
                    #signal.lfilter(c_nom,c_denom, z[i])
                    tmp = (vrms2dbm(USRP_calibration*np.mean( z[i]))) - off_data[i]
                    data[i] = np.append(data[i],  tmp )
                    z_accum[i] = np.append(z_accum[i],tmp)
                    #z[i] = np.angle(z[i])
                    #data[i] = np.append(data[i], (np.mean( z[i]) ) - off_data[i]  )

                    line[i].set_data(x_data,data[i])

                pl.axis([xmin, xmax, ymin, ymax])
                fig.canvas.blit(ax1.bbox)

            j +=  1#buffer_len
            jj +=  buffer_len
        
            acc_samp += len(z[0])
            
            #cut some stuff to keep fps stable
            '''
            if len(x_data) > -limit_axis:
                x_data = x_data[limit_axis:]
                for i in range(N_tones):
                    data[i] = data[i][limit_axis:]
            '''
            time_sig+=buffer_len/1.e6
            '''
            if int(time_sig)%20 == 0 and re_scale:
                print "Rescaling"
                for i in range(N_tones):
                    off_data[i] +=np.mean(data[i])
                re_scale =False

            if int(time_sig)%21 == 1:
                re_scale =True
            '''
            if n_frame%datachunk == (datachunk-1):
                ds = grp.create_dataset("dataset_"+str(int(kk)), data = z_accum )
                ds.attrs.create(name = "Nsamples", data = len(z_accum))
                print "dumping " +str(len(z_accum[0])) +" on file"
                print "single ch: " + str(len(z_accum[0]))+" all ch: " + str(len(z_accum))
                kk +=1
                #z_accum = [[]for t in range(N_tones)]
                
                z_accum = [np.empty([0,1])for i in range(N_tones)]
            
            n_frame +=1
            if(acc_samp%(20*buffer_len) == 1):
                x_data = x_data[limit_axis:]
                xmin, xmax, ymin, ymax = [((buffer_len+jj)-my_range)*decimation/float(rate), ( (buffer_len+jj))*decimation/float(rate) , np.asarray(data).min(),np.asarray(data).max()]
                xmin, xmax, ymin, ymax = [min(x_data), max(x_data) , np.asarray(data).min(),np.asarray(data).max()]
                
                for i in range(N_tones):
                    data[i] = data[i][limit_axis:]
                

                pl.pause(0.001)
                stop = time.time()
                acc_time += (-start+stop)
                print "Cumulated frame time "+str(acc_time)+" RT: "+str(0.01*n_frame)
                acc_time = 0
                n_frame = 0
    
    
        else:
            wait += 0.05
            time.sleep(0.05)
            if(wait > timeout):
                print "WARNING: Timeout reached."
                return 
        
        if(j*time_factor > max_time):
            USRP_user_typed("tx off rx off")
        
        if(USRP_endFlag):
            flag_raised = True
            
        if max_samp == sys.maxint:
            end_condition = not flag_raised or not USRP_data_queue.empty()
            #print "conditions -flag raised: " + str(flag_raised)+" - data empty: "+ str(USRP_data_queue.empty())
        else:
            end_condition = acc_samp < max_samp
            
        

    #create last dataset
    ds = grp.create_dataset("dataset_"+str(int(kk)), data = z_accum )
    ds.attrs.create(name = "Nsamples", data = len(z_accum[0]))

                
    #sUSRP_endFlag = False
    h5f.close()
    print "samples saved on file"
    if clear_queue:
        while(not USRP_data_queue.empty()):
            print "cleaning"
            USRP_data_queue.get()
            USRP_data_queue.task_done()
        USRP_data_queue.join()
    print "quitting live demo..."
    return

def USRP_noise_eval(complex_samples, sample_rate, NPERS_DIV = None):
    print "Running job"
    #os.system("taskset -p ffffffff %d" % os.getpid())
    z = np.asarray(complex_samples)/np.median(np.asarray(complex_samples))
    z -= np.median(z)
    if(NPERS_DIV!= None):
        N = nperseg=len(z)/NPERS_DIV
    else:
        N = len(z)
    Frequencies , RealPart      = signal.welch( z.real ,nperseg=N, fs=sample_rate ,detrend='linear',scaling='density',window = signal.get_window(('kaiser', 4.0), N))
    Frequencies , ImaginaryPart = signal.welch( z.imag ,nperseg=N, fs=sample_rate ,detrend='linear',scaling='density',window = signal.get_window(('kaiser', 4.0), N))
    
    return Frequencies, 10*np.log10(RealPart), 10*np.log10(ImaginaryPart)



#create a group in a noise measurement h5 file containing Im and Re spectra
#use parallel CPU analysis with whelch methods (1 core per channel)
#NPERS_DIV is used in whelch method to define the lenght of the segment to use. it will be lenght/NPERS_DIV
def USRP_analyze_noise(filename, NPERS_DIV = None):

    #open h5 file
    fv = h5py.File(filename+".h5",'r+')
    meas = fv["raw_data"]

    #retrive some attributess
    rate = meas.attrs.get("rate")
    decim = meas.attrs.get("decim")
    power = meas.attrs.get("power")

    freq = meas.attrs.get("freq")
    rf = meas.attrs.get("rf")
    
    
    #effective sample rate
    if(decim!=0):
        sample_rate  = rate/float(decim)
    else:
        sample_rate = rate
    
    Z = USRP_openH5file(filename)
    print "Running parallel FFT on CPU..."
    
    #run parallel analysis: returns when all channels have been analyzed
    Results =Parallel(n_jobs=-1,verbose=1 )(delayed(USRP_noise_eval)(np.asarray(Z[i][:]), sample_rate, NPERS_DIV) for i in range(len(Z)))

    

    
    print "Saving result on file "+filename+" ..."
    
    '''
    print "len result = "+str(len(Results))+" len result [0] = "+str(len(Results[0]))
    for i in range(len(z)):
        pl.plot(Results[i][0],Results[i][1])
        
    pl.xscale('log')

    '''
    spec = fv.create_group("noise")
    spec.attrs.create(name = "rate_eff", data=sample_rate)
    spec.attrs.create(name = "N_tones", data=len(Z))
    
    for i in range(len(Z)):
        spec.create_dataset("real_"+str(i), data = Results[i][1] )
        spec.create_dataset("imag_"+str(i), data = Results[i][2] )
        spec.create_dataset("freq_"+str(i), data = Results[i][0] )
    
    fv.close()
    print "Analysis done."

#plot the results of multiple VNA scan.
#filenames and gains are lists of filename containing pre analyzed VNA scans.
def USRP_plotpowerscan(filenames, out_filename = None, show_tones = True, expectedNpeaks = 125, deltaNpeaks = None):
    
    #output the htlm file
    if out_filename == None:
        out_filename = "PScan_"+filenames[0]
    
    fig = tools.make_subplots(rows=2, cols=1,subplot_titles=('Magnitude', 'Phase'))

    #setup color rotation
    color_base=iter(cm.rainbow(np.linspace(0,1,10)))
    color_ = []
    for d in range(10):
        color_.append(next(color_base))

    for i in range(len(filenames)):
        
        cl=color_[i%10]#next(color_)
        r = colorConverter.to_rgba(cl)
        c = "rgba("+str(r[0])+", "+str(r[1])+", "+str(r[2])+", 0.8)"
        
        fv = h5py.File(filenames[i]+".h5",'r+')
        measr = fv["raw_data"]
        tx_gain = measr.attrs.get("gain")

        #retrive some attributes

        meas = (fv["vna"])
        freq = np.asarray((meas["freq"]))[:]
        S21 = (meas["S21"])
        reso = meas.attrs.get("reso")
        
        phase = np.angle(S21)
        phase = np.unwrap(phase)
        m,q = np.polyfit((np.asarray(freq))/1.e6, phase, 1)
        linear_phase = m*(np.asarray(freq))/1.e6
        phase -= linear_phase
        phase -= np.mean(phase)
        
        magnitude = np.abs(S21)
        magnitudedb = vrms2dbm(magnitude)

        bad_index = []
        for x in range(len(magnitudedb)):
            if magnitudedb[x] < -1000:
                bad_index.append(x)
                print "appending "+str(x)

        print "WARNING: some point has not been correctly measures. Try to increase ppt in scan fcn"
        print "BAD frqs: (have been removed from measure)"
        print "the lenght is \t\t\t"+str(len(bad_index))
        if len(bad_index)>1:
            if max(bad_index)>len(freq):
                b_index_i = 0
            else:
                b_index_i = max(bad_index)+30

            if min(bad_index)<len(freq):
                b_index_f = 0
            else:
                b_index_f = min(bad_index)-30
        else:
            b_index_f = 0
            b_index_i = 0
                        

        godd_index_i = max(2, b_index_i)
        godd_index_f = min(-2, b_index_f)
        print "indexes    "+str(godd_index_i)+"         "+str(godd_index_f)

        freq=freq[godd_index_i:godd_index_f]
        phase=phase[godd_index_i:godd_index_f]
        magnitudedb=magnitudedb[godd_index_i:godd_index_f]

        #look for peaks
        dist = np.gradient(magnitudedb, 1)
        print "threshold: "+str(max(dist))

        if deltaNpeaks == None:
            deltaNpeaks = np.ceil(expectedNpeaks * 0.1)

        threshold = 1
        threshold_delta = 0.005
        indices = []
        while (len(indices) < expectedNpeaks - deltaNpeaks or len(indices) > expectedNpeaks + deltaNpeaks )and threshold > 0:
            indices = peakutils.indexes(dist, thres=max(dist)*threshold, min_dist=20)
            threshold -= threshold_delta
            print str(threshold)+" found "+str(len(indices))+" resonances"
            if threshold < 0:
                indices = []

        #print str(magnitudedb[0:10])
        print "Utility found " +str(len(indices))+ " resonators"
        meas.attrs.__setitem__("tones", [freq[j] for j in indices])
        #here there should be some feature analysis
        print str(indices)
        traceM = go.Scatter(
                            x = freq,
                            y = magnitudedb,
                            name = "TX: "+str(-6+tx_gain)+" dB",
                            legendgroup = "group" +str(i),
                            line = dict(color = c),
                            mode = 'lines'
                            )
        traceP = go.Scatter(
                            x = freq,
                            y = phase,
                            legendgroup = "group" +str(i),
                            showlegend = False,
                            line = dict(color = c),
                            mode = 'lines'
                            )
        traceSM = go.Scatter(
                             x=[freq[j] for j in indices],
                             y=[magnitudedb[j] for j in indices],
                             mode='markers',
                             legendgroup = "group" +str(i),
                             marker=dict(
                                         size=8,
                                         color = c,
                                         symbol='circle',
                                         opacity = 0.5,
                                         ),
                             name='Detected Peaks'
                             )
        traceSP = go.Scatter(
                             x=[freq[j] for j in indices],
                             y=[phase[j] for j in indices],
                             mode='markers',
                             showlegend = False,
                             legendgroup = "group" +str(i),
                             marker=dict(
                                         size=8,
                                         color = c,
                                         opacity = 0.5,
                                         symbol='circle'
                                         ),
                             name='Detected Peaks'
                             
                             )
        fig.append_trace(traceM, 1, 1)
        fig.append_trace(traceP, 2, 1)
        if(show_tones):
            fig.append_trace(traceSM, 1, 1)
            fig.append_trace(traceSP, 2, 1)
        
        fv.close()

    fig['layout'].update( title="Comparsion of VNA scans\nResolution: "+str(reso)+" Hz")
    fig['layout']['xaxis1'].update(title='Frequency [Hz]')
    fig['layout']['xaxis2'].update(title='Frequency [Hz]')

    fig['layout']['yaxis1'].update(title='Magnitude [dB]')
    fig['layout']['yaxis2'].update(title='Scaled Phase [Rad]')

    plotly.offline.plot(fig, filename=out_filename+".html",auto_open=True)

#generate plot of noise difference between all channels and a map of it in htlm/js format.
#ch_list is the list of channel to consider [1, 4, 32, ...]
def USRP_plotnoisediff(filename, ch_list, out_filename = None, max_freq = None):
    #open h5 file
    fv = h5py.File(filename+".h5",'r')
    meas = fv["noise"]
    N_tones = meas.attrs.get("N_tones")
    rate = meas.attrs.get("rate_eff")
    freq = fv["raw_data"].attrs.get("freq")
    rf = fv["raw_data"].attrs.get("rf")
    
    ch_range = len(ch_list)
    
    #plots = []
    color_=iter(cm.rainbow(np.linspace(0,1,ch_range*ch_range+1)))
    #labels = []
    
    #setup figure Immag (multiplot)
    #fig = pl.figure(figsize=(20, 10))
    #ax = fig.add_subplot(111)
    #ax.set_title("Noise correlation matrix. (Rate/decimation) = "+str(rate/1.e6)+" Msps From file: "+filename+"\n frequency intevall= ["+str(low_f)+","+str(hi_f)+"] Hz")
    #ax.set_xlabel('Frequency [Hz]')
    #ax.set_ylabel('PSD [dBc/sqrtHz]')
    #ax.set_xscale('log')
    
    color_base=iter(cm.rainbow(np.linspace(0,1,10)))
    color_ = []
    for d in range(10):
        color_.append(next(color_base))
    
    #pw1 - pw2
    def subtract(pw1, pw2):
        x = np.abs(10**(np.asarray(pw1)/10.)-10**(np.asarray(pw2)/10.))
        return 10*np.log10(x)
    
    dataPanda = []
    
    #loop over datasets and plot differencies
    for i in range(ch_range):
        for j in range(ch_range):
            cl=color_[(i*(j+1))%10]#next(color_)
            r = colorConverter.to_rgba(cl)
            c = "rgba("+str(r[0])+", "+str(r[1])+", "+str(r[2])+", 0.8)"
            
            if max_freq == None : 
                F = np.asarray(meas["freq_"+str(ch_list[i])])
                Ai = np.asarray(meas["imag_"+str(ch_list[i])])
                Bi = np.asarray(meas["imag_"+str(ch_list[j])])
                Ar = np.asarray(meas["real_"+str(ch_list[i])])
                Br = np.asarray(meas["real_"+str(ch_list[j])])
            else:
                max_freq = int(max_freq)
                F = np.asarray(meas["freq_"+str(ch_list[i])])[:max_freq]
                Ai = np.asarray(meas["imag_"+str(ch_list[i])])[:max_freq]
                Bi = np.asarray(meas["imag_"+str(ch_list[j])])[:max_freq]
                Ar = np.asarray(meas["real_"+str(ch_list[i])])[:max_freq]
                Br = np.asarray(meas["real_"+str(ch_list[j])])[:max_freq]

            if(i>j):

                x = subtract(meas["imag_"+str(ch_list[i])], meas["imag_"+str(ch_list[j])])
                traceA = go.Scatter(
                    x = F,
                    y=Ai, 
                    name = "Imag: "+str( ( rf+freq[int(ch_list[i])] ) /1.e6 ),
                    legendgroup = "groupi" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(dash = 'dash',color = c),
                    visible = "legendonly"
                )
                traceB = go.Scatter(
                    x = F,
                    y = Bi, 
                    name = "Imag: "+str( (rf + freq[ int(ch_list[j]) ])/1.e6 ) ,
                    legendgroup = "groupi" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(dash = 'dot',color = c),
                    visible = "legendonly"
                )
                traceX = go.Scatter(
                    x = F,
                    y = x, 
                    name = "Imag diff: "+str( ( rf+freq[int(ch_list[i])] ) /1.e6 )+ " - " +str( (rf + freq[ int(ch_list[j]) ])/1.e6 ) ,
                    legendgroup = "groupi" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(color = c),
                    visible = "legendonly"
                )
                dataPanda.append(traceA)
                dataPanda.append(traceB) 
                dataPanda.append(traceX)
                
            elif(j>i):
            

                x = subtract(meas["real_"+str(ch_list[i])], meas["real_"+str(ch_list[j])])   
                traceA = go.Scatter(
                    x = F,
                    y=Ar, 
                    name = "Real: "+str((rf+freq[int(ch_list[i])])/1.e6 ) ,
                    legendgroup = "groupr" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(dash = 'dash',color = c),
                    visible = "legendonly"
                )
                traceB = go.Scatter(
                    x = F,
                    y = Br, 
                    name = "Real: "+str((rf+freq[int(ch_list[j])])/1.e6 ),
                    legendgroup = "groupr" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(dash = 'dot',color = c),
                    visible = "legendonly"
                )
                traceX = go.Scatter(
                    x = F,
                    y = x, 
                    name = "Real diff: "+str((rf+freq[int(ch_list[j])])/1.e6 )+ " - " +str((rf+freq[int(ch_list[i])])/1.e6 ) ,
                    legendgroup = "groupr" +str(i*ch_range+j),  
                    mode = 'lines',
                    line = dict(color = c),
                    visible = "legendonly"
                )
                dataPanda.append(traceA)
                dataPanda.append(traceB) 
                dataPanda.append(traceX)
        
    #output the htlm file
    if out_filename == None:
        out_filename = "DIFF_"+filename
        
    fv.close()
    
    if max_freq==None :
        my_name = "Noise spectra correlation. (Rate/decimation) = "+str(rate/1.e6)+" Msps From file: "+filename
    else :
        my_name = "Noise spectra correlation. (Rate/decimation) = "+str(rate/1.e6)+" Msps From file: "+filename + " TRUNCATED: "+str(max_freq)+"/"+str(int(max(F)))+" Hz"
    
    fig0 = dict(data=dataPanda, layout=Layout(
        title=my_name,
        xaxis=dict(type='log',autorange=True,title='Frequency [Hz]'),
        yaxis=dict(title='PSD [dBc/sqrtHz]')
        )
    )
    
    plotly.offline.plot(fig0, filename=out_filename+".html",auto_open=True)


#open already analyzed noise file and plot result in HTLM page
#if ch_list is provided, plots only channels in list
#if max_freq is provided plot only until given freq
def USRP_plotnoise(filename, out_filename = None, max_freq = None, ch_list = None):

    #open h5 file
    fv = h5py.File(filename+".h5",'r')
    meas = fv["noise"]
    N_tones = meas.attrs.get("N_tones")
    rate = meas.attrs.get("rate_eff")
    freq = fv["raw_data"].attrs.get("freq")
    rf = fv["raw_data"].attrs.get("rf")

    color_base=iter(cm.rainbow(np.linspace(0,1,10)))
    color_ = []
    for d in range(10):
        color_.append(next(color_base))

    dataPanda = []
    
    if ch_list == None:
        ch_list = range(N_tones)
    
    #loop over datasets (one per channel)
    for i in range(len(ch_list)):
        cl=color_[i%10]#next(color_)
        r = colorConverter.to_rgba(cl)
        c = "rgba("+str(r[0])+", "+str(r[1])+", "+str(r[2])+", 0.8)"
        if(max_freq == None):

            traceR = go.Scatter(
                x = np.asarray(meas["freq_"+str(ch_list[i])]),
                y=np.asarray(meas["real_"+str(ch_list[i])]), 
                name = str( (rf+freq[i])/1.e6 )+ " MHz - Real" ,
                legendgroup = "group" +str(ch_list[i]),  
                mode = 'lines',
                line = dict(dash = 'dot',color = c),
                visible = "legendonly"
                )
            traceI = go.Scatter(
                x = np.asarray(meas["freq_"+str(ch_list[i])]),
                y=np.asarray(meas["imag_"+str(ch_list[i])]), 
                name = str( (rf+freq[ch_list[i]])/1.e6 )+ " MHz - Imag" ,
                legendgroup = "group" +str(ch_list[i]),  
                mode = 'lines',
                line = dict(color = c),
                visible = "legendonly"
                )
        else:
            arr = meas["freq_"+str(i)][:]
            max_f_ind = (np.abs(arr - max_freq)).argmin()
            traceR = go.Scatter(
                x = np.asarray(meas["freq_"+str(ch_list[i])][:max_f_ind]),
                y=np.asarray(meas["real_"+str(ch_list[i])][:max_f_ind]), 
                name = str( (rf+freq[ch_list[i]])/1.e6 )+ " MHz - Real" ,
                legendgroup = "group" +str(ch_list[i]), 
                mode = 'lines',
                line = dict(dash = 'dot',color = c),
                visible = "legendonly"
                )
            traceI = go.Scatter(
                x = np.asarray(meas["freq_"+str(ch_list[i])][:max_f_ind]),
                y=np.asarray(meas["imag_"+str(ch_list[i])][:max_f_ind]), 
                name = str( (rf+freq[i])/1.e6 )+ " MHz - Imag" ,
                legendgroup = "group" +str(ch_list[i]), 
                line = dict(color = c),
                mode = 'lines',
                visible = "legendonly"
                )
            
        #plots.append(single_plot)
        dataPanda.append(traceR)
        dataPanda.append(traceI)  
    
    #output the htlm file
    if out_filename == None:
        out_filename = filename
        

    fig0 = dict(data=dataPanda, layout=Layout(
        title="Noise spectra. (Rate/decimation) = "+str(rate/1.e6)+" Msps From file: "+filename,
        xaxis=dict(type='log',autorange=True,title='Frequency [Hz]'),
        yaxis=dict(title='PSD [dBc/sqrtHz]')
        )
    )
    
    plotly.offline.plot(fig0, filename=out_filename+".html",auto_open=True)
   
    fv.close()
    
#open already analyzed noise files and plot result in HTLM page
#filenames must be a list of filenames containing the same number of channels
#if ch_list is provided, plots only channels in list
#if max_freq is provided plot only until given freq
def USRP_plotnoise_fromfile(filenames, out_filename = None, max_freq = None, ch_list = None):
    j=0
    my_names =""
    powers = []
    my_scale = cl.scales['8']['qual']['Set1']
    dataPanda = []

        
    for filename in filenames:

        #open h5 file
        fv = h5py.File(filename+".h5",'r')
        meas = fv["noise"]
        grp = fv["raw_data"]
        tx_gain= grp.attrs.get("tx_gain")
        powers.append(tx_gain)
        N_tones = meas.attrs.get("N_tones")
        rate = meas.attrs.get("rate_eff")
        freq = fv["raw_data"].attrs.get("freq")
        rf = fv["raw_data"].attrs.get("rf")

        
        
        if ch_list == None:
            ch_list = range(N_tones)
        
        #loop over datasets (one per channel)
        for i in range(len(ch_list)):

            if(max_freq == None):

                traceR = go.Scatter(
                    x = np.asarray(meas["freq_"+str(ch_list[i])]),
                    y=np.asarray(meas["real_"+str(ch_list[i])]), 
                    name = ("%.2f" % ((rf+freq[i])/1.e6) )+ " MHz - Real\nMeas: "+str(int(j))+" Power: "+str(powers[j]) ,
                    legendgroup = "groupr" +str(ch_list[i]),  
                    mode = 'lines',
                    line = dict(dash = 'dot',color = my_scale[j%8]),
                    visible = "legendonly"
                    )
                traceI = go.Scatter(
                    x = np.asarray(meas["freq_"+str(ch_list[i])]),
                    y=np.asarray(meas["imag_"+str(ch_list[i])]), 
                    name = ("%.2f" % ((rf+freq[i])/1.e6) )+ " MHz - Imag\nMeas: "+str(int(j))+" Power: "+str(powers[j]),
                    legendgroup = "groupi" +str(ch_list[i]),  
                    mode = 'lines',
                    line = dict(color = my_scale[j%8]),
                    visible = "legendonly"
                    )
            else:
                arr = meas["freq_"+str(i)][:]
                max_f_ind = (np.abs(arr - max_freq)).argmin()
                traceR = go.Scatter(
                    x = np.asarray(meas["freq_"+str(ch_list[i])][:max_f_ind]),
                    y=np.asarray(meas["real_"+str(ch_list[i])][:max_f_ind]), 
                    name = ("%.2f" % ((rf+freq[i])/1.e6) )+ " MHz - Real\nMeas: "+str(int(j))+" Power: "+str(powers[j]) ,
                    legendgroup = "groupr" +str(ch_list[i]), 
                    mode = 'lines',
                    line = dict(dash = 'dot',color = my_scale[j%8]),
                    visible = "legendonly"
                    )
                traceI = go.Scatter(
                    x = np.asarray(meas["freq_"+str(ch_list[i])][:max_f_ind]),
                    y=np.asarray(meas["imag_"+str(ch_list[i])][:max_f_ind]), 
                    name = ("%.2f" % ((rf+freq[i])/1.e6) )+ " MHz - Imag\nMeas: "+str(int(j))+" Power: "+str(powers[j]) ,
                    legendgroup = "groupi" +str(ch_list[i]), 
                    line = dict(color = my_scale[j%8]),
                    mode = 'lines',
                    visible = "legendonly"
                    )
                
            #plots.append(single_plot)
            dataPanda.append(traceR)
            dataPanda.append(traceI) 
             
        j+=1
        my_names += " "+ filename
    #output the htlm file
    if out_filename == None:
        out_filename = filename
        

    fig0 = dict(data=dataPanda, layout=Layout(
        title="Noise spectra. (Rate/decimation) = "+str(rate/1.e6)+"\nMsps From file: "+my_names,
        xaxis=dict(type='log',autorange=True,title='Frequency [Hz]'),
        yaxis=dict(title='PSD [dBc/sqrtHz]')
        )
    )
    
    plotly.offline.plot(fig0, filename=out_filename+".html",auto_open=True)
   
    fv.close()
       
#plot timestreams from a noise file
#create a html file named out_filename.htlm containing a javascript implementation of the plot
#max_sample truncate the timestream at a certain sample
#ch_list is a list of channels. if not provvided, plots all
def USRP_plotallchannels(filename = "", out_filename = None, max_sample = None, ch_list = None, offline_decimation = None):

    #open h5 file
    fv = h5py.File(filename+".h5",'r')
    meas = fv["raw_data"]

    #retrive some attributes
    rate = meas.attrs.get("rate")

    decim = meas.attrs.get("decim")
    power = meas.attrs.get("power")
    freq = meas.attrs.get("freq")
    rf = meas.attrs.get("rf")
    
    meas_type = meas.attrs.get("meas_type")
    if meas_type == "Recorded Streams":
        time_const = meas.attrs.get("time_const")
    else:
        time_const = 1
    try:
        if rf == None:
            rf = 0
            print "WARNING legend will be wrong! plot is not a noise acquisition"
        if power == None:
            power = 0
            print "WARNING legend will be wrong! plot is not a noise acquisition"
        if freq == None:
            freq = []
            freq.append(0)
            print "WARNING legend will be wrong! plot is not a noise acquisition"

        if decim == None:
            decim = 0
            print "WARNING legend will be wrong! plot is not a noise acquisition"
        if rate == None:
            rate = 0
            print "WARNING legend will be wrong! plot is not a noise acquisition"
    except ValueError:
        rate = 0
        decim = 0
        #freq = []
        power = 0
        rf = 0
        #freq.append(0)

    #retrive samples
    z = USRP_openH5file(filename)
    

    dataPanda = []
    color_base=iter(cm.rainbow(np.linspace(0,1,10)))
    color_ = []
    for d in range(10):
        color_.append(next(color_base))
        
    if ch_list == None:
        ch_list = range(len(z))
    
    for i in ch_list:
        
        if offline_decimation != None:
            pass_band = max(1./(offline_decimation)-0.2,0.2)
            stop_band = pass_band+0.1
            c_nom,c_denom = signal.iirdesign(pass_band,stop_band,1,80,ftype="butter")
            z[i] = signal.lfilter(c_nom,c_denom, z[i])
            z[i] = signal.decimate(z[i] ,offline_decimation)
        
        cl=color_[i%10]
        r = colorConverter.to_rgba(cl)
        c = "rgba("+str(r[0])+", "+str(r[1])+", "+str(r[2])+", 0.8)"
        
        if max_sample == None:
            I = np.asarray(z[i]).imag
            Q = np.asarray(z[i]).real
            F = np.asarray(range(len(z[i])))*time_const

        else:
            max_sample = int(max_sample)
            
            I = np.asarray(z[i])[:max_sample].imag
            Q = np.asarray(z[i])[:max_sample].real
            F = np.asarray(range(len(z[i])))[:max_sample]*time_const


        if meas_type != "Recorded Streams" and meas_type != None:
            print "meas_type != Recorded Streams: it is: " + meas_type
            traceI = go.Scatter(
                    x = F,
                    y = I,
                    name = str( (rf+freq[i])/1.e6 )+ " MHz - I" ,
                    legendgroup = "group" +str(i),
                    mode = 'lines',
                    line = dict(color = c),
                    visible = "legendonly"
                    )
        traceQ = go.Scatter(
                x = F,
                y = Q, 
                name = str( (rf+freq[i])/1.e6 )+ " MHz - Q" ,
                legendgroup = "group" +str(i),  
                mode = 'lines',
                line = dict(color = c),
                visible = "legendonly"
        )
        
        dataPanda.append(traceQ)
        if meas_type != "Recorded Streams" and meas_type != None:
            dataPanda.append(traceI)


    if out_filename == None:
        out_filename = filename

    if meas_type != None:
        my_title = meas_type
    else:
        my_title = "Raw samples. "
    
    if max_sample == None:
        my_title += "rate: "+str(rate/1e6)+" Msps decimation: "+str(decim)+" \n"+("RF: %.3f" % (rf/1.e6))+" MHz from file: "+filename
    else:
        my_title += "rate: "+str(rate/1e6)+" Msps decimation: "+str(decim)+" \n"+("RF: %.3f" % (rf/1.e6))+" MHz from file: "+filename + " TRUNCATED: "+str(max_sample)+"/"+str(len(z[0]))
    if offline_decimation != None:
        my_title+="Offline decimation applied: "+str(int(offline_decimation))
    fig0 = dict(data=dataPanda, layout=Layout(
        title=my_title,
        xaxis=dict(title="sample #"),
        yaxis=dict(title="ADC units")
        )
    )
    if meas_type == "Recorded Streams":
        fig0 = dict(data=dataPanda, layout=Layout(
          title=my_title,
          xaxis=dict(title="Time [s]"),
          yaxis=dict(title="db relative to initial")
          )
        )
        
    plotly.offline.plot(fig0, filename=out_filename+".html",auto_open=True)  
    #mpld3.save_html(fig, out_filename+"htlm")
    #mpld3.show()

#matplotlib version of previews function
def USRP_plotchannel(ch, rate = None, filename = ""):
    global USRP_data_service_connected
    
    if(rate==None):
        print "Plotting channel without any setted rate. by default it will be 1 MHz"
        rate = 1e6
        
    
    ch = int(ch)
    if USRP_data_service_connected:
        if filename == "":
            print "if data streaming is on, filename must be specified!"
            return 
        z = USRP_openH5file(filename, channel = ch)
        
    else:
        z = USRP_loadfile(ch)


    IsSingleChannel = False
    try:
        len(z[0])
    except TypeError:
        IsSingleChannel = True

    if(IsSingleChannel):

        pfig = pl.figure()
        pl.subplot(211)
        Pxx,freqs = pl.psd(z,Fs= int(rate), NFFT=10000,linestyle='solid',color='b',label="channel"+str(ch))
        pl.legend(loc='lower left')

        ax = pl.subplot(212)
        pl.plot(z.real,color='b')
        pl.ylabel('Real part [ADC units]')
        pl.xlabel('Sample number')
        mean = float(np.mean(z.real))
        Mm = float(float(np.max(z.real)) - float(np.min(z.real)))

        textstr = '$\mathrm{Mean}=%.2f$\n$\mathrm{max-min}=%.2f$\n$\mathrm{rate[MHz]:}%.2f$'%(mean,Mm,1e-6*float(rate))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,verticalalignment='top', bbox=props)
        pl.show(block=False)
        
    else:
        print "Cannot plot multiple channels with this function or wrong data type."

#export pre-analyzed vna hdf5 file in csv with same name  
#filenames is a list of filenames      
def USRP_vna2csv(filenames, out_filename = None, timeanalysis = False):
    for filename in filenames:
        fv = h5py.File(filename+".h5",'r')
        
        measr = fv["raw_data"]
        tx_gain = measr.attrs.get("gain")
        n_points = measr.attrs.get("ppt")
        ppt = measr.attrs.get("n_points")           
        grp.attrs.create(name = "decim", data=decimation)

        meas = (fv["vna"])
        freq = np.asarray((meas["freq"]))[:]
        S21 = (meas["S21"])
        
        reso = meas.attrs.get("reso")
        

        phase = np.angle(S21)
        phase = np.unwrap(phase)
        m,q = np.polyfit((np.asarray(freq))/1.e6, phase, 1)
        linear_phase = m*(np.asarray(freq))/1.e6
        phase -= linear_phase
        phase -= np.mean(phase)
        
        magnitude = np.abs(S21)
        magnitudedb = vrms2dbm(magnitude)

        with open("USRPCSV_"+filename+".csv",'wb') as resultFile:
            wr = csv.writer(resultFile)
            x =("POWER:"+str(tx_gain))
            if timeanalysis:
                wr.writerow(["scantime:"+str,None,None]) 
            wr.writerow([x,None,None])
            wr.writerows(zip(freq,magnitudedb,phase))

#export pre-analyzed vna hdf5 file in csv with same name  
#filenames is a list of filenames
def USRP_samples2csv(filenames):
    for filename in filenames:
        fv = h5py.File(filename+".h5",'r')
        print "extracting "+filename
        measr = fv["raw_data"]
        #tx_gain = measr.attrs.get("gain")
        
        z_ch = USRP_openH5file(filename = filename, channel= None, start_sample = None, lengh = None)
        i = 0
        os.mkdir("CSV_"+filename)
        os.chdir("CSV_"+filename)
        for ch in z_ch:
            print "found "+str(len(ch))+" samples in channel "+str(i)
            I = np.asarray(ch).real
            Q = np.asarray(ch).imag
            with open("USRPCSV_ch"+str(int(i))+"_"+filename+".csv",'wb') as resultFile:
                wr = csv.writer(resultFile)
                wr.writerows(zip(I,Q))
            i+=1
        os.chdir("..")
    print "Export complete!"
            
#plot vna data from a pre-analyzed file, if plot = False return the data (magnitude in db and unwrapped phase) instead of plotting
def USRP_VNA_plot(filename, plot = True, out_filename = None, raw = False):

    fv = h5py.File(filename+".h5",'r')
    measr = fv["raw_data"]
    tx_gain = measr.attrs.get("gain")
    meas = (fv["vna"])
    freq = np.asarray((meas["freq"]))[:]
    S21 = (meas["S21"])
    
    reso = meas.attrs.get("reso")
    

    phase = np.angle(S21)
    magnitude = np.abs(S21)
    if not raw:
        phase = np.unwrap(phase)
        m,q = np.polyfit((np.asarray(freq))/1.e6, phase, 1)
        linear_phase = m*(np.asarray(freq))/1.e6
        phase -= linear_phase
        phase -= np.mean(phase)
        magnitude = vrms2dbm(magnitude)

    
    if(plot):

    
        fig = tools.make_subplots(rows=2, cols=1,subplot_titles=('Magnitude', 'Phase'))
        fig['layout'].update( title='VNA scan from file '+filename+"\nResolution: "+str(reso)+" Hz")
        traceM = go.Scatter(
                x = freq,
                y = magnitude, 
                name = "Magnitude",
                mode = 'lines'
                )
        traceP = go.Scatter(
                x = freq,
                y = phase, 
                name = "Phase",
                mode = 'lines'
                )
                
        fig.append_trace(traceM, 1, 1)
        fig.append_trace(traceP, 2, 1)
        
        fig['layout']['xaxis1'].update(title='Frequency [Hz]')
        fig['layout']['xaxis2'].update(title='Frequency [Hz]')
        
        fig['layout']['yaxis1'].update(title='Magnitude [dB]')
        fig['layout']['yaxis2'].update(title='Scaled Phase [Rad]')

        if out_filename == None:
            out_filename = "VNA_"+filename
            
        plotly.offline.plot(fig, filename=out_filename+".html",auto_open=True)

    else:
        if(raw):
            return tx_gain, freq, np.asarray(S21).real,np.asarray(S21).imag
        else:
            return freq, magnitude, phase

#open already analized (with USRP_VNA_plot_peak) file and output a list containing tone frequencies
def USRP_VNA_gettones(filename):
    fv = h5py.File(filename+".h5",'r')
    meas = (fv["vna"])
    tones = meas.attrs.get("tones")
    return np.asarray(tones)

#take in a preanalized file, fit the resonance in it, returns lapse power f0 qi qr
#fitting stuff has been directli ported from Bryan's code reso_fit
def USRP_vna_fit(filename, p0=None):

    def real_of_complex(z):
	    ''' flatten n-dim complex vector to 2n-dim real vector for fitting '''
	    r = np.hstack((z.real,z.imag))
	    return r
    def complex_of_real(r):
	    assert len(r.shape) == 1
	    nt = r.size
	    assert nt % 2 == 0
	    no = nt/2
	    z = r[:no] + 1j*r[no:]
	    return z

    def model_python(f,f0,A,B,D,Qr,Qe_re,Qe_im):
	    f0 = f0 * 1e6
	    cable_z = np.exp(2.j*pi*(1e-6*D*(f-f0)))
	    Qe = Qe_re + 1.j*Qe_im
	    x = (f - f0)/f0
	    s21 = (A+1.0j*B)*cable_z*(1. - (Qr/Qe)/(1. + 2.j*Qr*x))
	    return real_of_complex(s21)


    model = model_python

    def do_fit_quick(freq,re,im):
	    nt = len(freq)

	    D = -0.0087
	    Qe_re = 31500.
	    Qe_im = 12011.
	    def quick_model(f,f0,A,B,Qr):
		    return model(f,f0,A,B,D,Qr,Qe_re,Qe_im)

	    mag = np.sqrt(re*re+im*im)
	    phase = np.unwrap(np.arctan2(im,re))

	    p = np.polyfit(freq,phase,1)
	    phase -= np.polyval(p,freq)

	    z = mag*np.exp(1.j*phase)
	    re,im = z.real,z.imag

	    f0 = freq[np.argmin(mag)]*1e-6
	    scale = np.max(mag)
	    phi = 0.0
	    A = scale*np.cos(phi)
	    B = scale*np.sin(phi)
	    Qr = 12000
	    p0 = (f0,A,B,Qr)

	    ydata = np.hstack((re,im))

	    popt,pcov = optimize.curve_fit(quick_model,freq,ydata,p0=p0)
	    f0,A,B,Qr = popt

	    yfit = quick_model(freq,*popt)
	    zfit = complex_of_real(yfit)

	    Qi = 1.0/ (1./Qr - 1./Qe_re)

	    return f0,Qr,A,B,zfit

    def do_fit(freq,re,im,p0=None):
	    nt = len(freq)
	    
	    mag = np.sqrt(re*re+im*im)
	    phase = np.unwrap(np.arctan2(im,re))

	    #p = np.polyfit(freq,phase,1)
	    #phase -= np.polyval(p,freq)

	    #z = mag*np.exp(1.j*phase)
	    #re,im = z.real,z.imag
	    
	    if p0 is None:
		    f0 = freq[np.argmin(mag[20:-20])]/1.e6
		    print "f0 start is: "+str(freq[np.argmin(mag[20:-20])]/1.e6)
		    scale = np.max(mag)
		    phi = 1.0
		    A = scale*np.cos(phi)
		    B = scale*np.sin(phi)
		    D = 20e-9
		    Qr=10000
		    Qe_re = 11000
		    Qe_im =2000
		    p0 = (f0,A,B,D,Qr,Qe_re,Qe_im)

	    ydata = np.hstack((re,im))

	    bad_flag = False
	    try:
	        popt,pcov = optimize.curve_fit(model,freq,ydata,p0=p0)
	    except RuntimeError:
	        bad_flag = True
	    
	    if(bad_flag):
	        print "ERROR: fit function failed!"
	        pl.clf()
	        pl.subplot(211)
	        pl.plot(freq,re*re+im*im,label = "magnitude")

	        pl.legend()
	        
	        pl.subplot(212)
	        pl.plot(freq,np.unwrap(np.arctan2(im,re)), label = "phase")

	        pl.legend()
	        pl.show(block=False)
	        return 0,0,0
	    
	    f0,A,B,D,Qr,Qe_re,Qe_im = popt
	    yfit = model(freq,*popt)
	    zfit = complex_of_real(yfit)
	    #print p0
	    #print popt
	    #exit()

	    
	    pl.clf()
	    pl.subplot(211)
	    pl.plot(freq,re*re+im*im,label = "magnitude")
	    pl.plot(freq,np.abs(zfit)**2, label = "mag fit")
	    pl.legend()
	    
	    pl.subplot(212)
	    pl.plot(freq,np.unwrap(np.arctan2(im,re)), label = "phase")
	    pl.plot(freq,np.unwrap(np.angle(zfit)), label = "ph fit")
	    pl.legend()

	    #pl.clf()
	    zm = re + 1.j*im
	    resid = zfit - zm
	    #pl.plot(freq,resid)
	    #pl.hist(resid,bins=30)
	    pl.show(block=False)
	    #exit()
	    #pl.savefig('qpower_%dMHz.png'%(f0*1e6))

	    Qi = 1.0/ (1./Qr - 1./Qe_re)

	    #return f0,Qr,A,B,Qe_re,Qe_im,D,zfit,popt
	    return f0,Qi,Qr

    fv = h5py.File(filename+".h5",'r')

    measr = fv["raw_data"]

    lapse = measr.attrs.get("lapse")
    ppt = measr.attrs.get("n_points")           


    meas = (fv["vna"])
    freq = np.asarray((meas["freq"]))[:]
    S21 = (meas["S21"])

    reso = meas.attrs.get("reso")

    pwr, freq, re, im = USRP_VNA_plot(filename, plot = False, out_filename = None, raw = True)
    f0,Qi,Qr = do_fit(freq[20:-20],re[20:-20],im[20:-20],p0=p0)
    return lapse, pwr,f0,Qi,Qr

    
#take the second derivative of phase and amplitude and plot it
#plot vna data from a pre-analyzed file, return the data (magnitude in db and unwrapped phase) instead of plotting
def USRP_VNA_plot_peak(filename, plot = True, out_filename = None, expectedNpeaks = 115, deltaNpeaks = None):
    
    fv = h5py.File(filename+".h5",'r+')
    
    meas = (fv["vna"])
    freq = np.asarray((meas["freq"]))[:]
    S21 = (meas["S21"])
    
    reso = meas.attrs.get("reso")
    
    
    phase = np.angle(S21)
    phase = np.unwrap(phase)
    m,q = np.polyfit((np.asarray(freq))/1.e6, phase, 1)
    linear_phase = m*(np.asarray(freq))/1.e6
    phase -= linear_phase
    phase -= np.mean(phase)
    
    magnitude = np.abs(S21)
    magnitudedb = vrms2dbm(magnitude)
    
    #check for bad points and remove them
    bad_index = []
    for x in range(len(magnitudedb)):
        if magnitudedb[x] < -1000:
            bad_index.append(x)

    print "WARNING: some point has not been correctly measures. Try to increase ppt in scan fcn"
    print "BAD frqs: (have been removed from measure)"
    #for q in bad_index:
    #    print str(freq[q])+" Hz"

    if len(bad_index) > 0:
        if max(bad_index)>len(freq):
            b_index_i = 0
        else:
            b_index_i = max(bad_index)+30

        if min(bad_index)<len(freq):
            b_index_f = 0
        else:
            b_index_f = min(bad_index)-30
    else:
        b_index_f = 0
        b_index_i = 0

    godd_index_i = max(2, b_index_i)
    godd_index_f = min(-2, b_index_f)
    print "indexes    "+str(godd_index_i)+"         "+str(godd_index_f)

    freq=freq[godd_index_i:godd_index_f]
    phase=phase[godd_index_i:godd_index_f]
    magnitudedb=magnitudedb[godd_index_i:godd_index_f]
    
    #look for peaks
    dist = np.gradient(magnitudedb, 1)
    print "threshold: "+str(max(dist))

    if deltaNpeaks == None:
        deltaNpeaks = np.ceil(expectedNpeaks * 0.1)

    threshold = 1
    threshold_delta = 0.005
    indices = []
    while (len(indices) < expectedNpeaks - deltaNpeaks or len(indices) > expectedNpeaks + deltaNpeaks )and threshold > 0:
        indices = peakutils.indexes(dist, thres=max(dist)*threshold, min_dist=20)
        threshold -= threshold_delta
        #print str(threshold)+" found "+str(len(indices))+" resonances"
        if threshold < 0:
            indices = []

    #print str(magnitudedb[0:10])
    print "Utility found " +str(len(indices))+ " resonators"
    
    meas.attrs.__setitem__("tones", [freq[j] for j in indices])
    
    
    if(plot):
        fig = tools.make_subplots(rows=2, cols=1,subplot_titles=('Magnitude', 'Phase'))
        fig['layout'].update( title='VNA scan from file '+filename+"\nResolution: "+str(reso)+" Hz")
        traceM = go.Scatter(
                            x = freq,
                            y = magnitudedb,
                            name = "Magnitude grad",
                            mode = 'lines'
                            )
        traceP = go.Scatter(
                            x = freq,
                            y = phase,
                            name = "Phase grad",
                            mode = 'lines'
                            )
        traceSM = go.Scatter(
                            x=[freq[j] for j in indices],
                            y=[magnitudedb[j] for j in indices],
                            mode='markers',
                            marker=dict(
                                        size=8,
                                        color='rgb(255,0,0)',
                                        symbol='circle',
                                        opacity = 0.5,
                                        ),
                            name='Detected Peaks'
                            )
        traceSP = go.Scatter(
                             x=[freq[j] for j in indices],
                             y=[phase[j] for j in indices],
                             mode='markers',
                             showlegend = False,
                             marker=dict(
                                         size=8,
                                         color='rgb(255,0,0)',
                                         symbol='circle',
                                         opacity = 0.5,
                                         ),
                             name='Detected Peaks'
                             
                             )
        fig.append_trace(traceM, 1, 1)
        fig.append_trace(traceP, 2, 1)
        fig.append_trace(traceSM, 1, 1)
        fig.append_trace(traceSP, 2, 1)
        
        fig['layout']['xaxis1'].update(title='Frequency [Hz]')
        fig['layout']['xaxis2'].update(title='Frequency [Hz]')
        
        fig['layout']['yaxis1'].update(title='Magnitude [dB]')
        fig['layout']['yaxis2'].update(title='Scaled Phase [Rad]')
        
        if out_filename == None:
            out_filename = "GRAD_"+filename
                            
        plotly.offline.plot(fig, filename=out_filename+".html",auto_open=True)


    return [freq[j] for j in indices]

################################
#    MEASURE FUNCTIONS         #
################################

#returns the line delay of the TX/RX loop
def USRP_delay(tone, rate):
    print "Measuring TX/RX loop delay..."
    global USRP_socket
    global USRP_connected
    
    if(not USRP_connected):
        print "USRP_S is not connected."
        return 0
        
    
    
    tone = int(tone)
    rate = int(rate)
    
    sps = rate*5.
    lapse = 5.
    start_f = -int(rate/2-100)
    end_f = int(rate/2-100)
    coeff = lapse/(end_f-start_f)
    decimate = int(rate/1.e4)
    
    command = "tx on rx on"
    command += " tx_samples:" + str(sps)
    command += " rx_samples:" + str(sps)
    command += " rate:"+str(rate)
    command += " tone:"+str(tone)
    command += " rx_gain:"+ str(0)
    command += " tx_gain:"+ str(0)
    command += " tx_delay:" + str(1)
    command += " rx_delay:" + str(1)
    command += " bw:" + str(rate)

    command += " wave type:SWIPE"
    command += " ampl:"+str(1)
    command += " freq:"+str(start_f)
    command += " to:"+str(end_f)
    command += " lapse:"+str(lapse)
    command += " steps:"+str(1e20)

    command += " dem type:SWIPE"
    command += " decim:"+str(decimate)
    command += " freq:"+str(start_f)
    command += " to:"+str(end_f)
    command += " lapse:"+str(lapse)
    command += " steps:"+str(1e20)


    
    global USRP_data_service_connected
    
    if(USRP_data_service_connected):

        USRP_send_nowait(command)
        
        filename = "USRP_Delay_"+str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
        
        USRP_live_lan2file(filename,  wait4finish = True, max_samp = sps/float(decimate)-1) #max_samp = sps, timeout = 50,
        
        print "Samples saved. Opening..."
        
        z_data = USRP_openH5file(filename, channel= 0) 
        
        f = h5py.File(filename+".h5",'r+')
        grp = f["raw_data"]
        grp.attrs.create(name = "rate", data=rate)
        grp.attrs.create(name = "tone", data=tone)
        f.close()
        
        USRP_resetflags()
        time.sleep(1)
        print "resetting flags"
        
    else:
    
        USRP_send_wait(command)
    
        #this is due to the buffer of the disk....
        time.sleep(2)
        
        z_data = USRP_loadfile(0)
    
    freq,Pxx = signal.welch(z_data.real,nperseg=len(z_data),fs=int(float(rate)/decimate),detrend='linear',scaling='density')
    
    Pxx = np.asarray(Pxx)
    delay = freq[Pxx.argmax()] * coeff
    print "delay at rate "+str(float(rate)/1.e6)+" MHz is "+str(float(delay)*float(1e9))+" ns"
    return delay

#vna scan in a frequency range. saves on HDF5 file
#if analyze = True, returns a frequency axis and the real and immaginary part of the corresponding result as numpy array
#if analyze = False, return the filename associated with the measurement
#NOTE: it does not write the analysis on the HDF5 file. An other function in the library does it.
def USRP_vna(tx_gain, n_points, start_f, last_f, ppt = None, tone = None, analyze = True):
    print "VNA started"
    if ppt == None:
        print "WARNING: Autoset Point Per Tone to 100"
        ppt = 100
    global USRP_LineDelay
    global USRP_connected
    global USRP_LineDelay_Rate
    global USRP_socket
    
   
    
    if(not USRP_connected):
        print "ERROR: USRP_S is not connected!"
        return [(0,0),]
    
    n_points = int(n_points)

    delta_f = max(start_f,last_f)-min(start_f,last_f)
    
    if(tone == None):
        tone = int((min(start_f,last_f) + (delta_f)/2)/1.25e6)*1.25e6
        print "RF tone is: "+str(tone)+" Hz"
        
        start_f -= tone
        last_f -= tone 
    
    maxf = max(np.abs(start_f),np.abs(last_f))
    
    if maxf < 5e6:
        rate = 1e7
        decimation = 100
    elif maxf < 1e7:
        rate = 2e7
        decimation = 200
    elif maxf < 2e7:
        rate = 5e7
        decimation = 500
    elif maxf < 5e7:
        rate = 1e8
        decimation = 1000
    else:
        print "ERROR: Desired frequency interval bigger than bandwidth (100 MHz)!"
        return [(0,0),]
        
    print "Rate is: "+str(rate)+" MHz"
    
    USRP_delayset(tone, rate)
    
    lapse = float(n_points*ppt * decimation) / float(rate)
    print "meas time is: "+str(lapse)+" s"

    #before decimation
    expected_samples = lapse * rate
    
    command = "tx on rx on"
    command += " tx_samples:" + str(expected_samples)
    command += " rx_samples:" + str(expected_samples)
    command += " rate:"+str(rate)
    command += " tone:"+str(tone)
    command += " rx_gain:"+ str(0)
    command += " tx_gain:"+ str(tx_gain)
    command += " tx_delay:" + str(1-USRP_LineDelay)
    command += " rx_delay:" + str(1)
    command += " bw:" + str(max(start_f,last_f)+20e6)

    command += " wave type:SWIPE"
    command += " ampl:"+str(1)
    command += " freq:"+str(start_f)
    command += " to:"+str(last_f)
    command += " lapse:"+str(lapse)
    command += " steps:"+str(int(n_points))

    command += " dem type:SWIPE"
    command += " freq:"+str(start_f)
    command += " to:"+str(last_f)
    command += " lapse:"+str(lapse)
    command += " steps:"+str(n_points)
    command += " decim:"+str(decimation)


    global USRP_data_service_connected
    
    if(USRP_data_service_connected):

        USRP_send_nowait(command)
        
        filename = "USRP_VNA_"+str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
        
        USRP_live_lan2file(filename, wait4finish = True,max_samp = expected_samples/decimation-1) #max_samp = expected_samples/decimation
        
        USRP_resetflags()

        f = h5py.File(filename+".h5",'r+')
        grp = f["raw_data"]
        grp.attrs.create(name = "rate", data=rate)
        grp.attrs.create(name = "lapse", data=lapse)
        grp.attrs.create(name = "decim", data=decimation)
        grp.attrs.create(name = "start_f", data=tone+start_f)
        grp.attrs.create(name = "last_f", data=tone+last_f)
        grp.attrs.create(name = "n_points", data=n_points-1)
        grp.attrs.create(name = "ppt", data=ppt)
        grp.attrs.create(name = "gain", data=tx_gain)
        grp.attrs.create(name = "Ntones", data=1)
        f.close()
        
        print "samples saved on file."
        if(analyze): 
            samples = USRP_openH5file(filename)


    else:

        USRP_send_wait(command)
        
        USRP_filewait(5, USRP_rxfilename(0))
        
        #this is due to the buffer of the disk....?
        time.sleep(3)
    
        print "samples saved on file"
    
    if(analyze):
    
        #the following loop interpret 0 as an element
        n_points -=1
        
        #build the frequency axis
        step_hz = (last_f-start_f)/float(n_points)
        frequency = np.asarray([tone + start_f + f*step_hz for f in range(0, n_points)])
        
        data_z = [0 for f in range(0, n_points)]
        cutted = 0
        #f
        for f in range(0, n_points):
        
            if(USRP_data_service_connected):
            
                z = samples[0][ppt*f:ppt*f+ppt]#channel= 0
                
            else:
            
                z = USRP_loadfile(0, start_sample = ppt*f, lengh = ppt)
            
            
            if(ppt>59):
                low_cutted_index, high_cutted_index = USRP_cut_delay_fx(z, 1, ppt/20.)
                cutted += (float(high_cutted_index-low_cutted_index)/float(len(z)))
                data_z[f] = np.mean(z[low_cutted_index:high_cutted_index])
            else:
                data_z[f] = np.mean(z)
        if(ppt>59):
            print "Average cut un single tones: "+str(cutted/float(n_points))+"% of the samples"    

        
        return frequency,np.asarray(data_z)*USRP_calibration/dbm2vrms(USRP_outpower+tx_gain)
    else:
        time.sleep(0.4)
        return filename
        
#returns S21 samples and write vna group to HDF5 file
def USRP_VNA_analysis(filename):

    global USRP_calibration
    global USRP_outpower

    samples = USRP_openH5file(filename)

    #recover attributes from VNA file
    #raw_data contains untouched data from USRP server
    fv = h5py.File(filename+".h5",'r+')
    try:
        meas = fv["raw_data"]
    except KeyError:
        print "ERROR: Cannot find the raw_data group in the HDF5 file."
        return

    rate = meas.attrs.get("rate")
    start_f = meas.attrs.get("start_f")
    last_f = meas.attrs.get("last_f")
    n_points = meas.attrs.get("n_points")
    ppt = meas.attrs.get("ppt")
    tx_gain = meas.attrs.get("gain")

    if rate == None or start_f == None or n_points == None:
        print "ERROR: the HDF5 file is missing VNA scan attributes."   
        return
        
    #build the frequency axis
    step_hz = (last_f-start_f)/float(n_points)
    frequency = np.asarray([start_f + f*step_hz for f in range(0, n_points)])
    
    data_z = [0 for f in range(0, n_points)]
    cutted = 0
    
    for f in range(0, n_points):
        z = samples[0][ppt*f:ppt*f+ppt]
        if(ppt>1000):
            low_cutted_index, high_cutted_index = USRP_cut_delay_fx(z, 1, ppt/20.)
            cutted += (float(high_cutted_index-low_cutted_index)/float(len(z)))
            data_z[f] = np.mean(z[low_cutted_index:high_cutted_index])
        else:
            data_z[f] = np.mean(z)
            
    S21_data = np.asarray(data_z)*USRP_calibration/dbm2vrms(-USRP_outpower+tx_gain)
    
    try:
        vna = fv.create_group("vna")
    except ValueError:
        print "File was already analyzed. Overwrite was not implemented."
        return
    vna.attrs.create(name = "reso", data=np.abs(start_f-last_f)/float(n_points))
    S21 = vna.create_dataset("S21", data = S21_data )
    freq = vna.create_dataset("freq", data = frequency )
    
    fv.close()
    print "VNA analysis complete."
    return freq, S21

    
#Acquire signal from a list o tones
#returns the filename of the acquisition
def USRP_getnoise(tones, rate, time, decimation, tx_gain, rx_gain, RFtone = -1, ampl = [-1,]):
    global USRP_LineDelay
    global USRP_connected
    global USRP_LineDelay_Rate
    global USRP_socket
    #global USRP_data_queue
   

    if(not USRP_connected):
        print "ERROR: USRP_S is not connected!"
        return 0,[(0,0),]

    try:
        print "demodulating " + str(len(tones)) + "tones"
    except:
        tones = [tones,]
        
    #normalize the amplitude array if given
    if ampl[0] != -1:
        if len(tones) != len(ampl):
            print "WARNING: incorrect number of amplitudes ("+str(len(ampl))+") given!"
            print "Amplitudes will be flattened to 1/"+str(len(tones))
            ampl = np.ones(len(tones))

        else:
            print "Amplitudes will be normalized."
        ampl_norm = 0
        for i in range(len(tones)):
            ampl_norm += ampl[i]
        for i in range(len(tones)):
            ampl[i] /= ampl_norm
        
    #or equally divide the power by default
    else:
        ampl = np.zeros(len(tones))
        ampl.fill(1./len(tones))
        
    #by default the RFtone is the average of the tones
    if RFtone == -1:
        RFtone = np.mean(tones)
        RFtone = int(RFtone/1.25e6)*1.25e6
        print "RF tone is: "+str(RFtone)+" Hz"
        for i in range(len(tones)):
            tones[i] -= RFtone
    elif(int(RFtone/1.25e6)*1.25e6 != RFtone):
        print "WARNING: tone selected may (will) create a frequency shift in demodulation."

    expected_samples = float(rate)*time

    command = "tx on rx on"
    command += " tx_samples:" + str(expected_samples)
    command += " rx_samples:" + str(expected_samples)
    command += " rate:"+str(rate)
    command += " tone:"+str(RFtone)
    command += " rx_gain:"+ str(rx_gain)
    command += " tx_gain:"+ str(tx_gain)
    command += " tx_delay:" + str(1-USRP_LineDelay)
    command += " rx_delay:" + str(1)
    
    for i in range(len(tones)):
        command += " wave type:SINE"
        command += " ampl:"+str(ampl[i])
        command += " freq:"+str(tones[i])

    for i in range(len(tones)):
        command += " dem type:SINE"
        command += " freq:"+str(tones[i])
        command += " decim:"+str(decimation)
    
    global USRP_data_service_connected
    
    if(USRP_data_service_connected):


        USRP_send_nowait(command)
        USRP_resetflags()
        
        filename = "USRP_Noise_"+str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
        USRP_live_lan2file(filename, wait4finish = True, max_samp = expected_samples/max(decimation,1))
        
        fv = h5py.File(filename+".h5",'r+')
        try:
            meas = fv["raw_data"]
        except KeyError:
            print "ERROR: Cannot find the raw_data group in the HDF5 file."
            return
            
        grp = fv["raw_data"]
        grp.attrs.create(name = "rate", data=rate)
        grp.attrs.create(name = "decim", data=decimation)
        grp.attrs.create(name = "power", data=ampl)
        grp.attrs.create(name = "freq", data=tones)
        grp.attrs.create(name = "rf", data=RFtone)
        grp.attrs.create(name = "tx_gain", data=tx_gain)
        grp.attrs.create(name = "rx_gain", data=rx_gain)
        
        fv.close()

        #samples = USRP_openH5file(filename)
        
        USRP_resetflags()

        print "Samples stored and tagged."
        
        return filename
             
    else:
    
        
        print "ERROR: No connection with USRP_S. the offline functionality has to be implemented"
        return "ERROR"


#set the RX calibration constant so that a VNA scan with the same configuration should return (attenuator) dB
def USRP_calibrate(RFtone, tx_gain, attenuator = 0, USRP_maxpower = USRP_outpower):
    print "Calibrating: assuming USRP is in loopback."
    global USRP_LineDelay
    global USRP_connected
    global USRP_socket
    global USRP_calibration

    

    if(not USRP_connected):
        print "ERROR: USRP_S is not connected!"
        return -1
    
    #a reasonable ammount
    expected_samples = 1e6
    
    command = "tx on rx on"
    command += " tx_samples:" + str(expected_samples)
    command += " rx_samples:" + str(expected_samples)
    command += " rate:"+str(1e7)
    command += " tone:"+str(RFtone)
    command += " rx_gain:"+ str(0)
    command += " tx_gain:"+ str(tx_gain)
    command += " tx_delay:" + str(1-USRP_LineDelay)
    command += " rx_delay:" + str(1)
    
    command += " wave type:SINE"
    command += " ampl:"+str(1)
    command += " freq:"+str(1e5)


    #cut evetual delay/transient
    cutoff = int(2e5)

    global USRP_data_service_connected
    
    if(USRP_data_service_connected):

        USRP_send_nowait(command)
        
        USRP_live_lan2file("calibration", max_samp = expected_samples, wait4finish = True)
        
        data = (USRP_openH5file("calibration", channel= 0))[cutoff:-cutoff]

	USRP_resetflags()
        
    else:
    
        USRP_send_wait(command)
        
        data = (USRP_loadfile(0))[cutoff:-cutoff]
    
    ADC_ampl_rms = np.mean(np.absolute(data)/np.sqrt(2))
    
    USRP_calibration =  dbm2vrms(USRP_maxpower - attenuator)/ADC_ampl_rms
    
    print "Calibearion constant set to: " + str(USRP_calibration) + " Vrms/ADCunit"
    
    return USRP_calibration
    
def USRP_fast_chirp(tx_gain, n_points, start_f, last_f, ppt = None, tone = None, ringdown = 0.1):
    print "CHIRP started"
    if ppt == None:
        print "WARNING: Autoset Point Per Tone to 100"
        ppt = 100
    global USRP_LineDelay
    global USRP_connected
    global USRP_LineDelay_Rate
    global USRP_socket
    
   
    
    if(not USRP_connected):
        print "ERROR: USRP_S is not connected!"
        return [(0,0),]
    
    n_points = int(n_points)

    delta_f = max(start_f,last_f)-min(start_f,last_f)
    
    if(tone == None):
        tone = int((min(start_f,last_f) + (delta_f)/2)/1.25e6)*1.25e6
        print "RF tone is: "+str(tone)+" Hz"
        
        start_f -= tone
        last_f -= tone 
    
    maxf = max(np.abs(start_f),np.abs(last_f))
    
    if maxf < 5e6:
        rate = 1e7
        decimation = 100
    elif maxf < 1e7:
        rate = 2e7
        decimation = 200
    elif maxf < 2e7:
        rate = 5e7
        decimation = 500
    elif maxf < 5e7:
        rate = 1e8
        decimation = 1000
    else:
        print "ERROR: Desired frequency interval bigger than bandwidth (100 MHz)!"
        return [(0,0),]
        
    print "Rate is: "+str(rate)+" MHz"
    
    USRP_delayset(tone, rate)
    
    lapse = float(n_points*ppt * decimation) / float(rate)
    print "meas time is: "+str(lapse)+" s"

    #before decimation
    expected_samples = lapse * rate
    
    command = "tx on rx on"
    command += " tx_samples:" + str(expected_samples)
    command += " rx_samples:" + str(expected_samples+rate*ringdown)
    command += " rate:"+str(rate)
    command += " tone:"+str(tone)
    command += " rx_gain:"+ str(0)
    command += " tx_gain:"+ str(tx_gain)
    command += " tx_delay:" + str(1-USRP_LineDelay)
    command += " rx_delay:" + str(1)
    command += " bw:" + str(max(start_f,last_f)+20e6)

    command += " wave type:SWIPE"
    command += " ampl:"+str(1)
    command += " freq:"+str(start_f)
    command += " to:"+str(last_f)
    command += " lapse:"+str(lapse)
    command += " steps:"+str(int(n_points))


    global USRP_data_service_connected
    
    if(USRP_data_service_connected):

        USRP_send_nowait(command)
        
        filename = "USRP_VNA_"+str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
        
        USRP_live_lan2file(filename, wait4finish = True,max_samp = expected_samples/decimation-1) #max_samp = expected_samples/decimation
        
        USRP_resetflags()

        f = h5py.File(filename+".h5",'r+')
        grp = f["raw_data"]
        grp.attrs.create(name = "rate", data=rate)
        grp.attrs.create(name = "lapse", data=lapse)
        grp.attrs.create(name = "decim", data=decimation)
        grp.attrs.create(name = "start_f", data=tone+start_f)
        grp.attrs.create(name = "last_f", data=tone+last_f)
        grp.attrs.create(name = "n_points", data=n_points-1)
        grp.attrs.create(name = "ppt", data=ppt)
        grp.attrs.create(name = "gain", data=tx_gain)
        grp.attrs.create(name = "Ntones", data=1)
        f.close()
        
        print "samples saved on file."


    else:

        USRP_send_wait(command)
        
        USRP_filewait(5, USRP_rxfilename(0))
        
        #this is due to the buffer of the disk....?
        time.sleep(3)
    
        print "samples saved on file"
    

    return filename
            






