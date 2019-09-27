#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import serial
import io
import time
import pdb

class LEDMapper():
    Nrows = 6
    Ncols = 6
    def __init__(self, serial_port, baud_rate):
        self.arduino = serial.Serial(serial_port, baud_rate, timeout=10, write_timeout=10)
        if not self.arduino.is_open:
            raise IOError("No device found on %s"%(serial_port))
    
    # Need to initialize the device because opening the serial connection causes the arduino to reset.
    # This makes the program halt until it receives its first input from the arduino which tells the system
    # that it is ready to send and receive data. Seems to be system dependent. Happens on my macbook. Not sure
    # about linux yet.
    def init(self):
        while True:
            if (self.arduino.in_waiting > 0):
                #state = self.arduino.read_until(terminator=b'\r\n').decode()
                state = self.arduino.readline().decode()
                print ("Received back the string: ", state)
                break

    def led_on(self, row, col):
        cmd = b"CONST %d %d"%(row, col)
        num_bytes = self.arduino.write(cmd)
        print ("Number of bytes written ", num_bytes)
        self.arduino.flush()

    def get_state(self):
        state = self.arduino.read_until(terminator=b'\r\n').decode('ascii')
        return state

    def led_off(self, row, col):
        cmd = b"OFF %d %d"%(row, col)
        num_bytes = self.arduino.write(cmd)
        print ("Number of bytes written ", num_bytes)
        self.arduino.flush()

    def pulse_led(self, row, col):
        pass


if __name__=="__main__":
    mapper = LEDMapper('/dev/cu.usbmodem14101', 19200)
    mapper.init()
    #pdb.set_trace()
    for row in range(LEDMapper.Nrows):
        for col in range(LEDMapper.Ncols):
            mapper.led_on(row, col)
            time.sleep(1)
            print (mapper.get_state())
            mapper.led_off(row, col)
            print (mapper.get_state())

