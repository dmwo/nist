# -*- coding: utf-8 -*-
# Function: FieldFox Network Analyser Communication
# Author: Dylan Oh

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import visa

visa_name = "TCPIP::%s::INSTR" % '192.168.6.100'

class AgilentN9912A(object):
    def __init__(self, visa_name):
        self.rm = visa.ResourceManager()
        self.pyvisa = self.rm.open_resource(visa_name)
        self.pyvisa.timeout = 5000 # Set response timeout (in milliseconds)

    def write(self, string):
        self.pyvisa.write(string)

    def query(self, string):
        return self.pyvisa.query(string)

    def reset(self):
        write('*RST')

    def identify(self):
        ID = self.query('*IDN?')
        print(ID)

    def close(self):
        self.pyvisa.close()

    def set_mode(self, mode):
        if mode is 'S11':
            # Set to Network Analyser Mode
            self.write('INST:SEL "NA"')
            # Disables hold if one exists
            self.write('INIT:CONT 1')
            # Set mode to S11
            self.write('CALC:PAR:DEF S11')
            
        elif mode is 'S21':
            # Set to Network Analyser Mode
            self.write('INST:SEL "NA"')
            # Disables hold if one exists
            self.write('INIT:CONT 1')
            # Set mode to S11
            self.write('CALC:PAR:DEF S21')
        
    # Sets frequency range
    def set_freq_range(self,
                       f_start = 10e6,
                       f_stop = 3e9,
                       num_points = 201,
                       f_centre = None,
                       f_span = None):

        # Sets frequency range by start and stop if centre and span are not defined
        if f_centre is None and f_span is None:
            self.write('FREQ:STAR %0.6e' % f_start)
            self.write('FREQ:STOP %0.6e' % f_stop)
            self.write('SWE:POIN %i' % num_points)
            
        # Sets frequency range by centre and span if defined
        elif f_centre is not None and f_span is not None:
            self.write('FREQ:CENT %0.6e' % f_centre)
            self.write('FREQ:SPAN %0.6e' % f_span)
            
        # Returning data
        return f_start, f_stop, num_points
        
    # Changes the power (0 to -31 dBm in 1 dBm steps)
    def change_power(self, power_dBm = -31):
        self.write('SOUR:POW %i' % power_dBm)

    # Changes format. Input 'MLOG' for dB, 'MLIN' for relative
    def change_format(self, measure_format):
        self.write('CALC:FORM %s' % measure_format)

    # Pull whatever S21 or S11 data is on the screen - in terms of dB or absolute magnitude
    def measure(self, frequencies, num_points, trace = 1):
        num_points = frequencies[2]
        mags = [None] * num_points
        
        # Setup Network Analyser
        # Disables hold if one exists
        self.write('INIT:CONT 1')
        # Select trace to measure
        self.write('CALC:PAR%i:SEL' % trace)
            
        # Read data and parse into a list of floats
        mag_data = self.query('CALC:DATA:FDATa?')
        mag_data = mag_data.split(",")
        for i in range(0, num_points):
            mags[i] = float(magnitude_data[i])
                
        # Creating a vector of frequencies to correspond to magnitude data
        freq_range = frequencies[1] - frequencies[0]
        freq_increment = freq_range / (num_points - 1)
        
        # Initialising variable "freqs" with num_points number of elements
        freqs = [None] * num_points
        for i in range(0, num_points):
            freqs[i] = frequencies[0] + freq_increment * i

        # Return values
        return freqs, mags

# Plot data
def plot_data(data, mode):
    # Extracting frequencies and magnitudes from 'data' tuple variable
    freqs = data[0]; mags = data[1]
    
    # Plotting
    plt.plot(freqs, mags)
    
    # Setting window
    xmin = min(freqs) * 0.98
    xmax = max(freqs) * 1.02
    ymin = min(mags) * 0.98
    ymax = max(mags) * 1.02
    plt.axis([xmin, xmax, ymin, ymax])
    
    # Titling and labelling
    plt.title('Magnitude vs Frequency in %s Mode' % mode)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBm)')
    plt.show()
