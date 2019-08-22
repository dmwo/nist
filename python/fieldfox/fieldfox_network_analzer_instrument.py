# -*- coding: utf-8 -*-
# Function: FieldFox Network Analyser Communication
# Author: Dylan Oh

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import visa
import pickle
import time

visa_name = "TCPIP::%s::INSTR" % '169.254.112.226'

class AgilentN9912A(object):
    def __init__(self):
        visa_name = "TCPIP::%s::INSTR" % '169.254.112.226'
        self.rm = visa.ResourceManager()
        self.pyvisa = self.rm.open_resource(visa_name)
        self.pyvisa.timeout = 5000 # Set response timeout (in milliseconds)

    def write(self, string):
        self.pyvisa.write(string)

    def query(self, string):
        return self.pyvisa.query(string)

    def reset(self):
        self.write('*RST')

    def identify(self):
        return self.query('*IDN?')

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
        return (f_start, f_stop, num_points)
        
    # Changes the power (0 to -31 dBm in 1 dBm steps)
    def change_power(self, power_dBm = -31):
        self.write('SOUR:POW %i' % power_dBm)

    # Changes format. Input 'MLOG' for dB, 'MLIN' for relative
    def change_format(self, measure_format):
        self.write('CALC:FORM %s' % measure_format)

    # Pull whatever S21 or S11 data is on the screen - in terms of dB or absolute magnitude
    def measure(self, frequencies, trace = 1):
        num_points = frequencies[2]
        mags = []
        
        # Setup Network Analyser
        # Disables hold if one exists
        self.write('INIT:CONT 1')
        # Select trace to measure
        self.write('CALC:PAR%i:SEL' % trace)
            
        # Read data and parse into a list of floats
        mag_data = self.query('CALC:DATA:FDATa?')
        mag_data = mag_data.split(",")
        for i in range(0, num_points):
            mags.append(float(mag_data[i]))
                
        # Creating a vector of frequencies to correspond to magnitude data
        freq_range = frequencies[1] - frequencies[0]
        freq_increment = freq_range / (num_points - 1)
        
        # Initialising variable "freqs" with num_points number of elements
        freqs = []
        for i in range(0, num_points):
            freqs.append(frequencies[0] + freq_increment * i)

        # Return values
        return (freqs, mags)

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
    
def plot_setup(mode, fname):
    plt.title('Magnitude vs Frequency in {0} Mode'.format(mode))
    plt.legend(loc = 'best')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBm)')
    plt.show()
    plt.savefig('{0}.png'.format(fname), format = 'png')

#%%

na = AgilentN9912A()
na.reset()
time.sleep(5)
na.set_mode('S11')
na.change_format('MLOG')
time.sleep(3)
freqs = na.set_freq_range()

dataS1133 = na.measure(freqs)
na.set_mode('S21')
time.sleep(3)
dataS2133 = na.measure(freqs)

input('Change trace connection to crosstalk')
dataS2134 = na.measure(freqs)
na.set_mode('S11')
time.sleep(3)
dataS1134 = na.measure(freqs)

input('Change trace connection to 4-4')
dataS1144 = na.measure(freqs)
na.set_mode('S21')
time.sleep(3)
dataS2144 = na.measure(freqs)

figS11direct = plt.figure()
freqS1133 = dataS1133[0]; magS1133 = dataS1133[1]
freqS1144 = dataS1144[0]; magS1144 = dataS1144[1]
trace1 = plt.plot(freqS1133, magS1133, label = 'P3 to P3')
trace2 = plt.plot(freqS1144, magS1144, label = 'P4 to P4')
plot_setup('S11', 'S11_direct')

figS21direct = plt.figure()
freqS2133 = dataS2133[0]; magS2133 = dataS2133[1]
freqS2144 = dataS2144[0]; magS2144 = dataS2144[1]
trace1 = plt.plot(freqS2133, magS2133, label = 'P3 to P3')
trace2 = plt.plot(freqS2144, magS2144, label = 'P4 to P4')
plot_setup('S21', 'S21_direct')

figS11cross = plt.figure()
freqS1134 = dataS1134[0]; magS1134 = dataS1134[1]
trace1 = plt.plot(freqS1134, magS1134, label = 'P3 to P4')
plot_setup('S11', 'S11_cross')

figS21cross = plt.figure()
freqS2134 = dataS2134[0]; magS2134 = dataS2134[1]
trace1 = plt.plot(freqS2134, magS2134, label = 'P3 to P4')
plot_setup('S21', 'S21_cross')

figS11both = plt.figure()
trace1 = plt.plot(freqS1133, magS1133, label = 'P3 to P3')
trace2 = plt.plot(freqS1144, magS1144, label = 'P4 to P4')
trace3 = plt.plot(freqS1134, magS1134, label = 'P3 to P4')
plot_setup('S11', 'S11_both')

figS21both = plt.figure()
trace1 = plt.plot(freqS2133, magS2133, label = 'P3 to P3')
trace2 = plt.plot(freqS2144, magS2144, label = 'P4 to P4')
trace3 = plt.plot(freqS2134, magS2134, label = 'P3 to P4')
plot_setup('S21', 'S21_both')

pickle.dump(dataS1133, open('S11_3-3.pickle', 'wb'))
pickle.dump(dataS1144, open('S11_4-4.pickle', 'wb'))
pickle.dump(figS11direct, open('S11_direct.fig.pickle', 'wb'))
pickle.dump(dataS2133, open('S21_3-3.pickle', 'wb'))
pickle.dump(dataS2144, open('S21_4-4.pickle', 'wb'))
pickle.dump(figS21direct, open('S21_direct.fig.pickle', 'wb'))
pickle.dump(dataS1134, open('S11_3-4.pickle', 'wb'))
pickle.dump(figS11cross, open('S11_cross.fig.pickle', 'wb'))
pickle.dump(dataS2134, open('S21_3-4.pickle', 'wb'))
pickle.dump(figS21cross, open('S21_cross.fig.pickle', 'wb'))
pickle.dump(figS11both, open('S11_both.fig.pickle', 'wb'))
pickle.dump(figS21both, open('S21_both.fig.pickle', 'wb'))


#%%

dataS1133 = pickle.load(open('S11_3-3.pickle', 'rb'))
dataS1144 = pickle.load(open('S11_4-4.pickle', 'rb'))
dataS2133 = pickle.load(open('S21_3-3.pickle', 'rb'))
dataS2144 = pickle.load(open('S21_4-4.pickle', 'rb'))
dataS1134 = pickle.load(open('S11_3-4.pickle', 'rb'))
dataS2134 = pickle.load(open('S21_3-4.pickle', 'rb'))

#%%

figS11direct = pickle.load(open('S11_direct.fig.pickle', 'wb'))
figS21direct = pickle.load(open('S21_direct.fig.pickle', 'wb'))
figS11cross = pickle.load(open('S11_cross.fig.pickle', 'wb'))
figS21cross = pickle.load(open('S21_cross.fig.pickle', 'wb'))
figS11both = pickle.load(open('S11_both.fig.pickle', 'wb'))
figS21both = pickle.load(open('S21_both.fig.pickle', 'wb'))