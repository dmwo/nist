# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:50:41 2019

@author: dmo2
"""

#%%

import math
import numpy as np
import pyvisa as visa
import matplotlib.pyplot as plt

res_meas = True

class Agilent89410a:
    def __init__(self, visa_name, trace=1):
        if not 1<=int(trace)<=4:
            raise Exception('Trace must be integer between 1 and 4')
        self.rm = visa.ResourceManager()
        self.visa = self.rm.open_resource(visa_name)
        self.trace = str(int(trace))
        self.timeout = 5000 # Set response timeout (in milliseconds)

    def read(self):
        return self.visa.read()
    
    def write(self, string):
        self.visa.write(string)

    def query(self, string):
        return self.visa.query(string)
    
    def identify(self):
        # expected output: HEWLETT-PACKARD,89410A,3416A01211,A.09.01
        return self.query('*IDN?')
    
    def close(self):
        self.close()

    def cal(self):
        # calibrates analyzer and returns pass/fail result
        return self.query('*CAL?')

    def reset(self):
        # executes a device reset and cancels pending commands/queries
        return self.write('*RST')

    def trigger(self):
        # triggers analyzer when GPIB is designated trigger source and analyzer
        # is waiting to trigger
        return self.write('*TRG')

    def cont(self):
        # continues a paused measurement
        return self.write('CONTINUE')

    def pause(self):
        # Pauses measurement in progress
        return self.write('PAUSE')

    def abort(self):
        # aborts measurement in progress and resets trigger system
        return self.write('ABORT')

    def wait(self):
        # device will not process subsequent commands until preceding ones are 
        # complete
        return self.write('*WAI')

    def preset(self):
        # resets most device parameters to default state
        return self.write('SYSTEM:PRESET')

    def set_free_run_mode(self):
        # sets triggering to free-run mode
        return self.write('TRIG:SOUR IMM')

    def set_detector(self, f='SIGNAL'):
        # determines what analyzer displays when number of data points exceeds
        # displayed points
        if f not in ['SIGNAL','SAMPLE','POSITIVE']:
            raise Exception('invalid input')
        return self.write('SENSE:DETECTOR:FUNCTION '+f)

    def get_ccdf_count(self):
        # returns current number of data samples in CCDF measurement
        return self.query('CALCULATE'+self.trace+':CCDF:COUNT?')

    def get_ccdf_power(self):
        # returns average signal power used to compute CCDF measurement
        return self.query('CALCULATE'+self.trace+':CCDF:POWER?')

    def set_fast_average(self, val=True):
        # Turns fast averaging on or off
        return self.write('AVERAGE:IRES '+'1' if val else '0')

    def get_average_progress(self):
        # returns current number of averages
        return self.query('SENSE:AVERAGE:COUNT:INTERMEDIATE?')

    def get_average_total(self):
        # returns number of averages to be done
        return self.query('SENSE:AVERAGE:COUNT?')

    def set_measurement(self, mode='S', measurement='P'):
        # Sets measurement data to be displayed.
        # Modes:
            # Scalar mode: 'S'
            # Vector mode: 'V'
        # Measurements:
            # Spectrum: 'S'
            # PSD: 'P'
            # Frequency response: 'F'
        s = ''
        if mode is 'S':
            if measurement is 'S':
                s = 'XFR:POW 1'
            elif measurement is 'P':
                s = 'XFR:POW:PSD 1'
            else:
                raise Exception('Invalid measurement in scalar mode')
        elif mode is 'V':
            if measurement is 'S':
                s = 'XFR:POW 1'
            elif measurement is 'P':
                s = 'XFR:POW:PSD 1'
            elif measurement is 'F':
                s = 'XFR:POW:RAT 2,1'
            else:
                raise Exception('Invalid measurement in vector mode')
        else:
            raise Exception('Invalid mode')
        return self.write('CALCULATE'+self.trace+':FEED \''+s+'\'')

    def get_mode(self):
        # returns string representing current measurement and mode
        return self.query('CALCULATE'+self.trace+':FEED?')

    def set_trace_coords(self, coords):
        # sets coordinates of selected trace
        # Coordinate settings:
            # 'MLIN': linear magnitude
            # 'MLOG': log-magnitude
            # 'PHAS': wrapped phase coordinates
            # 'UPH': unwrapped phase coordinates
            # 'REAL': real part of waveform
            # 'IMAG': imaginary part of waveform
            # 'GDEL': group delay
            # 'COMP': complex polar vector diagram
            # 'CONS': constellation polar vector diagram
            # 'IEYE': in-phase eye diagram
            # 'QEYE': quadrature-phase eye diagram
            # 'TEYE': trellis eye diagram
        if coords not in ['MLIN','MLOG','PHAS','UPH','REAL','IMAG', 'DGEL',\
                          'COMP','CONS','IEYE','QEYE','TEYE']:
            raise Exception('Invalid coordinates')
        return self.write('CALCULATE'+self.trace+':FORMAT '+coords)

    def get_trace_coords(self):
        # returns string representing current trace coordinates
        return self.query('CALCULATE'+self.trace+':FORMAT?')

    def set_marker_band_position(self, mkr, pos):
        # sets left or rignt band marker position
        # mkr: 'L' for left. 'R' for right marker
        # pos: x-axis units, e.g. Hz or s
        if mkr not in ['L','R']:
            raise Exception('Invalid mkr')
        return self.write('CALCULATE'+self.trace+':MARKER:BAND:'+\
                         ('START ' if pos is 'L' else 'STOP ')+str(pos))

    def set_marker_coupling(self, val):
        # turn on/off marker coupling; val should be boolean
        return self.write('CALCULATE'+self.trace+'MARKER:COUPLED '+\
                         ('ON' if val else 'OFF'))

    def set_freq_center(self, val):
        # turn on/off marker frequency counter; val should be boolean
        return self.write('CALCULATE'+self.trace+':MARKER:FCOUNT '+\
                         ('ON' if val else 'OFF'))
    def get_freq_center(self):
        # returns frequency counter measurement
        return self.query('CALCULATE'+self.trace+':MARKER:FCOUNT:RESULT?')

    def move_marker(self, move):
        # Moves marker based on value of move variable
        # 'MAX': moves marker to maximum value in trace
        # 'MIN': moves marker to minimum value in trace
        # 'LEFT': moves marker to nearest local maximum to left
        # 'RIGHT': moves marker to nearest local maximum to right
        # 'NEXT': moves marker to next-highest value in trace
        if move not in ['MAX','MIN','LEFT','RIGHT','NEXT']:
            raise Exception('Invalid movement')
        return self.write('CALCULATE'+self.trace+':MARKER:'+\
                                ('MAXIMUM' if move is not 'MIN' else 'MINIMUM')+\
                                (move if move not in ['MIN','MAX'] else ''))
    def set_tracking(self, val):
        # Turns marker peak-tracking on or off. val should be boolean.
        return self.write('CALCULATE'+self.trace+':MARKER:MAXIMUM:TRACK '+\
                                ('ON' if val else 'OFF'))

    def set_coupling(self, coupling = 'AC', channel=1):
        # Sets selected channel to AC (coupling='AC') or DC (coupling='DC') coupling
        if channel not in [1,2] or coupling not in ['AC','DC']:
            raise Exception('Invalid input')
        return self.write('INPUT'+str(channel)+':COUPLING '+coupling)

    def set_input_impedance(self, ohms = 50, channel=1):
        # Sets input impedance. imp should be integer 50, 75, or 1e6 and channel
        # should be 1 or 2.
        if channel not in [1,2] or ohms not in [50,75,1e6]:
            raise Exception('Invalid input. ohms should be integer 50, 75, or 1e6 and channel')
        return self.write('INPUT'+str(channel)+':IMPEDANCE '+str(ohms))

    def set_channel(self, val, channel=1):
        # Turns selected channel on or off. val should be boolean and channel
        # should be 1 or 2.
        if channel not in [1,2]:
            raise Exception('Invalid input')
        return self.write('INPUT'+str(channel)+' '+\
                                ('ON' if val else 'OFF'))

    def set_mode(self, mode):
        # Sets instrument mode. mode should be 'SCALAR' or 'VECTOR'.
        if mode not in ['SCALAR','VECTOR']:
            raise Exception('Invalid input')
        return self.write('INSTRUMENT '+mode)

    def set_average(self, val, tp='VRMS', count=100):
        # Turns averaging on or off based on boolean val. If val is True, 
        # type (tp) can be :
            # 'VRMS' for rms(video)
            # 'ERMS' for rms(video) exponential
            # 'TIME' for time
            # 'TEXP' for time exponential
            # 'CPH' for continuous peak hold
        # Count is number of averages to take; should be integer between 1 and 99999. 
        if not val:
            return self.write('AVERAGE OFF')
        if tp not in ['VRMS','ERMS','TIME','TEXP','CPH'] or not 1<=count<=99999:
            raise Exception('invalid input')
        s = ''
        if tp is 'VRMS':
            s = 'RMS;TCON NORM'
        elif tp is 'ERMS':
            s = 'RMS;TCON EXP'
        elif tp is 'TIME':
            s = 'COMP;TCON NORM'
        elif tp is 'TEXP':
            s = 'COMP;TCON EXP'
        elif tp is 'CPH':
            s = 'MAX'
        return (self.write('AVERAGE ON'),\
                self.write('AVERAGE:TYPE '+s),\
                self.write('AVERAGE:COUNT '+str(count)))

    def set_freq(self, span, center):
        # Specifies resolution span and center. Both should be integers [Hz].
        return (self.write('FREQUENCY:SPAN '+str(span)),\
                self.write('FREQUENCY:CENTER '+str(center)))

    def set_channel_gain(self, gain, channel=1):
        # Sets gain for selected channel; gain should be number between 1e-6 and 1e6
        if channel not in [1,2] or not 1e-6<=gain<=1e6:
            raise Exception('Invalid input')
        return self.write('CORRECTION'+str(channel)+':LOSS:MAGNITUDE '+str(gain))

    def set_dc_offset(self, offset, channel=1):
        # Sets DC offset for selected channel; channel should be 1 or 2 and 
        # offset should be between -20 and 20
        if channel not in [1,2] or not -20<=offset<=20:
            raise Exception('Invalid input')
        return self.write('CORRECTION'+str(channel)+':OFFS '+str(offset))

    def get_x_data(self):
        # Returns x-axis units and data
        ret = []
        ret.append(self.query('TRACE:X:UNIT? TRACE'+self.trace).strip().replace('"',''))
        x = [float(n) for n in self.query('TRACE:X? TRACE'+self.trace).split(',')]
        ret.append([n for n in x if 0 <= n <= float(self.query('SENS:FREQ:STOP?')) + 1e-3])
        return ret

    def get_y_data(self):
        # Returns y-axis units and data
        ret = []
        ret.append(self.query('CALCULATE'+str(self.trace)+':UNIT:POWER?').strip().replace('"',''))
        ret.append([float(n) for n in self.query('CALCULATE'+str(self.trace)+':DATA?').split(',')])
        return ret

def plot(x, y):
    plt.plot(x[1], y[1])
    
    # Titling and labelling
    plt.title('Power Spectral Density')
    plt.xlabel('Frequency ({})'.format(x[0]))
    plt.ylabel('Magnitude ({})'.format(y[0]))
    plt.show()
    
def impedance_plot(R, C, x, y):
    f = x[1]; y = y[1]; kb = 1.38e-23; T = 300
    zero_C = math.sqrt(4 * kb * T * R)
    Z = [[R / math.sqrt(1 + (2*math.pi*n*R*i) ** 2) for n in f] for i in C]
    Z = [[zero_C * n / max(Z[i]) for n in Z[i]] for i in range(len(Z))]
    idx = next(n[0] for n in enumerate(y) if n[1] < zero_C)
#    for i in range(len(Z)):
#      for j in range(idx-5):
#        Z[i].insert(0, zero_C)
#      Z[i] = Z[i][:-(idx-5) or None]
    y[:idx] = [zero_C] * idx
    plt.figure()
    for i in range(len(Z)): plt.plot(f, Z[i], label = 'C = {} nF'.format(C[i] * 1e9))
    plt.scatter(f, y, s = 2**1, label = 'PSD Curve')
    plt.legend()
    plt.title('Effect of Various Capacitances on Parallel RC Total Impedance')
    plt.xlabel('f (Hz)')
    plt.ylabel('Z (Ohms)')
    plt.show()
    
try: sa.close()
except: pass
sa = Agilent89410a(visa_name = 'GPIB0::19::INSTR')
sa.identify() # The spectrum analyzer should respond now

if res_meas:
    x = sa.get_x_data()
    y = sa.get_y_data()
    plot(x, y)

#%%
kb = 1.38e-23
T = 300
R = 1e6