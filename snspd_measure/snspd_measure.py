# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:21:23 2018

@author: Dylan Oh
"""
#%%
from __future__ import division

import time
import visa
import pickle
import datetime
import numpy as np
import matplotlib.pyplot as plt

#%%

"""
Importing Instrument Classes

"""

from instruments.agilent_53131a import Agilent53131a
from instruments.srs_sim928     import SIM928
from instruments.srs_sim970     import SIM970
from instruments.jds_ha9        import JDSHA9
from instruments.lecroy_620zi   import LeCroy620Zi
from instruments.rigol_dg5000   import RigolDG5000

#%%

"""
Measurement Functions

"""


def iv_curve(v_min = 0.0,
             v_max = 0.5,
             Z_in = 1e5,
             num_meas = 501):

    """ Calculates values germane to producing an i-V curve of an SNSPD using
    the SRS SIM928 voltage source and SIM970 digital voltmeter.

    Parameters
    ----------
    v_min    : float, startpoint of voltage ramp (generally 0)
    v_max    : float, peak of voltage ramp
    Z_in     : float, value of current-limiting resistor
    num_meas : int  , number of values between start and peak (actual
                      experiment will run 4 times the length of this number)

    Returns
    -------
    (V,I) : tuple, voltage across DUT and bias current
    """

    # Initialising
    s28 = SIM928('GPIB0::4::INSTR', 4)
    s70 = SIM970('GPIB0::4::INSTR', 7)

    v_in = np.linspace(v_min, v_max, num_meas)
    v_in = np.append(v_in, np.flip(v_in, 0))     # Sweep voltages (triangle wave)
    v_in = np.append(v_in, -v_in)                # Adding negative voltages
    I = [];  V = []                              # Bias current and V device
    s28.reset();  s70.reset()
    s70.set_impedance(gigaohm = True, channel = 2)

    # Measuring
    for i in range(len(v_in)):
        s28.set_voltage(voltage = v_in[i])       # Setting sweep voltage [V]
        time.sleep(1)
        v = s70.read_voltage(channel = 1)        # Reading DVM 1 [V]
        V.append(s70.read_voltage(channel = 2))  # Reading DVM 2 [V]
        I.append((v - V[i]) * 1e6 / Z_in)        # Calculating current [μA]

    s28.close(); s70.close()
    return (V,I)


def dark_counts(v_min = 0.35,
                v_max = 0.5,
                Z_in = 1e5,
                num_meas = 50,
                t_int = 60,
                trig_volt = 0.03,
                attenuation = 0,
                show_plot = False,
                data_name = 'dark'):

    """ Measures dark counts detected for a ramp of I bias values using the SRS
    SIM928 voltage source and SIM970 digial voltmeter and the Agilent 53131A
    counter.

    Parameters
    ----------
    v_min       : float, startpoint of voltage ramp
    v_max       : float, peak of voltage ramp
    Z_in        : float, value of current-limiting resistor
    num_meas    : int  , number of values between start and peak
    t_int       : int  , integration time of the counter    
    trig_volt   : float, trigger voltage of the Agilent counter
    attenuation : int  , (optional) attenuation through JDSHA9 if using laser
    show_plot   : bool , show data plot within this function or return data

    Returns
    -------
    (I,counts) : tuple, bias current and number of counts
    """

    # Setting up instruments
    s28 = SIM928('GPIB0::4::INSTR', 4)
    s70 = SIM970('GPIB0::4::INSTR', 7)
    ag = Agilent53131a('GPIB0::5::INSTR')
    att = JDSHA9('GPIB0::15::INSTR')

    s28.reset();  s70.reset()
    ag.basic_setup()
    ag.set_trigger(trigger_voltage = trig_volt, slope_positive = True)
    if attenuation != 0:
        att.set_beam_block(beam_block = False)
        att.set_attenuation_db(attenuation)
    att.close()
    s28.set_voltage(0)
    s28.set_output(True)
    ag.pyvisa.timeout = t_int*1000 + 5000 # Stop timeout on longer integrations

    # Initialising vectors
    v_in = np.linspace(v_min, v_max, num_meas)      # Voltage ramp
    I = [];  counts = []                            # Bias current and counts

    # Counting
    for i in range(len(v_in)):
        s28.set_voltage(voltage = v_in[i])          # Setting ramp voltage
        time.sleep(0.5)                             # Letting voltage settle
        v = s70.read_voltage(channel = 1)           # Reading source voltage
        V = s70.read_voltage(channel = 2)           # Reading DUT voltage
        I.append((v - V) * 1e6 / Z_in)              # Calculating I bias
        counts.append(ag.timed_count(t_int))         # Number of counts
        s28.set_voltage(0)                          # Resetting I bias
        time.sleep(0.5)                             # Letting voltage settle

    s28.close(); s70.close(); ag.close()
    
    if show_plot:
        make_plot((I, counts),
                  xlab = 'Bias Current (uA)',
                  ylab = 'Countrate (Hz)',
                  title = 'Dark Counts',
                  d_name = 'dark')
    else:
        return (I, counts)


def ic_dist(wf_pts_t = [0.0, 20e-3, 100e-3, 101e-3],
            wf_pts_v = [0.0, 0.0, 5.0, 0.0],
            v_peak = 0.8,
            trig_source = 'C2',
            aux_source = 'C1',
            data_source = 'F1',
            trig_volt = 0.1,
            time_div = 1e-6,
            samples = 10000):

    """ Samples many switching currents of the SNSPD using the Rigol DG5000
    arbitrary waveform generator and the LeCroy 620Zi oscilloscope.

    Parameters
    ----------
    wf_pts_t    : list  , vector of time values of arbitrary waveform
    wf_pts_v    : list  , vector of voltage values of arbitrary waveform
    v_peak      : float , peak-to-peak voltage of waveform
    trig_source : string, trigger source of the oscilloscope
    data_source : string, data source of the oscilloscope
    trig_volt   : float , trigger voltage of the oscilloscope
    samples     : int   , number of samples to take

    Returns
    -------
    data : list, contains voltages from trend
    """

    # Setting up arbitrary waveform
    awg = RigolDG5000('USB0::0x1AB1::0x0640::DG5T171200124::INSTR')
    time.sleep(0.5)
    awg.reset()
    time.sleep(1)
    awg.setup_arb_wf(t = wf_pts_t, v = wf_pts_v)
    awg.set_vhighlow(vlow = 0, vhigh = v_peak)
    awg.set_freq(14)
    awg.set_impedance()
    awg.set_output(output = True)

    # Pulling trace from oscilloscope
    setup_lecroy(trig_source,
                 'DC1M',
                 aux_source,
                 'DC1M',
                 'F1',
                 'P1',
                 trig_volt,
                 1e-6,
                 trend = True)
    awg.close()


def setup_lecroy(trig_source = 'C1',
                 trig_couple = 'DC50',
                 aux_source = 'C2',
                 aux_couple = 'DC50',
                 data_source = 'F1',
                 meas_source = 'P1',
                 trig_volt = 0.03,
                 time_div = 1e-6,
                 trend = True):

    lecroy = LeCroy620Zi("TCPIP::%s::INSTR" % '192.168.1.100')
    lecroy.reset()
    lecroy.set_trigger(source = trig_source, volt_level = trig_volt)
    lecroy.set_coupling(channel = trig_source, coupling = trig_couple)
    lecroy.set_coupling(channel = aux_source, coupling = aux_couple)
    lecroy.set_horizontal_scale(time_per_div = time_div)
    lecroy.set_trigger_mode('Auto')

    if trend == True:
        lecroy.set_parameter(param_engine = 'Mean')
        lecroy.setup_math_trend(math_channel = data_source,
                                source = meas_source,
                                num_values = 10000)
    lecroy.close()


def get_trace(data_source = 'F1',
              meas_source = 'P1'):

    lecroy = LeCroy620Zi("TCPIP::%s::INSTR" % '192.168.1.100')
    lecroy.set_trigger_mode('Single')

    if data_source[0] == 'C':
        t, volt = (lecroy.get_wf_data(channel = data_source))
        t = t * 1e9;  volt = volt * 1e3
        data = (t, volt)
    elif data_source[0] == 'F':
        t, data = lecroy.get_wf_data(channel = data_source)

    lecroy.set_trigger_mode('Auto')
    lecroy.close()
    return data


#%%

"""
Data Visualisation Functions

"""

def make_plot(data,
         xlab = 'Voltage (V)',
         ylab = 'Current (uA)',
         title = 'i-V Curve of SNSPD',
         d_name = '',
         lty = [],
         bins = 1000):

    # Plotting
    fig = plt.figure()
    if type(data) == tuple:
        plt.plot(data[0], data[1], lty)
        plt.ylabel('%s' % ylab)
    else:
        plt.hist(data, bins)
        plt.ylabel('Frequency')
    plt.xlabel('%s' % xlab)
    plt.title('%s' % title)

    filename = datetime.datetime.now().strftime('%Y-%m-%d %H-%M-%S')
    filename = filename + ' ' + d_name
    pickle.dump({'%s' % d_name:data}, open(filename + '.pickle', 'wb'))
    pickle.dump(fig, open(filename + '.fig.pickle','wb'))
    fig.savefig(filename)

#%%

"""
Simulation Functions (Used with LTSpice)

"""

def sim_iv_curve(fname = 'SNSPD-stochastic-gamma.txt'):

    V = [];  I = []
    data = open(fname, 'r')

    for line in data:
        try:
            type(int(line[0])) == int
            V.append(float(line.split()[1]) * 1e6)
            I.append(float(line.split()[2]) * 1e6)
        except: pass
    return (V,I)

def sim_ic_dist(fname = 'SNSPD-stochastic-gamma.txt',
                num_bins = 1000,
                num_samples = 10000):

    I = [];  flag = 1
    raw = open(fname, 'r')

    for line in raw:
        if line[0] == 'S':
            flag = 0
        elif line[0] != 'S' and flag == 0:
            try:
                type(int(line[0])) == int
                if float(line.split()[1]) > 10e-6:
                    I.append(float(line.split()[2]) * 1e6)
                    flag = 1
            except: pass
        else: pass

    make_plot(I,
         bins = num_bins,
         xlab = 'Bias Current (uA)',
         ylab = 'Frequency',
         title = 'Distribution of Switching Currents Simulated, samples = %s bins = %s' 
                 % (str(num_samples), str(num_bins)))

    return I
