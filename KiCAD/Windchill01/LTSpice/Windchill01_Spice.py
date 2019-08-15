# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:08:54 2019

@author: dmo2
"""

import matplotlib.pyplot as plt

def sim_iv_curve(fname = 'Windchill01.txt'):

    t = []; Vth = []; Vthrc = []; Vouta = []; Voutb = []; Vrc = []; Vsig = []; I = []
    data = open(fname, 'r')

    for line in data:
        try:
            type(int(line[0])) == int
            t.append(float(line.split()[0]) * 1e3)
            Vthrc.append(float(line.split()[1]))
            Vth.append(float(line.split()[2]) * 1e3)
            Vouta.append(float(line.split()[3]))
            Voutb.append(float(line.split()[4]))
            Vrc.append(float(line.split()[5]))
            Vsig.append(float(line.split()[6]) * 1e3)
            I.append(abs(float(line.split()[7]) * 1e3))
        except: pass
    return (t, Vthrc, Vth, Vouta, Voutb, Vrc, Vsig, I)

data = sim_iv_curve()
t = data[0]
Vthrc = data[1]
Vth = data[2]
Voutb = data[3]
Vouta = data[4]
Vrc = data[5]
Vsig = data[6]
I = data[7]

#%%
    
plt.figure()
Vthline = plt.plot(t, Vth, label = 'Threshold Voltage A')
Vsigline = plt.plot(t, Vsig, label = 'Signal')
Voutaline = plt.plot(t, Vouta, label = 'Comparator A Output')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('Comparator A Response')
plt.legend(loc = 'best')
plt.show()

#%%

plt.figure()
Voutaline = plt.plot(t, Vouta, label = 'Comparator A Output')
Vrcline = plt.plot(t, Vrc, label = 'RC Capacitor Voltage')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (V)')
plt.title('RC Charging from Output A Voltage')
plt.legend(loc = 'best')
plt.show()

#%%

plt.figure()
Vthrcline = plt.plot(t, Vthrc, label = 'Threshold Voltage B')
Vrcline = plt.plot(t, Vrc, label = 'RC Capacitor Voltage')
Voutbline = plt.plot(t, Voutb, label = 'Comparator B Output')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (V)')
plt.title('Comparator B Response')
plt.legend(loc = 'best')
plt.show()

#%%

plt.figure()
Vsigline = plt.plot(t, Vsig, label = 'Signal')
Voutaline = plt.plot(t, Vouta, label = 'Comparator A Output')
Voutbline = plt.plot(t, Voutb, label = 'Comparator B Output')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('Signal Switching Behaviour')
plt.legend(loc = 'best')
plt.show()

#%%

plt.figure()
Voutaline = plt.plot(t, I, label = 'Power Supply Current')
plt.xlabel('Time (ms)')
plt.ylabel('Current (mA)')
plt.title('Total Current Draw')
plt.legend(loc = 'best')
plt.show()
