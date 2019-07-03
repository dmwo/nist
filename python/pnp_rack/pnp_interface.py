# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:21:23 2018

@author: Jimmy Gammell
"""
#%%
import serial
import time
import math

class PNPRack:
    def __init__(self, port):
        """ Initialisation for the PNPRack class. Sets slot location and I2C
        addresses for connected modules, establishes timeout period, and
        connects to an Arduino. Sets DAC voltage to 0.

        Parameters
        ----------
        port    : string, Arduino USB port (in the form 'COMXX')
        """

        # Declaring variables
        self.slot_DAC = (1,3) # DAC PNP slot
        self.slot_ADC = (2,3) # ADC PNP slot
        self.addr_DAC = 0x60  # DAC I2C address
        self.addr_ADC = 0x14  # ADC I2C address
        self.timeout  = 5     # Response timeout [sec]

        # Initialising Arduino as a serial object with matching baud 19200
        self.Arduino = serial.Serial(port, 19200, timeout = self.timeout)

        # Arduino returns 's' when finished initialising
        if not 's' in self.read():
            raise Exception('Failed to connect to Arduino')

        # Resetting DAC
        self.setDACV(0)
    
    def close(self):
        self.Arduino.close()

    def read(self):
        return self.Arduino.readline().decode().strip()

    def write(self, data):
        self.Arduino.write(data.encode())

    def setSlot(self, slot):
        """ Communicates with the Arduino and provides information needed to
        select the active slot of the PNP rack.

        Parameters
        ----------
        slot : tuple, pnp slot of target device (of form (row, col))
        """

        # Writing message in an Arduino-interpretable format
        msg = 's{0}'.format(''.join([chr(8-n) for n in slot]))
        self.write(msg)
        self.read()
 
    def writeBytes(self, addr, data):
        """ Communicates with the Arduino and provides information needed to
        write byte data to an I2C device.

        Parameters
        ----------
        addr : int  , target device address on I2C bus
        data : tuple, data to send to target device (truncated to bytes)
        """

        # Writing message in an Arduino-interpretable format
        msg = 'w{0}{1}\n'.format(chr(addr),
                                 ''.join([chr(b) for b in data]))
        self.write(msg)
        self.read()

    def readBytes(self, addr, num_bytes):
        """ Communicates with the Arduino and provides information needed to
        request and read a specified number of bytes from the target device.

        Parameters
        ----------
        addr      : int, target device address on I2C bus
        num_bytes : int, number of bytes to request

        Returns
        -------
        in_data : list, data returned by the device
        """

        # Requesting read data from target device
        msg = 'r{0}{1}'.format(chr(addr), chr(num_bytes))
        self.write(msg)
        
        timeout = time.time() + self.timeout
        in_data = []
        while len(in_data) < num_bytes:
            # Case where Arduino buffer is empty and timeout has occurred
            if not self.Arduino.in_waiting and timeout < time.time():
                raise Exception('Failed to read device register')
            # Case where Arduino buffer is not empty / data available
            elif self.Arduino.in_waiting:
                val = ''.join([c if ord('0') <= ord(c) <= ord('9')\
                                 else '' for c in self.read()])
                if val != '': in_data.append(int(val))
        self.read()
        return in_data

    def setDACV(self, volt, slot = None):
        """ Sets voltage of the MCP4706 8-bit DAC.

        Parameters
        ----------
        volt : float, voltage setting of the DAC (range 0-3.288[V])
        slot : tuple, pnp slot of target device (of form (row, col))
        """        

        # Check that input conditions are met, set default values otherwise
        if not slot: slot = self.slot_DAC
        if not 0 <= volt <= 3.288:
            raise Exception('Voltage must be in [0, 3.288]')

        # Representing voltage as an 8-bit integer value
        volt = int(math.floor(255 * (float(volt) / 3.2871)))
        self.writeBytes(self.addr_DAC, (0x70, volt, 0x00))

    def readDAC(self, slot = None):
        """ Requests a DAC read and retrieves both the integer value at the
        DAC register and the current voltage of the DAC. Prints both values and
        returns the voltage.

        Parameters
        ----------
        slot : tuple, pnp slot of target device (of form (row, col))
        
        Returns
        -------
        volt : float, voltage value of the DAC
        """ 
        
        # Accessing the DAC
        if not slot: slot = self.slot_DAC
        self.setSlot(slot)
        self.writeBytes(self.addr_DAC, ()) # must write before read to wake DAC
        out = []
        timeout = time.time() + self.timeout
        
        # Continue to read bytes until 2 bytes have been received or until
        # timeout occurs
        while len(out) != 2 and timeout > time.time():
            out = self.readBytes(self.addr_DAC, 2)
        # Case where 2 bytes haven't been received and timeout has occurred
        if len(out) != 2:
            raise Exception('Failed to read DAC register')
        out = out[1] # second byte contains relevant voltage data
        volt = 3.288 * out / 255.
        
        print('The integer DAC value is: {}'.format(str(out)))
        print('The DAC voltage is: {} V'.format(str(round(volt, 3))))
        # return volt

    def readADC(self, slot = None, ref = 5):
        """ Requests a read of the LTC2451 16-bit ADC and retrieves both the 
        integer value at the ADC register and the current voltage of the ADC.
        Prints both values and returns the voltage.
        
        Parameters
        ----------
        slot : tuple, pnp slot of target device (of form (row, col))
        ref  : int  , reference value
        
        Returns
        -------
        volt : float, voltage value of the ADC
        """         

        # Accessing the ADC
        if not slot: slot = self.slot_ADC
        self.setSlot(slot)
        self.writeBytes(self.addr_ADC, (0x00,)) # wake and set up device
        
        time.sleep(.05)
        # ADC reads before measuring--must read twice
        self.readBytes(self.addr_ADC, 2)
        out = self.readBytes(self.addr_ADC, 2) # returns two output bytes
        out = (out[0] << 8) | out[1] # combine bytes into single integer
        volt = float(out) * float(ref) / 65535.
        
        print('The integer ADC value is: {}'.format(str(out)))
        print('The ADC voltage is: {} V'.format(str(round(volt, 3))))
        # return volt

port = 'COM11'
try: pnp.close()
except: pass
pnp = PNPRack(port) 
