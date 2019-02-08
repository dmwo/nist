#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" jdi File Modifier

Takes in a jdi file that contains information about the percent dosage to be
applied to resists in order to compensate for various production anomalies,
factors in the nominal dosage value in uC/cm^2, and sweeps a variable amount
of dosages to output to a modified jdi file that includes the percent dosages
for each sweep dosage.
"""

import numpy as np
import math

__author__ = 'Dylan Oh'

def jdiSweep(
        dose_nominal = 350,
        dose_min = 300,
        dose_max = 400,
        dose_step = 10,
        shot_pitch = 4e-9,
        current = 2e-9,
        jdi_in = 'HBT_ebeamSleevev5',
        jdi_out = 'test',
        jdf_out = 'test2',
        v30_in = 'HBT_ebeamSleevev5',
        top_left = (-4850,4850),
        bottom_right = (4850,-4850),
        size_x = 10000,
        size_y = 10000
    ):

    """ Modifies a .jdi file to sweep a range of doses.

    Keyword arguments:
    dose_nominal -- The nominal dosage (default 350)
    dose_min     -- The min sweep dosage (default 300)
    dose_max     -- The max sweep dosage (default 400)
    dose_step    -- The stepping dosage (default 10)
    shot_pitch   -- The shot pitch of the JEOL EBW in metres (default 4e-9)
    current      -- The current of the JEOL EBW in amps (default 2e-9)
    jdi_in       -- The name of the input file (no extension)
    jdi_out      -- The name of the output file (no extension)
    """

    #-------------------------------------------------------------------------#
    #                       Extracting Data from File                         #
    #-------------------------------------------------------------------------#
    f = open('%s.jdi' % jdi_in, 'r')
    jdi_data = f.read()
    f.close()

    #-------------------------------------------------------------------------#
    #                   Cleaning File and Extracting Values                   #
    #-------------------------------------------------------------------------#
    # Removing unnecessary characters
    char_delete = [' ', '\n-', '\n', ')', '(']
    for i in char_delete:
        jdi_data = jdi_data.replace(i, '')

    # Splitting string by ';', which separates data from trailer information
    temp = jdi_data.split(';')

    # Splitting data string by ',', which yields individual values
    temp = temp[0].split(',')
    list_data = temp[1:len(temp)]

    # Pulling relevant numbers from data string and converting to float
    vals = [None] * (int((len(list_data) / 2) + 1))
    for i in range(0, len(list_data), 2):
        vals[int(i / 2)] = float(list_data[i])

    # Calculating minimum dose allowable
    abs_min = (current / (shot_pitch ** 2 * 50e6) / 350) - 100

    #-------------------------------------------------------------------------#
    #       Sweeping Doses and Creating Matrix of Adjusted Dose Values        #
    #-------------------------------------------------------------------------#    
    # Initialising sweep variables
    dose_adj = np.zeros((dose_step, len(vals)))
    sweep = np.linspace(dose_min, dose_max, dose_step)

    for i in range(0, dose_step):
        for j in range(0, len(vals)):
            # Calculating the adjusted dose relative to nominal
            temp = ((100 + vals[j]) * sweep[i] / dose_nominal) - 100
            if temp < abs_min:
                # Calculating the lowest possible sweep dose under conditions
                min_sweep = ((abs_min + 100) * dose_nominal)/(100 + vals[j])
                raise ValueError(
                        """Adjusted dose falls below minimum allowable under 
                           specified conditions. Change minimum dose sweep to
                           %s.""" % min_sweep)
            else:
                dose_adj[i,j] = temp

    #-------------------------------------------------------------------------#
    #                     Calculations for the .jdi File                      #
    #-------------------------------------------------------------------------#
    # Calculating the chip mark positions
    M = []
    M.append((top_left[0],top_left[1]))
    M.append((bottom_right[0],top_left[1]))
    M.append((bottom_right[0],bottom_right[1]))
    M.append((top_left[0],bottom_right[1]))
    
    # Calculating the optimal number of cells in the x and y directions
    n = math.floor(math.sqrt(dose_step))
    d = math.sqrt(dose_step) - n
    if math.floor(n) % 2 == 0:
        x = n + 1
        y = x
    else:
        if d != 0:
            y = n
            x = y + 2
        else:
            x = n
            y = n
        
    #-------------------------------------------------------------------------#
    #                          Writing New .jdi File                          #
    #-------------------------------------------------------------------------#
    # Remainder determines whether the last line will be complete (3 data
    # points) or not
    remainder = len(vals) % 3

    # max_val evaluates what the last number in the for loop will be, as it
    # iterates by 3
    max_val = int(len(vals) / 3) * 3 - 1

    f = open('%s.jdi' % jdi_out, 'w')
    g = open('%s.jdf' % jdf_out, 'w')

    MOD_vec = []

    # Iterating through MOD types
    for i in range(0, dose_step):
        # Calculating number of zeros needed to make MOD number 3 digits
        num_zeros = 3 - len(str(i + 1))
        MOD_vec.append(('0' * num_zeros) + str(i + 1))
        # Writing header
        f.write('MOD%s: MODULAT (' % MOD_vec[i])
        # Finishing first line
        f.write('( %s, %s ) , ( %s, %s ) , ( %s, %s )\n' 
                % (str(0), '%.1f' % dose_adj[i,0],
                   str(1), '%.1f' % dose_adj[i,1],
                   str(2), '%.1f' % dose_adj[i,2]))
        # Iterating through data points in groups of 3 (3 points per line)
        for j in range(3, len(vals), 3):
            # Compensating for incomplete last line (only one data point)
            if remainder == 1 and max_val <= j:
                f.write('-     , ( %s, %s ))\n' 
                        % (str(j), '%.1f' % dose_adj[i,j]))
            # Compensating for incomplete last line (only two data points)
            elif remainder == 2 and max_val <= j:
                f.write('-     , ( %s, %s ) , ( %s, %s ))\n'
                        % (str(j), '%.1f' % dose_adj[i,j],
                           str(j + 1), '%.1f' % dose_adj[i,j + 1]))
            # Adding extra closing parenthesis if last line is complete
            elif remainder == 0 and max_val == j:
                f.write('-     , ( %s, %s ) , ( %s, %s ) , ( %s , %s ))\n'
                        % (str(j), '%.1f' % dose_adj[i,j],
                           str(j + 1), '%.1f' % dose_adj[i,j + 1],
                           str(j + 2), '%.1f' % dose_adj[i,j + 2]))
            # Standard data point writing
            else:
                f.write('-     , ( %s, %s ) , ( %s, %s ) , ( %s , %s )\n'
                        % (str(j), '%.1f' % dose_adj[i,j],
                           str(j + 1), '%.1f' % dose_adj[i,j + 1],
                           str(j + 2), '%.1f' % dose_adj[i,j + 2]))

    f.close()
    
    #-------------------------------------------------------------------------#
    #                          Writing New .jdf File                          #
    #-------------------------------------------------------------------------#
    g.write('JOB/W   ,3,3\n\n')
    g.write('GLMPOS P=(-32500,2000),Q=(32500,3000)\n')
    g.write('GLMP 3.0,1000.0,0,0\n')
    g.write('GLMQRS 3.0,1000.0,0,0\n\n')
    g.write('PATH DIRE01\n\n')
    
    # Writing size of the chip
    g.write('1: ARRAY   (-%s,%s,%s)/(%s,%s,%s)\n\n' 
            % (str(size_x * math.floor(x / 2)), str(x), str(size_x),
               str(size_y * math.floor(y / 2)), str(y), str(size_y)))
    
    # Writing positions of the chip mark
    g.write('CHMPOS  M1=%s, M2=%s, M3=%s, M4=%s\n\n' 
            % (M[0],M[1],M[2],M[3]))
    
    g.write('CHMARK  3.0,50.0,0,0\n\n')
    
    # Writing indices for assign command
    count = 0
    for j in range(1, y + 1):
        for i in range(1, x + 1):
            if (i,j) == (round((x + 1) / 2), round((y + 1) / 2)):
                pass
            else:
                if count == dose_step:
                    pass
                else:
                    g.write('ASSIGN P(1) -> ((%s,%s),MOD%s)\n' 
                            % (str(i), str(j), MOD_vec[count]))
                    count += 1
        g.write('\n')
            
    g.write('AEND\n\nPEND\n\n;--- layer definitions ---;\n\nLAYER 1\n\n')
    
    # Inserting .v30 file name
    g.write("P(1) '%s.v30'\n\n" % v30_in)
    
    # Specifying shot pitch
    g.write('SHOT A,%s\n' % str(int(shot_pitch * 10 ** 9)))
    g.write("RESTYP POSI,'ZEP-520A'\n")
    
    # Specifying nominal dose
    g.write('RESIST %s,10000\n' % str(dose_nominal))
    
    # Specifying current
    g.write('STDCUR %s\n\n' % str(int(current * 10 ** 9)))
    g.write('SHOT1:  	   MODULAT ((0,0));\n')
    
    # Inserting .jdi file name
    g.write("@ '%s.jdi'\n\n" % jdi_out)
    g.write('END')
    
    g.close()
#%%

