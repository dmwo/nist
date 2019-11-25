# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:09:45 2019

@author: dmo2
"""

def fan(c = 12.7, x = 140, y = 100, n = 10):
    import math
    a = math.pi / (n - 1)
    r = c / 2 / math.sin(a / 2)
    xvec = []; yvec = []
    for i in range(n):
        xvec.append(round((x - r * math.cos(a * i)), 3))
    for i in range(round(n / 2)):
        yvec.append(round((y - 5 - r * math.sin(a * i)), 3))
    yvec += list(reversed(yvec))
    for i in range(n):
        print('{0}, {1}\n'.format(xvec[i], yvec[i]))