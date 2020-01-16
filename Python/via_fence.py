# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math
x0 = 142.75
y0 = 95.5
n = 10
r = 0.8
l = 33.64
a = 45.1

#%%
def via_fence(x0, y0, n, r, l, a):
    l_v = l / (n + 1)
    dy_top = r * math.sin(math.radians(90 + a))
    dx_top = r * math.cos(math.radians(90 + a))
    dy_bot = -dy_top
    dx_bot = -dx_top
    
    print("Topside vias\n")
    for i in range(1, n + 1):
        dx = x0 + l_v * i * math.cos(math.radians(a))
        dy = y0 - l_v * i * math.sin(math.radians(a))
        print("(" + str(round(dx + dx_top, 3)) + 
              ", " + str(round(dy - dy_top, 3)) + ")\n")
    
    print("Bottomside vias\n")
    for i in range(1, n + 1):
        dx = x0 + l_v * i * math.cos(math.radians(a))
        dy = y0 - l_v * i * math.sin(math.radians(a))
        print("(" + str(round(dx + dx_bot, 3)) + 
              ", " + str(round(dy - dy_bot, 3)) + ")\n")