# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 17:42:13 2019

@author: dmo2
"""

import random
import numpy as np

#%%
def small_world(p = 0.5, n = 10, k = 2):
    # constructing matrix
    mat = [0] * n
    for i in range(n):
        temp = [0] * n
        # connect each member with nearest to kth neighbours
        for j in range(1, k+1):
            temp[i-j] = 1
            temp[i+j-n] = 1
            
        mat[i] = temp
    mat = np.matrix(mat)
        
    # rearranging matrix randomly
    for i in range(k):
        for j in range(n):
            # roll figurative dice
            if random.random() < p:
                # break connections to kth clockwise neighbour
                mat[j, j+i-n] = 0
                mat[j+i-n, j] = 0
                # make a new random connection
                r = random.randint(0, n-1)
                mat[j, r] = 1
                mat[r, j] = 1
    return mat