# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:58:33 2020

@author: Dylan
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc

def cumtrapz(x):
    """ Numpy-based implementation of the cumulative trapezoidal integration 
    function usually found in scipy (scipy.integrate.cumtrapz) """
    return np.cumsum((x[1:] + x[:-1]) / 2)

def pp(points):
    """ Plots a list of points """
    plt.plot(points[:,0], points[:,1], '.--')
    plt.gca().set_aspect('equal', 'box')

def _rotate_points(points, angle = 45, center = (0,0)):
    """ Rotates points around a centerpoint defined by ``center``.  ``points`` may be
    input as either single points [1,2] or array-like[N][2], and will return in kind
    """
    if angle == 0:
         return points
    angle = np.radians(angle)
    ca = np.cos(angle)
    sa = np.sin(angle)
    sa = np.array((-sa, sa))
    c0 = np.array(center)
    if np.asarray(points).ndim == 2:
        return (points - c0) * ca + (points - c0)[:,::-1] * sa + c0
    if np.asarray(points).ndim == 1:
        return (points - c0) * ca + (points - c0)[::-1] * sa + c0

def partial_euler(Reff = 3, a = 90, p = 0.2, num_pts = 4000):
    """ Creates a partial Euler bend formed by two clothoid curves and one normal arc.
    
    Parameters
    ----------
    R0 : int or float
        Radius of the clothoid curves.
    p  : float
        Bend parameter. Expressed as the decimal percentage of the curve that is clothoid.
    a : int or float
        Total angle of the complete curve in degrees.
    num_pts : int
        The number of points of the final curve (must be divisible by 4). Actual number of points will be ``num_pts`` + 1

    Returns
    -------
    x : ndarray
        Array of the x-coordinate values of the partial Euler curve.
    y : ndarray
        Array of the y-coordinate values of the partial Euler curve.
    """
    a   = np.radians(a)
    asp = p * a / 2
    Rp  = 1 / 2 / np.sqrt(asp)
    sp  = np.sqrt(2 * asp)

    s0  = 2 * sp + Rp * a * (1 - p)
    scale = a / (2 * sp * (s0 - sp))
    if p == 0: 
        s0 = a
        scale = a / s0

    s = np.linspace(0, s0, num_pts)
    K = np.zeros(num_pts)
    if p == 0: K += 1
    else:
        for i in range(len(K)):
            if   s[i] <= sp         : K[i] = 2 * s[i]
            elif sp < s[i] < s0 - sp: K[i] = 2 * sp
            elif s0 - sp < s[i] < s0: K[i] = 2 * (s0 - s[i])
    K *= scale / Reff
    s *= Reff

    ds = s[1] - s[0]
    φ = cumtrapz(K * ds)
    x = np.cumsum(ds * np.cos(φ))
    y = np.cumsum(ds * np.sin(φ))
    x = np.concatenate([[0], x])
    y = np.concatenate([[0], y])

    return x, y, s, K

def curvature(x, y):
    dxdt = np.gradient(x, edge_order = 1); dydt = np.gradient(y, edge_order = 1)
    d2xdt2 = np.gradient(dxdt, edge_order = 2); d2ydt2 = np.gradient(dydt, edge_order = 2)
    dφdt = (dxdt * d2ydt2 - dydt * d2xdt2) / (dxdt**2 + dydt**2)
    dsdt = np.sqrt((dxdt)**2 + (dydt)**2)
    s = cumtrapz(dsdt)
    s = np.concatenate([[0], s])
    K = dφdt / dsdt
    return s, K

# %%
x, y, s, K = partial_euler(Reff = 3, a = 90, p = 0.2, num_pts = 4000)

middle = int((num_pts - 1) / 2)
η = Reff / (y[middle] + x[middle] / np.tan(a / 2))

Rp *= Reff * η
sp *= Reff * η
s0 = 2 * sp + Rp * a * (1 - p)
scale = a / (2 * sp * (s0 - sp))
s = np.linspace(0, s0, num_pts)
K = np.zeros(num_pts)
if p == 0: K += Reff
else:
    for i in range(len(K)):
        if   s[i] <= sp         : K[i] = 2 * s[i]
        elif sp < s[i] < s0 - sp: K[i] = 2 * sp
        elif s0 - sp < s[i] < s0: K[i] = 2 * (s0 - s[i])
K *= scale
ds = s[1] - s[0]
φ = cumtrapz(K * ds)
x = np.cumsum(ds * np.cos(φ))
y = np.cumsum(ds * np.sin(φ))
x = np.concatenate([[0], x])
y = np.concatenate([[0], y])

s1, K1 = curvature(x, y)

pp(np.array([x, y]).T)
plt.show()
plt.plot(s, K)
plt.plot(s1, K1)
plt.show()

Reff = 3
print('Target Reff = {1}, Actual Reff = {0}'.format(np.sqrt((y[-1] - Reff)**2 + (x[-1])**2), Reff))


 # %%

# %%
