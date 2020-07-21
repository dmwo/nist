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

def fresnel(R0, s = 6, num_pts = 1000):
    t = np.linspace(0, s, num_pts)
    dt = t[1] - t[0]
    dxdt = np.cos(t**2 / (2 * R0**2))
    dydt = np.sin(t**2 / (2 * R0**2))
    x = cumtrapz(dxdt)*dt
    y = cumtrapz(dydt)*dt
    x = np.concatenate([[0], x])
    y = np.concatenate([[0], y])
    # Fixme prepend zeros to x and y?  YES
    return x,y

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

def partial_euler(R0 = 3, p = 0.2, a = 90, num_pts = 1000):
    """ Creates a partial Euler bend formed by two clothoid curves and one normal arc 
    
    Parameters
    ----------
    R0
    p
    a
    num_pts

    Returns
    -------
    """
    a = np.radians(a)
    asp = p * a / 2
    Rp = R0 / 2 / np.sqrt(asp)
    sp = R0 * np.sqrt(2 * asp)
    s0 = 2 * sp + Rp * a * (1 - p)
    if p == 0:
        s0 = R0 * a
        Rp = R0

    # Clothoid curve
    xbend1, ybend1 = fresnel(R0, sp, 1000)

    # Normal curve
    s = np.linspace(sp, s0 / 2, 1000)
    xbend2 = Rp * np.sin(((s - sp) / Rp + asp))
    ybend2 = Rp * (1 - np.cos(((s - sp) / Rp + asp)))
    x_offset = xbend1[-1] - xbend2[0]
    y_offset = ybend1[-1] - ybend2[0]
    xbend2 += x_offset
    ybend2 += y_offset

    # Mirroring and combining into complete curve
    xbend = np.concatenate([xbend1, xbend2[1:]])
    ybend = np.concatenate([ybend1, ybend2[1:]])
    pts = list(zip(np.flip(xbend[:-1]), np.flip((2 * max(ybend) - ybend)[:-1])))
    rot_pts = _rotate_points(pts, np.degrees(a) - 180, (xbend[-1], ybend[-1]))
    xbend_flip = np.array(([a for a,b in rot_pts], [b for a,b in rot_pts])[0])
    ybend_flip = np.array(([a for a,b in rot_pts], [b for a,b in rot_pts])[1])
    x = np.concatenate([xbend, xbend_flip])
    y = np.concatenate([ybend, ybend_flip])

    return x, y

# %%
x, y = partial_euler(a=370)
pp(np.array([x, y]).T)


# %%
