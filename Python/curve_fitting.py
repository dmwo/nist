# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:58:33 2020

@author: Dylan
"""
#%%
import numpy as np
from numpy import pi, sin, cos, sqrt
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

# def partial_euler(Reff = 3, a = 90, p = 0.2, num_pts = 4000, *args, **kwargs):
#     """ Creates a partial Euler bend formed by two clothoid curves and one 
#     normal arc.
    
#     Parameters
#     ----------
#     Reff : float
#         Total effective radius of the partial Euler bend.
#     a : float
#         Total angle of the partial Euler bend in degrees.
#     p : float
#         Bend parameter. Expressed as the decimal percentage of the curve that 
#         is clothoid.
#     num_pts : int
#         The number of points in the final curve.

#     Returns
#     -------
#     x : ndarray
#         Array of the x-coordinate values of the partial Euler curve.
#     y : ndarray
#         Array of the y-coordinate values of the partial Euler curve.
#     """
#     η = kwargs.get('η', 1)
#     j = kwargs.get('j', 0)
#     if a <= 0 or a > 180: 
#         raise ValueError("'a' must be a float such that 0 < a ≤ 180.")

#     # Overhead calculations
#     a = np.radians(a)
#     asp = p*a / 2
#     Rp = (j*Reff*η + 1 - j) / (2*sqrt(asp))
#     sp = (j*Reff*η + 1 - j) * sqrt(2*asp)
#     s0 = 2*sp + Rp*a*(1 - p)
#     scale = a / (2*sp*(s0 - sp))
#     if p == 0:
#         s0 = a * (j*Reff*η + 1 - j)
#         scale = a / s0

#     # Constructing s and K arrays
#     s = np.linspace(0, s0, num_pts)
#     K = np.zeros(num_pts)
#     if p == 0: K += 1
#     else:
#         i1 = np.argmax(s > sp)
#         i2 = np.argmax(s >= s0 - sp)
#         K = np.concatenate([np.multiply(np.ones(i1), 2*s[:i1]),
#                             np.multiply(np.ones(i2-i1), 2*sp),
#                             np.multiply(np.ones(num_pts-i2), 
#                                         2 * (s0 - s[i2:num_pts]))])
#     K *= scale * ((1 - j)/Reff + j)
#     s *= Reff * (1 - j) + j

#     # Integrating to find x and y
#     ds = s[1] - s[0]
#     φ = cumtrapz(K*ds)
#     x = np.cumsum(ds*cos(φ))
#     y = np.cumsum(ds*sin(φ))
#     x = np.concatenate([[0], x])
#     y = np.concatenate([[0], y])

#     # Calculating η rescaling factor
#     middle = int((num_pts - 1) / 2)
#     η = Reff / (y[middle] + x[middle] / np.tan(a/2))

#     if j == 1: return x, y
#     else: return partial_euler(Reff, np.degrees(a), p, num_pts, η = η, j = 1)

def partial_euler(Rmin = 3, a = 90, p = 0.2, num_pts = 4000):
    """ Creates a partial Euler bend formed by two clothoid curves and one 
    normal arc.
    
    Parameters
    ----------
    Rmin : float
        Radius of the normal portion of the Euler bend.
    a : float
        Total angle of the partial Euler bend in degrees.
    p : float
        Bend parameter. Expressed as the decimal percentage of the curve that 
        is clothoid.
    num_pts : int
        The number of points in the final curve.

    Returns
    -------
    x : ndarray
        Array of the x-coordinate values of the partial Euler curve.
    y : ndarray
        Array of the y-coordinate values of the partial Euler curve.
    """
    # Overhead calculations
    a = np.radians(a)
    sp = sqrt(p*a)      # Clothoid-to-normal transition point s value
    s0 = a*Rmin + sp    # Total path length derived from curvature integral = a
    c = 1 / (2*sp*Rmin) # Scaling factor to enforce Rmin

    # Constructing s and K arrays
    s = np.linspace(0, s0, num_pts)
    K = np.zeros(num_pts)
    if p == 0: K += 1/Rmin
    else:
        i1 = np.argmax(s > sp)
        i2 = np.argmax(s >= s0 - sp)
        K = c * np.concatenate([np.multiply(np.ones(i1), 2*s[:i1]),
                                np.multiply(np.ones(i2-i1), 2*sp),
                                np.multiply(np.ones(num_pts-i2), 
                                            2*(s0 - s[i2:num_pts]))])

    # Integrating to find x and y
    ds = s[1] - s[0]
    p = cumtrapz(K*ds)
    x, y = np.concatenate([np.array([[0],[0]]), 
                           np.cumsum([ds*cos(p), ds*sin(p)], axis = 1)],
                          axis = 1)

    return x, y

# def euler(s, K, num_pts):



def spiral(num_turns = 3.25, gap = 3, inner_gap = 9, num_pts_half = 3000):
    num_turns1 = np.floor(num_turns)
    # num_pts2 = 2*(num_turns - num_turns1)*num_pts_half / num_turns
    K_start = 2 / (inner_gap + 2*gap*num_turns1)
    K_end = 2 / (inner_gap + 2*gap*(num_turns1 - 1))
    K_arc = 4 / inner_gap
    a = pi * num_turns1
    s01 = 2*a / (K_start + K_end)
    s02 = pi / K_arc
    s = np.linspace(0, s01 + s02, num_pts_half)
    ds = s[1] - s[0]
    K = np.concatenate([(K_end - K_start)/s01*s[:np.argmax(s > s01)] + K_start,
                        [K_arc] * len(s[np.argmax(s >= s01):])])
    
    p = cumtrapz(K*ds)
    x, y = np.concatenate([np.array([[0],[0]]), 
                           np.cumsum([ds*cos(p), ds*sin(p)], axis = 1)],
                          axis = 1)
    # x, y = np.concatenate([(x, y), (-x + gap/2, -y + y.max())], axis = 1)
    return x, y

def curvature(x, y):
    """ Calculates the curvature vs path length for a curve defined by 
    coordinates (``x``, ``y``).

    Parameters
    ----------
    x : nd.array
        Array of x-coordinates of the input curve.
    y : nd.array
        Array of y-coordinates of the input curve.

    Returns
    -------
    s : np.array
        Array of path lengths.
    K : np.array
        Array of curvatures.
    """
    dxdt = np.gradient(x, edge_order = 1)
    dydt = np.gradient(y, edge_order = 1)
    d2xdt2 = np.gradient(dxdt, edge_order = 2)
    d2ydt2 = np.gradient(dydt, edge_order = 2)
    dφdt = (dxdt*d2ydt2 - dydt*d2xdt2) / (dxdt**2 + dydt**2)
    dsdt = sqrt((dxdt)**2 + (dydt)**2)
    s = cumtrapz(dsdt)
    s = np.concatenate([[0], s])
    K = dφdt / dsdt
    return s, K

# %%
pp(np.array([partial_euler(Rmin = 3, a = 90, p = 1, num_pts = 4000)]).T)
plt.show()
plt.plot(curvature(partial_euler(Rmin, a = 90, p = 1, num_pts = 4000)))
plt.show()

print('Target Rmin = {0}, Actual Rmin = {1}'.format(Rmin, 1/K[2000]))

#%%
num_turns = 3.25; gap = 1; inner_gap = 9; num_pts_half = 5000
# num_turns must be greater than or equal to 1
num_turns1 = np.floor(num_turns)
if (num_turns % 2) == 0: num_turns1 -= 1

# Creating angle array
# a[0] represents the angle covered by the normal curve,
# a[1] represents the angle covered by the normal + stub curve
a1 = pi*num_turns1 + pi/2
a2 = pi*num_turns + pi/2
a = np.array([np.linspace(0, a1, num_pts_half),
              np.concatenate([np.linspace(0, a1, num_pts_half),
                              np.arange(a1, a2, a1/(num_pts_half-1))[1:]])])

# Calculating relevant indices
# i1 is the transition point between the centre arc and the spiral
# i2 is the end of the spiral section for the normal curve and the stub
i1 = np.argmax(a[0] > pi/2)
i2 = [len(x) for x in a]


r = np.array([np.ones(i2[0]), np.ones(i2[1])])
r[0][:i1] = np.multiply(r[0][:i1], inner_gap/2 * sin(a[0][:i1]))
r[1][:i1] = np.multiply(r[1][:i1], inner_gap/2 * sin(a[0][:i1]))
t = inner_gap/2 + (a[0][i1:i2[0]] - pi/2)/pi*gap
r[0][i1:i2[0]] = np.multiply(r[0][i1:i2[0]], t)
r[1][i1:i2[0]] = np.multiply(r[1][i1:i2[0]], t)
if i2[0] == 0 or i2[1] != 0:
    t = inner_gap/2 + (a[1][i2[0]:] - pi/2)/pi*gap
    r[1][i2[0]:] = np.multiply(r[1][i2[0]:], t)
else: pass

a, r = np.concatenate([[np.flip(a[1]), -np.flip(r[1])], [a[0], r[0]]],
                      axis = 1)
x = r * cos(a)
y = r * sin(a)

s, K = curvature(x, y)
pp(np.array([x, y]).T)
plt.show()
plt.plot(s, K)
plt.show()
# print((Kend-Kstart)/s0)
# %%
