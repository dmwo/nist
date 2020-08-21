# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:58:33 2020

@author: Dylan
"""
#%%
import numpy as np
from numpy import pi, sin, cos, tanh, sqrt
from numpy.linalg import norm
import matplotlib.pyplot as plt
import scipy.optimize as opt

def _rotate_points(points, angle = 45, center = (0,0)):
    """ Rotates points around a centerpoint defined by ``center``.  ``points`` may be
    input as either single points [1,2] or array-like[N][2], and will return in kind
    """
    if angle == 0:
         return points
    angle = angle*pi/180
    ca = cos(angle)
    sa = sin(angle)
    sa = np.array((-sa, sa))
    c0 = np.array(center)
    if np.asarray(points).ndim == 2:
        return (points - c0) * ca + (points - c0)[:,::-1] * sa + c0
    if np.asarray(points).ndim == 1:
        return (points - c0) * ca + (points - c0)[::-1] * sa + c0

def _reflect_points(points, p1 = (0,0), p2 = (1,0)):
    """ Reflects points across the line formed by p1 and p2.  ``points`` may be
    input as either single points [1,2] or array-like[N][2], and will return in kind
    """
    # From http://math.stackexchange.com/questions/11515/point-reflection-across-a-line
    points = np.array(points); p1 = np.array(p1); p2 = np.array(p2)
    if np.asarray(points).ndim == 1:
        return 2*(p1 + (p2-p1)*np.dot((p2-p1),(points-p1))/norm(p2-p1)**2) - points
    if np.asarray(points).ndim == 2:
        return np.array([2*(p1 + (p2-p1)*np.dot((p2-p1),(p-p1))/norm(p2-p1)**2) - p for p in points])

def cumtrapz(x):
    """ Numpy-based implementation of the cumulative trapezoidal integration 
    function usually found in scipy (scipy.integrate.cumtrapz) """
    return np.cumsum((x[1:] + x[:-1]) / 2)

def pp(points):
    """ Plots a list of points """
    plt.plot(points[:,0], points[:,1], '.--')
    plt.gca().set_aspect('equal', 'box')

def arc(radius = 10, angle = 90, num_pts = 720):
    """ Produces an arc of points with `num_pts` per 360 degrees.  An extra point is
    tacked on each end to ensure that the numerical gradient is accurate """
    t = np.linspace(0, angle*np.pi/180, abs(int(num_pts*angle/360))-2)
    x = radius*np.cos(t)
    y = radius*np.sin(t)
    points = np.array((x,y)).T
    start_angle = 90*np.sign(angle)
    end_angle = start_angle + angle
    return points, start_angle, end_angle

def partial_euler_Reff(Reff = 3, a = 90, p = 0.2, num_pts = 4000, 
                       *args, **kwargs):
    """ Creates a partial Euler bend formed by two clothoid curves and one 
    normal arc.
    
    Parameters
    ----------
    Reff : float
        Total effective radius of the partial Euler bend.
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
    η = kwargs.get('η', 1)
    j = kwargs.get('j', 0)
    if a <= 0 or a > 180: 
        raise ValueError("'a' must be a float such that 0 < a ≤ 180.")

    # Overhead calculations
    a = np.radians(a)
    asp = p*a / 2
    Rp = (j*Reff*η + 1 - j) / (2*sqrt(asp))
    sp = (j*Reff*η + 1 - j) * sqrt(2*asp)
    s0 = 2*sp + Rp*a*(1 - p)
    scale = a / (2*sp*(s0 - sp))
    if p == 0:
        s0 = a * (j*Reff*η + 1 - j)
        scale = a / s0

    # Constructing s and K arrays
    s = np.linspace(0, s0, num_pts)
    K = np.zeros(num_pts)
    if p == 0: K += 1
    else:
        i1 = np.argmax(s > sp)
        i2 = np.argmax(s >= s0 - sp)
        K = np.concatenate([np.multiply(np.ones(i1), 2*s[:i1]),
                            np.multiply(np.ones(i2-i1), 2*sp),
                            np.multiply(np.ones(num_pts-i2), 
                                        2 * (s0 - s[i2:num_pts]))])
    K *= scale * ((1 - j)/Reff + j)
    s *= Reff * (1 - j) + j

    # Integrating to find x and y
    ds = s[1] - s[0]
    φ = cumtrapz(K*ds)
    x = np.cumsum(ds*cos(φ))
    y = np.cumsum(ds*sin(φ))
    x = np.concatenate([[0], x])
    y = np.concatenate([[0], y])

    # Calculating η rescaling factor
    middle = int((num_pts - 1) / 2)
    η = Reff / (y[middle] + x[middle] / np.tan(a/2))

    if j == 1: return x, y
    else: return partial_euler_Reff(Reff, np.degrees(a), p, num_pts, η = η, j = 1)

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

    FIXME add reference
    """
    # Overhead calculations
    a = np.radians(a)
    sp = sqrt(p*a)      # Clothoid-to-normal transition point s value
    s0 = 2*sp + a*(1-p)/(2*sqrt(p*a/2))
    c = 1 / (2*sp*Rmin) # Scaling factor to enforce Rmin
    print(sp)

    # Constructing s and K arrays
    s = np.linspace(0, s0, num_pts)
    if p == 0: K = np.array([[1/Rmin] * len(s)])
    else:
        i1 = np.argmax(s > sp)
        i2 = np.argmax(s >= s0 - sp)
        print(i1)
        print(i2)
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

def spiral(num_turns = 3.25, gap = 3, inner_gap = 9, num_pts_half = 3000):
    """ Creates a 

    Parameters
    ----------
    num_turns : int or float
        The number of turns in the spiral. Must be greater than or equal to 1. 
        If ``num_turns`` is an even integer, the output port of the spiral 
        will be extended to be on the same side as the input port. If 
        ``num_turns`` is an odd integer, the spiral will be rotationally 
        symmetrical. If ``num_turns`` is a float, the output port of the 
        spiral will be extended, with a 0.5 decimal value resulting in a 90 
        degree increase in the arm length.
    gap : int or float
        The distance between any point on one arm of the spiral and the adjacent point on the opposite arm of the spiral.
    inner_gap : int or float
        The inner size of the spiral, equal to the distance between the two 
        spiral-to-S-curve transition points.
    num_pts_half: int
        The number of points in one arm of the normal curve. If ``num_turns`` 
        is int, the total number of points in the spiral will be 
        2 * ``num_pts_half``; if ``num_turns`` is float, the total number of 
        points will be 2 * ``num_pts_half`` + num_pts_stub.

    Returns
    -------
    x : ndarray
        Array of the x-coordinates of the spiral.
    y : ndarray
        Array of the y-coordinates of the spiral.
    """
    if num_turns < 1: raise ValueError('"num_turns" must be greater than or '
                                       'equal to 1')
    num_turns1 = np.floor(num_turns)
    num_turns2 = num_turns1 + 2*(num_turns - num_turns1)
    if (num_turns % 2) == 0: num_turns1 -= 1
    num_pts1 = int(num_pts_half/2)

    a1 = pi/2
    a2 = pi*num_turns1 + a1
    a3 = pi*num_turns2 + a1
    n_range = (a2-a1) / (num_pts1-1)

    a = inner_gap/2 - gap/2
    b = gap/pi
    a_spiral = np.array([np.linspace(a1, a2, num_pts1),
                        np.concatenate([np.linspace(a1, a2, num_pts1),
                                        np.arange(a2, a3, n_range)[1:]])])
    r_spiral = a + b * a_spiral
    x_spiral = np.array([np.zeros(num_pts1), np.zeros(len(a_spiral[1]))])
    y_spiral = np.array([np.zeros(num_pts1), np.zeros(len(a_spiral[1]))])
    for i in range(2):
        x_spiral[i] = r_spiral[i]*cos(a_spiral[i])
        y_spiral[i] = r_spiral[i]*sin(a_spiral[i])

    a_centre = np.linspace(0, a1, num_pts1)
    m = (y_spiral[0][0] - y_spiral[0][1])/-x_spiral[0][1]
    # h = inner_gap/2
    phi = a_centre[-2]
    d = opt.brentq(
        lambda z: tanh(z*phi)/tanh(z*pi/2) - 1/(sin(phi) - m*cos(phi)),
        1e-6, 10
        )
    # d = opt.brentq(lambda z: (h*sin(phi)*tanh(z*phi)/tanh(z*pi/2) - h)/(h*cos(phi)*tanh(z*phi)/tanh(z*pi/2)) - m, 1e-6, 10)
    c = inner_gap / (2*tanh(d*pi/2))
    r_centre = c * tanh(d*a_centre)
    x_centre = r_centre*cos(a_centre); y_centre = r_centre*sin(a_centre)

    x1, y1 = np.concatenate([[x_centre, y_centre],
                            [x_spiral[0][1:], y_spiral[0][1:]]], axis = 1)
    x2, y2 = np.concatenate([[x_centre, y_centre],
                            [x_spiral[1][1:], y_spiral[1][1:]]], axis = 1)
    x, y = np.concatenate([[np.flip(-x2[1:]), np.flip(-y2[1:])],
                        [x1, y1]], axis = 1)

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
    s : ndarray
        Array of path lengths.
    K : ndarray
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
Rmin = 3; Reff = 3
# x, y = partial_euler(Rmin = Rmin, a = 90, p = 1, num_pts = 4000)
x, y = partial_euler_Reff(Reff = Reff, a = 180, p = 0.1, num_pts = 4000)
pp(np.array([x, y]).T)
plt.show()
s, K = curvature(x, y)
plt.plot(s, K)
plt.show()

# print('Target Rmin = {0}, Actual Rmin = {1}'.format(Rmin, 1/K[2000]))
print('Target Reff = {1}, Actual Reff = {0}'.format(np.sqrt((y[-1] - Reff)**2 + (x[-1])**2), Reff))

# %%
num_turns = 3.5; gap = 1; inner_gap = 6; num_pts_half = 3000
num_turns1 = np.floor(num_turns)
if (num_turns % 2) == 0: num_turns1 -= 1

# Creating angle array
# a[0] represents the angle covered by the normal curve,
# a[1] represents the angle covered by the normal + stub curve
a1 = pi*num_turns1 + pi/2
a2 = pi*num_turns + pi/2
a_centre = np.concatenate()
a = np.array([np.linspace(0, a1, num_pts_half),
                np.concatenate([np.linspace(0, a1, num_pts_half),
                                np.arange(a1,a2, a1/(num_pts_half-1))[1:]])])

# Calculating relevant indices
# i1 is the transition point between the centre arc and the spiral
# i2 is the end of the spiral section for the normal curve and the stub
i1 = np.argmax(a[0] > pi/2)
i2 = [len(x) for x in a]

# Forming radius array
r = np.array([np.ones(i2[0]), np.ones(i2[1])])
for i in range(2):
    r[i][:i1] = inner_gap/2 * sin(a[0][:i1])
    r[i][i1:i2[0]] = inner_gap/2 + (a[0][i1:i2[0]] - pi/2)/pi*gap
if i2[0] == 0 or i2[1] != 0:
    r[1][i2[0]:] = inner_gap/2 + (a[1][i2[0]:] - pi/2)/pi*gap
else: pass

s1 = (a[0]*r[0])[i1]

# Combining both arms of the spiral together and converting to Cartesian
# a, r = np.concatenate([[np.flip(a[1]), -np.flip(r[1])], [a[0], r[0]]],
#                       axis = 1)
x = r[0][i1:i2[0]] * cos(a[0][i1:i2[0]]); y = r[0][i1:i2[0]] * sin(a[0][i1:i2[0]])

#%%
x, y = spiral()
plt.plot(x, y)
plt.show()
s, K = curvature(x, y)
plt.plot(s, K)