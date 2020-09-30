# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:58:33 2020

@author: Dylan
"""
#%%
import numpy as np
from numpy import sqrt
from numpy.linalg import norm
import matplotlib.pyplot as plt
import scipy.optimize as opt

def _rotate_points(points, angle = 45, center = (0,0)):
    """ Rotates points around a centerpoint defined by ``center``.  ``points`` may be
    input as either single points [1,2] or array-like[N][2], and will return in kind
    """
    if angle == 0:
         return points
    angle = angle*np.pi/180
    ca = np.cos(angle)
    sa = np.sin(angle)
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
    """ Produces an arc of points with `num_pts` per 360 degrees.  An extra 
    point is tacked on each end to ensure that the numerical gradient is 
    accurate """
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
    Rp = (j*Reff*η + 1 - j) / (2*np.sqrt(asp))
    sp = (j*Reff*η + 1 - j) * np.sqrt(2*asp)
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
    x = np.cumsum(ds*np.cos(φ))
    y = np.cumsum(ds*np.sin(φ))
    x = np.concatenate([[0], x])
    y = np.concatenate([[0], y])

    # Calculating η rescaling factor
    middle = int((num_pts - 1) / 2)
    η = Reff / (y[middle] + x[middle] / np.tan(a/2))

    if j == 1: return x, y
    else: return partial_euler_Reff(Reff, np.degrees(a), p, num_pts, η = η, j = 1)


def partial_euler(angle = 90, Rmin = 3, Reff = None, p = 0.2, num_pts = 720):
    """ Creates a partial Euler bend formed by two clothoid curves and one 
    normal arc.
    
    Parameters
    ----------
    a : float
        Total angle of the partial Euler bend in degrees.
    Rmin : float
        Radius of the normal portion of the Euler bend.
    Reff : float
        Effective radius of the total Euler bend.
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
    num_pts = abs(int(num_pts * angle/360))
    angle = np.radians(angle)
    sp = np.sqrt(p*angle)      # Clothoid-to-normal transition point s value
    s0 = 2*sp + angle*(1-p)/(2*np.sqrt(p*angle/2))
    c = 1 / (2*sp*Rmin) # Scaling factor to enforce Rmin
    print(sp)

    # Constructing s and K arrays
    s = np.linspace(0, s0, num_pts)
    if p == 0: K = np.array([[1/Rmin] * len(s)])
    else:
        i1 = np.argmax(s > sp)
        i2 = np.argmax(s >= s0 - sp)
        K = c * np.concatenate([np.multiply(np.ones(i1), 2*s[:i1]),
                                np.multiply(np.ones(i2-i1), 2*sp),
                                np.multiply(np.ones(num_pts-i2), 
                                            2*(s0 - s[i2:num_pts]))])

    # Integrating to find x and y
    ds = s[1] - s[0]
    phi = cumtrapz(K*ds)
    x, y = np.concatenate([np.array([[0],[0]]), 
                           np.cumsum([ds*np.cos(phi), ds*np.sin(phi)], axis = 1)],
                          axis = 1)

    return x, y


def spiral(num_turns = 5, gap = 1, inner_gap = 2, num_pts = 10000):
    """ Creates a spiral geometry consisting of two oddly-symmetric 
    semi-circular arcs in the centre and two Archimedean spiral arms extending 
    outward from the ends of both arcs.

    Parameters
    ----------
    num_turns : int or float
        The number of turns in the spiral. Must be greater than 1. A full 
        spiral rotation counts as 1 turn, and the center arcs will together 
        always be 0.5 turn.
    gap : int or float
        The distance between any point on one arm of the spiral and a point 
        with the same angular coordinate on an adjacent arm.
    inner_gap : int or float
        The inner size of the spiral, equal to twice the chord length of the 
        centre arcs.
    num_pts: int
        The number of points in the entire spiral. The actual number of points 
        will be slightly different than the specified value, as they are 
        dynamically allocated using the path lengths of the spiral.

    Returns
    -------
    x : ndarray
        Array of the x-coordinates of the spiral.
    y : ndarray
        Array of the y-coordinates of the spiral.

    Notes
    -----
    ``num_turns`` usage (x is any whole number):
        - ``num_turns = x.0``: Output arm will be extended 0.5 turn to be on 
        the same side as the input.
        - ``num_turns < x.5``: Input arm will be extended by the fractional 
        amount.
        - ``num_turns = x.5``: Both arms will be the same length and the input 
        and output will be on opposite sides.
        - ``num_turns > x.5``: Output arm will be extended by the fractional 
        amount.
    """
    # Establishing number of turns in each arm
    if num_turns <= 1: raise ValueError('num_turns must be greater than 1')
    diff = num_turns - np.floor(num_turns)
    if diff < 0.5: num_turns1 = np.floor(num_turns) - 1 + 2*diff
    else: num_turns1 = np.floor(num_turns)
    if diff > 0.5: num_turns2 = np.floor(num_turns) - 1 + 2*diff
    else: num_turns2 = np.floor(num_turns)

    # Establishing relevant angles and spiral/centre arc parameters
    a1 = np.pi/2
    a2 = np.array([np.pi*num_turns1 + a1, np.pi*num_turns2 + a1])
    a = inner_gap/2 - gap/2
    b = gap/np.pi
    Rc = inner_gap*np.sqrt(1 + (b/(a+b*a1))**2) / 4
    theta = np.degrees(2*np.arcsin(inner_gap/4/Rc))

    # Establishing number of points in each arm
    s_centre = Rc*np.radians(theta)
    s_spiral = ((a + a2*b)**2 + b**2)**(3/2) / (3*(a*b + (a2*b**2)))
    z = num_pts / (s_spiral[0] + s_spiral[1] + 2*s_centre)
    num_pts0 = int(z*s_centre)
    num_pts1 = int(z*s_spiral[0])
    num_pts2 = int(z*s_spiral[1]) - num_pts1

    # Forming both spiral arms
    arm1 = np.linspace(a1, a2[0], num_pts1)
    arm2 = np.linspace(a2[0], a2[1], num_pts2)[1:]
    a_spiral = np.array([arm1, np.concatenate([arm1, arm2])])
    r_spiral = a + b*a_spiral
    x_spiral = np.array([np.zeros(num_pts1), np.zeros(len(a_spiral[1]))])
    y_spiral = np.array([np.zeros(num_pts1), np.zeros(len(a_spiral[1]))])
    for i in range(2):
        x_spiral[i] = r_spiral[i]*np.cos(a_spiral[i])
        y_spiral[i] = r_spiral[i]*np.sin(a_spiral[i])

    # Forming centre arcs
    pts = _rotate_points(arc(Rc, theta, 360*num_pts0/theta)[0], -theta/2)
    x_centre = pts[:,0] + x_spiral[0][0] - pts[:,0][-1]
    y_centre = pts[:,1] + y_spiral[0][0] - pts[:,1][-1]
    x_centre = np.concatenate([-np.flip(x_centre), x_centre])
    y_centre = np.concatenate([-np.flip(y_centre), y_centre])

    # Combining into final spiral
    x = np.concatenate([-np.flip(x_spiral[1]), x_centre, x_spiral[0]])
    y = np.concatenate([-np.flip(y_spiral[1]), y_centre, y_spiral[0]])
    return x, y


def curvature(pts):
    x = pts[:,0]
    y = pts[:,1]
    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt((dx)**2 + (dy)**2)
    s = np.cumsum(ds)
    theta = np.arctan2(dy,dx)
    # Fix discontinuities arising from np.arctan2
    dtheta = np.diff(theta)
    dtheta[np.where(dtheta > np.pi)] += -2*np.pi
    dtheta[np.where(dtheta < -np.pi)] += 2*np.pi
    theta = np.concatenate([[0], np.cumsum(dtheta)]) + theta[0]
    K = np.gradient(theta, s, edge_order = 2)
    return s, K

#%%
# test_angles = [1.125, 1.25, 1.5, 1.625, 1.75, 1.875, 2]
# test_angles = [2.125, 2.25, 2.5, 2.625, 2.75, 2.875, 3]
test_angles = [3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]
# test_angles = [5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
# test_angles = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
# test_angles = [100, 200, 300, 400, 500, 600, 700, 800]
# test_angles = [1000, 10000, 100000, 1000000, 10000000]
# test_pts = [1000, 5000, 10000, 50000, 100000, 500000]
# test_pts = [10000]
for turns in test_angles:
    x, y = spiral(gap = 0.1, inner_gap = 1, num_turns = turns, num_pts = 50000)
    print('Number of turns: %s' % str(turns))
    plt.plot(x, y)
    plt.show()

# for pts in test_pts:
#     x, y = spiral(num_turns = 4.5, num_pts = pts)
#     print('Number of points in half the spiral: %s' % str(pts))
#     plt.plot(x, y)
#     plt.show()

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