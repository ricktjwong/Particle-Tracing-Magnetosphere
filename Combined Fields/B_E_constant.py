#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 22:10:40 2017

@author: ricktjwong
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as spi
from mpl_toolkits.mplot3d import Axes3D

""" Global quantities """

mu = 4*np.pi*1e-7
q = 1.6 * 1e-19
mp = 1.67 * 1e-27
r0 = np.array([0.0, 0.0, 0.0])

""" Initial conditions"""
B = [0, 0, 100]
E = np.array([0, 0, 100])
xinit = [-10, 0.0, 0.0, 100.0, 0.0, 0.0]

def deriv(x,t):
    xx, xy, xz = x[0], x[1], x[2]   # Initial conditions position
    vx, vy, vz = x[3], x[4], x[5]   # Initial conditions velocity
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    a = q * (np.cross(v,B) + E) / mp
    return (vx, vy, vz, a[0], a[1], a[2])

binit = np.linalg.norm(B)
r = mp*xinit[3]/(q*binit)   # Larmar radius
T = 2*np.pi*r/xinit[3]      # Gyroperiod particle drift
print r
print T

t = np.linspace(0,2,1000)*10*T

soln = spi.odeint(deriv,xinit,t)    # Solve ODE

x, y, z = soln[:,0], soln[:,1], soln[:,2]
vx, vy, vz = soln[:,3], soln[:,4], soln[:,5]

plt.figure(1)
plt.plot(x,y)
#plt.xlim([-600-3.04e7,-3.04e7])
#plt.ylim([0,600])
plt.xlabel("position, x")
plt.ylabel("position, y")

plt.figure(2)
plt.plot(x,z)
plt.xlabel("position, x")
plt.ylabel("position, z")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
plt.show()