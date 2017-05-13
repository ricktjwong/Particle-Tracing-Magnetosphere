#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:06:40 2017

@author: ricktjwong
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as spi

""" Global quantities """

mu = 1e-7
q = -1.6 * 1e-19
me = 9.11 * 10e-31
m = np.array([0.0, 1.0, 0.0]) 
r0 = np.array([0.0, 0.0, 0.0])

def b_field(m, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(m, r_diff))/r_mag**5
    term2 = m/r_mag**3
    B = mu*(term1 - term2) 
    return B

def deriv(x,t):
    xx, xy, xz = x[0], x[1], x[2]
    vx, vy, vz = x[3], x[4], x[5]
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    B = b_field(m, x, r0)
    a = (q * np.cross(v,B) / me)
    return (vx, vy, vz, a[0], a[1], a[2])
    
xinit = [0.0,1.0,0.0,200.,0.0,0.0]
t = np.linspace(0.,0.01,1000)

soln = spi.odeint(deriv,xinit,t)
print np.shape(soln)

x = soln[:,0]  
y = soln[:,1]   
z = soln[:,2]
vx = soln[:,3]  
vy = soln[:,4] 
vz = soln[:,5]

plt.figure(1,)
plt.plot(x,y)
plt.xlabel("position, x")
plt.ylabel("position, y")

plt.figure(2,)
plt.plot(x,z)
plt.xlabel("position, x")
plt.ylabel("position, z")