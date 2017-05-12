#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:06:40 2017

@author: ricktjwong
"""

import numpy as np
from matplotlib import pyplot as plt

""" Global quantities """

mu = 1e-7
q = 1.6 * 1e-19
me = 9.11 * 10e-31

def b_field(m, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(m, r_diff))/r_mag**5
    term2 = m/r_mag**3
    B = mu*(term1 - term2) 
    return B

def electron_path(m, r0):
    dt = 1e-5
    steps = 10000
    count = 0
    
    pos = np.array([0., 1.0, 0.])
    vel = np.array([200., 0., 0.])
    
    # Initial conditions
  
    x = np.array([0., 1.0, 0.])
    v = np.array([200., 0., 0.])
    
    while count<steps:
        B = b_field(m, x, r0)
        dv = (q * np.cross(v,B) / me) * dt
        dx = (q * np.cross(v,B) / me) * (dt**2)
        v = v + dv
        x = x + dx
        pos = np.vstack((pos, x))
        vel = np.vstack((vel, v))
        count += 1
        
    return pos, vel

m=np.array([0.0, 1.0, 0.0]) 
r0=np.array([0.0, 0.0, 0.0])

pos,vel = electron_path(m,r0)

totalt=10000*1e-5
t=np.linspace(0,totalt,10001)
print t

#plt.figure(1)
#plt.title("x against t")
#plt.plot(t, pos[:,0])

plt.figure(1)
plt.title("pos z against x")
plt.plot(pos[:,0], pos[:,2])

plt.figure(2)
plt.title("x against v")
plt.plot(vel[:,2], pos[:,2])
        