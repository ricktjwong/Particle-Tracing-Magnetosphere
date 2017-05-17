#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 02:55:47 2017

@author: ricktjwong
"""

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

q = -1.602e-19
m = 9.11e-31
B = 0.2

""" 

Second Order DE, x''=qvb/m
We change this to v'=qvb/m

"""

def dv_dx(x,t):
    xx,xy,xz=x[0],x[1],x[2]
    vx,vy,vz=x[3],x[4],x[5]
    ax,ay,az=0.0,q*vx*B/m,0.0
    return [vx, vy, vz, ax, ay, az]

xinit = [0.0, 0.0, 0.0, 2.e2, 0.0, 0.0] # x0, y0, z0, vx0, vy0, vz0
t = np.linspace(0., 0.01, 1000)
soln = odeint(dv_dx,xinit,t)
print soln
xx, xy, xz = soln[:,0], soln[:,1], soln[:,2]
vx, vy, vz = soln[:,3], soln[:,4], soln[:,5]

plt.figure(1)
plt.plot(t,vx)
plt.xlabel("t")
plt.ylabel("v")
plt.title("x - B field")

plt.figure(2)
plt.plot(t,vy)
plt.xlabel("t")
plt.ylabel("v")
plt.title("y - B field")

plt.figure(3)
plt.plot(t,vz)
plt.xlabel("t")
plt.ylabel("v")
plt.title("z - B field")