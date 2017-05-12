#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:14:03 2017

@author: ricktjwong
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
#import sys  
#
#reload(sys)  
#sys.setdefaultencoding('utf8')

q = -1.602e-19
m = 9.11e-31
E = 10.

def deriv(x,t):     # Define a function which calculates the derivatives
    xx = x[0] 
    xy = x[1] 
    vx = x[2] 
    vy = x[3]
    ax = 0.
    ay = q*E/m
    return [vx,vy,ax,ay]

xinit = [0.0,0.0,1.e4,0.] # x0, y0, vx0, vy0
t = np.linspace(0.,0.01,1000)

soln = spi.odeint(deriv,xinit,t)
print soln
print type(soln)

x = soln[:,0]   # first column
y = soln[:,1]   # second column
vx = soln[:,2]  # third column
vy = soln[:,3]  # fourth column
print x
print y
print vx
print vy

print soln[2]


plt.figure(1,)
plt.plot(x,y)
plt.xlabel("position, x")
plt.ylabel("position, y")

plt.figure(2,)
plt.plot(t,vx)
plt.xlabel("time (s)")
plt.ylabel("x‐velocity (m/s)")

plt.figure(3)
plt.plot(t,vy)
plt.xlabel("time (s)")
plt.ylabel("y‐velocity (m/s)")

plt.show()