#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 02:37:41 2017

@author: ricktjwong
"""


from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

q = -1.602e-19
m = 9.11e-31
E = 10.

""" 

Second Order DE, x''=qE/m
We change this to v'=qE/m

"""

def dv_dx(x,t):
    xx=x[0]
    xy=x[1]
    vx=x[2]
    vy=x[3]
    ax=0.
    ay=q*E/m
    return [vx, vy, ax, ay]

xinit = [0.0,0.0,1.e4,0.] # x0, y0, vx0, vy0
t = np.linspace(0.,0.01,1000)
soln = odeint(dv_dx,xinit,t)
print soln
vx = soln[:,0]
vy = soln[:,1]

plt.figure(1)
plt.plot(t,vx)
plt.plot(t,vy)
plt.xlabel("t")
plt.ylabel("v")
plt.title("E field")

