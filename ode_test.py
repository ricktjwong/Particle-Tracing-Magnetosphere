#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 01:56:02 2017

@author: ricktjwong
"""

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

""" First Order DE, y'+y=x, y(0)=1  """

def dy_dx(y,x):
    return x-y

xs = np.linspace(0,5,100)
y0 = 1.0
ys = odeint(dy_dx, y0, xs)

y_output = ys[:,0]
print y_output

plt.figure(1)
plt.plot(xs, y_output)
plt.xlabel("x")
plt.ylabel("y")

""" 

Second Order DE, y''+2y'+2y=cos(2x), y(0)=0, y'(0)=0 
We change this to z'+2z+2y=cos(2x), z(0)=y(0)=0, where z=y'

"""
def dU_dx(U,x): # Here U is a vector such that y=U[0] and z=U[1]. This function should return [y', z']
    return (U[1], -2*U[1]-2*U[0]+np.cos(2*x))

U = [0,0]
xs = np.linspace(0,10,200)
Us = odeint(dU_dx,U,xs)
ys = Us[:,0]
zs = Us[:,1]

plt.figure(2)
plt.plot(xs,ys)
plt.plot(xs,zs)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Damped harmonic oscillator")
