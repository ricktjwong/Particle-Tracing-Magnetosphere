# -*- coding: utf-8 -*-
"""
Spyder Editor

First Order DE, y'+y=x, y(0)=1

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

plt.figure(1)
x=np.linspace(0, 10, 1000)
y=x-1+2*np.exp(-x)
plt.title("The Truth")
plt.plot(x,y)

# ode solver

def ode_Solver():
    x=0
    y=1
    steps=1000
    count=0
    dx=0.01
    x_values=[0]
    y_values=[1]
    while count<steps:
        x=x+dx
        y=y+(x-y)*dx
        x_values.append(x)
        y_values.append(y)
        count=count+1
    return x_values, y_values

x_values, y_values =ode_Solver()
print len(x)
plt.figure(2)
plt.title("Ode Solver")
plt.plot(x_values, y_values)

"""Ode int"""

def dy_dx(y,x):
    return x-y

xs = np.linspace(0,10,1000)
y0 = 1.0
ys = odeint(dy_dx, y0, xs)

y_output = ys[:,0]
print y_output

plt.figure(3)
plt.plot(xs, y_output)
plt.title("Scipy's odeint")
plt.xlabel("x")
plt.ylabel("y")


""" Deviation from the truth """
def truth_extractor(y,y_truth):
    dev_y=[]
    for y1,y2 in zip(y,y_truth):
        dev_y.append(np.abs(y1-y2))
    return dev_y

ode_solver_y = truth_extractor(y_values,y)
plt.figure(4)
plt.semilogy(x,ode_solver_y)

y_output_dev = truth_extractor(y_output,y)
plt.figure(5)
plt.semilogy(x,y_output_dev)
    