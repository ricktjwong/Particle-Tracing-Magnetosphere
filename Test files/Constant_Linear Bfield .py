# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:25:00 2017

@author: Danielsrq
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

'''Global Variables / Constants '''
epsilon0 = 8.85e-12         ## permitivity constant
k = 1/(4*math.pi*epsilon0)  ## grouped electric constant
k_2=1
Q = 5e-17                       ## relatively strong electric pole
q = 1.6e-19                       ## elementary charge
mp = 1.67e-27               ## mass of proton
s = 0.01                    ## seperation from origin
E_acc = -0.1
B = 1

''' try using self-made iteration method '''
## underscore refer to vector
B_ = np.array([0,0,B])
r_ = np.array([0,0,0])
v_ = np.array([200,0,0])
delta = 1e-10 ## step size
pos_ = np.array([0,0,0])
vel_ = np.array([200,0,0])
N = 0 ## counter
steps = 400 ## number of steps

## iteration to find position
while N < steps:
    r_ = r_ + v_*delta
    v_ = v_ + (q/mp)*(np.cross(v_,B_))*delta
    pos_=np.vstack((pos_, r_))
    vel_=np.vstack((vel_, v_))
    N += 1

print pos_

plt.figure(3)
x=pos_[:,0]
y=pos_[:,1]
plt.plot(x,y, 'b')
plt.xlabel("x /m")
plt.ylabel("y /m")
plt.show()
