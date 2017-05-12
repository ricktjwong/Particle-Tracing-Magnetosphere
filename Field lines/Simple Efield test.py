# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:13:09 2017

@author: Danielsrq
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import scipy.integrate as spi
import scipy as sp

'''Global Variables / Constants '''
epsilon0 = 8.85e-12         ## permitivity constant
k = 1/(4*math.pi*epsilon0)  ## grouped electric constant
k_2=1
Q = 5e-17                       ## relatively strong electric pole
q_ = 1.6e-19                       ## elementary charge
mp = 1.67e-27               ## mass of proton
s = 0.01                    ## seperation from origin
E_acc = -0.1

''' Define E field function at a position r0 ; charges align on x axis '''
def E (q, r0, x, y):
    ''' Need to obtain value of E at any postition r(x,y) '''
    denom = ((x-r0[0])**2 + y**2)**1.5
    ## returns compenents in i j
    return k*q*(x-r0[0]) / denom , k*q*y / denom

## initial conditions (x,y,Vx,Vy)
delta = 1e-3
#V0 = (-0.5,1,1,1)
x=0
y=0.1
vx= 0
vy= 0
N=0
x0= []
y0= []
vx0= []
vy0= []
ax0= []
ay0= []
step = np.arange(0,100,delta)
while N < 500:
    x = x + vx*delta
    y = y + vy*delta
    vx = vx + (q_/mp)*E(-5e-17, (-0.1,0), x, y)[0]*delta + (q_/mp)*E(5e-17, (0.1,0), x, y)[0]*delta
    vy = vy + (q_/mp)*E(-5e-17, (-0.1,0), x, y)[1]*delta + (q_/mp)*E(5e-17, (0.1,0), x, y)[1]*delta
    ax0.append((q_/mp)*E(-5e-17, (-0.1,0), x, y)[0] + (q_/mp)*E(-5e-17, (0.1,0), x, y)[0])
    x0.append(x)
    y0.append(y)
    vx0.append(vx)
    vy0.append(vy)
    N += 1

### (0,0.1,0,0) ; Q = +-5e-17 , separation = 0.2 gives sickle shape plot    
plt.figure(1)
plt.plot(x0, y0, "o")
#plt.plot(t_array/200., soln_u, "b")
plt.xlabel("x")
plt.ylabel("y")

plt.figure(2)
plt.plot(x0,vx0,'b')
plt.plot(y0,vy0,'r')
plt.xlabel("position x or y")
plt.ylabel("velocity x or y")
plt.show()

'''      SPACE         '''

def E_dipole(x,y,z):
    R = (x**2 + y**2 + z**2)**0.5
    factor = (1/(k*R**5))
    
    E_f = factor*(3*x**2 - Q*0.2*R**2) , factor*y*x , factor*z*x
    return E_f
    
nx, ny = 64, 64
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)

##Bx, By = dipole(m=[0, 0.1], r=np.meshgrid(X, Y), r0=[0.0,0.0])
Ex , Ey , Ez = E_dipole(X,Y,0)
plt.figure(figsize=(5, 5))
plt.streamplot(X, Y, Ex, Ey)
#plt.margins(0, 0)
plt.show()

'''
t = np.linspace(1.E-3,0.1,500)

def deriv(r,t):
    rpos_x = r[0]
    rpos_y = r[1]
    rpos_z = r[2]
    rvel_x = r[3]
    rvel_y = r[4]
    rvel_z = r[5]
    
    R = math.sqrt(rpos_x**2 + rpos_y**2 + rpos_z**2)
    
    racc_x = (k*q_/(mp*R**5)) * (3*rpos_x**2 - Q*0.2*R**2)
    racc_y = (k*q_/(mp*R**5)) * (rpos_y*rpos_x)
    racc_z = (k*q_/(mp*R**5)) * (rpos_z*rpos_x)
    
    return [rvel_x,rvel_y,rvel_z,racc_x,racc_y,racc_z]

soln=[]
r0 = [ 0.0, 0.1, 0., 0, 0, 0]
soln.append(spi.odeint(deriv,r0,t))

pos_x = soln[0][:,0]
pos_y = soln[0][:,1]
#pos_z = soln[2]
#vel_x = soln[3]
#vel_y = soln[4]
#vel_z = soln[5]
'''
    
