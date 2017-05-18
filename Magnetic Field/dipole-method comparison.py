# -*- coding: utf-8 -*-
"""
Created on Wed May 17 00:29:06 2017

@author: Danielsrq
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as spi
from mpl_toolkits.mplot3d import Axes3D

## tagged with C (for computational) to avoid confusion with analytic attempt
mC = np.array([0.0, 0.0, 1e3]) 
m = 1e4
q = 1.6e-19
mp = 1.67e-27
mu = 4*np.pi*(10**-7)
k = mu*q*m/(4*np.pi)

## T = 2.15 for m = 0,0,1e4

t = np.linspace(0,4,600)
r0 = np.array([0.0, 0.0, 0.0])

def b_field(mC, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(mC, r_diff))/r_mag**5
    term2 = mC/r_mag**3
    B = mu*(term1 - term2) 
    uniform = [0,0,5e-9]    # option for constant B field
    return B

def derivC(x,t):
    xx, xy, xz = x[0], x[1], x[2]
    vx, vy, vz = x[3], x[4], x[5]
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    B = b_field(mC, x, r0)
    a = (q * np.cross(v,B) / mp)
    return (vx, vy, vz, a[0], a[1], a[2])    
    
## Differential Eqns from Lagrangian Mechanics
def deriv(r,t):
    x,y,z = r[0] , r[1] , r[2]
    vx,vy,vz = r[3], r[4] , r[5]
    pos = np.sqrt(x**2 + y**2 + z**2)
    
    xacc = (k/(mp*pos**5)) * ((2*z**2-x**2-y**2)*vy -3*y*z*vz)
    yacc = (-k/(mp*pos**5)) * ((2*z**2-x**2-y**2)*vx -3*x*z*vz)
    zacc = (-3*k/(mp*pos**5)) * z*(x*vy - y*vx)
    
    return [vx,vy,vz,xacc,yacc,zacc]

init = [-10,0,0,100,0,0]
initC = [-10,0,0,100,0,0]
solnC = spi.odeint(derivC,initC,t)
soln = spi.odeint(deriv,init,t)

xtrajC = []
ytrajC = []
ztrajC = []
for i in solnC:
    xtrajC.append(i[0])

for i in solnC:
    ytrajC.append(i[1])
    
for i in solnC:
    ztrajC.append(i[2])

xtraj=[]
ytraj=[]
ztraj=[]
for i in soln:
    xtraj.append(i[0])

for i in soln:
    ytraj.append(i[1])

for i in soln:
    ztraj.append(i[2])

plt.figure(1)
plt.plot(xtrajC,ytrajC,'o')
plt.plot(xtraj,ytraj,'r')
plt.xlabel("x / m")
plt.ylabel("y / m")
plt.show()

#plt.figure(2)
#plt.plot(xtraj,ztraj,'r')
#plt.plot(xtrajC,ztraj, 'o')
#plt.xlabel("x / m")
#plt.ylabel("z / m")
#plt.show()
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot(xtraj,ytraj,ztraj)
#plt.show()

### semi-analytical solution of DE differs from purely computational calculation 
### by 1 order of magnitude. Purely computational solution appears to be more accurate