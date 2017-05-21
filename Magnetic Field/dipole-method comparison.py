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
mC = np.array([0.0, 0.0, 10000]) 
m = 1e4
q = 1.6e-19
mp = 1.67e-27
mu = 4*np.pi*(10**-7)
k = mu*q*m/(4*np.pi)

## T = 2.15 for m = 0,0,1e4

t = np.linspace(0,4.0,1000)
r0 = np.array([0.0, 0.0, 0.0])

def b_field(mC, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(mC, r_diff))/r_mag**5
    term2 = mC/r_mag**3
    B = (mu/(4*np.pi))*(term1 - term2) 
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
    
    xacc = (k/(mp*pos**5)) * ((2*z**2 - x**2 - y**2)*vy - 3*y*z*vz)
    yacc = (-k/(mp*pos**5)) * ((2*z**2 - x**2 - y**2)*vx - 3*x*z*vz)
    zacc = (-3*k/(mp*pos**5)) * z*(x*vy - y*vx)
    
    return [vx,vy,vz,xacc,yacc,zacc]

init = [-10,0,0,100,0,0]
initC = [-10,0,0,100,0,0]
solnC = spi.odeint(derivC,initC,t)
soln = spi.odeint(deriv,init,t)

xtrajC = []
ytrajC = []
ztrajC = []
vxtrajC = []
vytrajC = []
vztrajC = []
for i in solnC:
    xtrajC.append(i[0])

for i in solnC:
    ytrajC.append(i[1])
    
for i in solnC:
    ztrajC.append(i[2])
    
for i in solnC:
    vxtrajC.append(i[3])

for i in solnC:
    vytrajC.append(i[4])
    
for i in solnC:
    vztrajC.append(i[5])    

xtraj=[]
ytraj=[]
ztraj=[]
vxtraj=[]
vytraj=[]
vztraj=[]
for i in soln:
    xtraj.append(i[0])

for i in soln:
    ytraj.append(i[1])

for i in soln:
    ztraj.append(i[2])
    
for i in soln:
    vxtraj.append(i[3])

for i in soln:
    vytraj.append(i[4])

for i in soln:
    vztraj.append(i[5])

plt.figure(1)
plt.plot(xtrajC,ytrajC,'b')
plt.plot(xtraj,ytraj,'r')
plt.xlim([-12,12])
plt.ylim([-12,12])
plt.xlabel("x / m")
plt.ylabel("y / m")
plt.show()

## Vc is a measure of speed from paper by E H Avrett
VcC = q*mu*np.linalg.norm(mC)/(mp*np.abs(initC[0])**2)
Vc = q*mu*m/(mp*np.abs(init[0])**2)

print "for purely computational method Vc is " + str(VcC)
print "for Langrangian DE Vc is " + str(Vc)

#plt.figure(2)
#plt.plot(vxtrajC,vytrajC,'b')
#plt.show()
#
#plt.figure(3)
#plt.plot(t,vxtrajC,'y')
#plt.xlim([0,0.5])
#plt.show()
#
#plt.figure(4)
#plt.plot(vxtraj,vytraj, 'r')
#plt.show()
#
#plt.figure(5)
#plt.plot(t,vxtraj,'y')
#plt.xlim([0,0.5])
#plt.show()

plt.figure(6)
plt.plot(xtraj,ztraj,'r')
#plt.plot(xtrajC,ztrajC, 'o')
plt.xlabel("x / m")
plt.ylabel("z / m")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xtraj,ytraj,ztraj)
plt.show()

