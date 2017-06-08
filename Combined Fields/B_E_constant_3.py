#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 22:10:40 2017

@author: ricktjwong
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as spi
from mpl_toolkits.mplot3d import Axes3D

#plt.rcParams["figure.figsize"] = (4,3)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['lines.linewidth'] = 0.8

""" Global quantities """

mu = 4*np.pi*1e-7
q = 1.6 * 1e-19
mp = 1.67 * 1e-27
r0 = np.array([0.0, 0.0, 0.0])

""" Initial conditions"""
B = [0, 0, 1e-5]
E = np.array([0.0, 1e-4, 0.0])
xinit = [-10, 0.0, 0.0, 0.0, 0.0, 0.0]
x0 = np.array([xinit[0], xinit[1], xinit[2]])
v0 = np.array([xinit[3], xinit[4], xinit[5]])

def deriv(x,t):
    xx, xy, xz = x[0], x[1], x[2]   # Initial conditions position
    vx, vy, vz = x[3], x[4], x[5]   # Initial conditions velocity
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    a = q * (np.cross(v,B) + E) / mp
    return (vx, vy, vz, a[0], a[1], a[2])

binit = np.linalg.norm(B)
r = mp*xinit[3]/(q*binit)   # Larmar radius
T = 2*np.pi*r/xinit[3]      # Gyroperiod particle drift
print r
print T

t = np.linspace(0,0.1,1000)

soln = spi.odeint(deriv,xinit,t)    # Solve ODE

x, y, z = soln[:,0], soln[:,1], soln[:,2]
vx, vy, vz = soln[:,3], soln[:,4], soln[:,5]

plt.figure(1)
plt.plot(x,y)
#plt.xlim([-600-3.04e7,-3.04e7])
#plt.ylim([0,600])
plt.xlabel("position, x")
plt.ylabel("position, y")

plt.figure(2)
plt.plot(x,z)
plt.xlabel("position, x")
plt.ylabel("position, z")

fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()

""" Energy Calculations """

v_xyz = np.vstack((vx,vy,vz)).T
v_ = []
for i in v_xyz:
    v_.append(np.linalg.norm(i))
v_ = np.array(v_)
KE = 0.5*mp*v_**2
KE_ave = np.average(KE)
KE_var = np.var(KE)

posr = np.vstack((x,y,z)).T
PE = []
for i in posr:
    PE.append(q*np.dot((x0-i), E))

TE = []
for i,j in zip(KE, PE):
    TE.append(i+j)

plt.figure(4)
TE, = plt.plot(t, TE, linestyle='--', color='black')
PE, = plt.plot(t, PE, color='red')
KE, = plt.plot(t, KE, color='blue')
plt.legend([TE, PE, KE], ['TE', 'PE', 'KE'])
plt.xlabel("Time (s)")
plt.ylabel("Energy")

""" Parallel and Perpendicular velocities """

v_para_ = []
for i in v_xyz:
    B_hat = np.array(B) / np.linalg.norm(B)
    v_para_.append( np.dot(i,B_hat) * B_hat )

v_para_mod = []
for i in v_para_:
    v_para_mod.append(np.linalg.norm(i))

v_perp_ = []
for i in v_xyz:
    B_hat = np.array(B) / np.linalg.norm(B)
    v_perp_.append(np.cross(i, B_hat))
v_perp_mod =[]
for i in v_perp_:
    v_perp_mod.append(np.linalg.norm(i))

v_perpx = []
v_perpy = []
for i in v_perp_:
    v_perpx.append(i[0])
    v_perpy.append(i[1])
    
plt.figure(5)
v_perpx, = plt.plot(x,v_perpx, color='red')
v_perpy, = plt.plot(x,v_perpy, color='blue')
plt.legend([v_perpx, v_perpy], ['v_perpx', 'v_perpy'])
plt.xlabel("x-Distance (m)")
plt.ylabel("Velocity (ms-1)")