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

def rec_spherical(pos_vel):
    xx, xy, xz = pos_vel[0], pos_vel[1], pos_vel[2]
    vx, vy, vz = pos_vel[3], pos_vel[4], pos_vel[5]
    r_ = np.array([xx,xy,xz])
    r = np.linalg.norm(r_)
    rho = np.sqrt(r_[0]**2 + r_[1]**2)
    v_ = np.array([vx,vy,vz])
    
    r_hat = np.array([r_[0] , r_[1], r_[2]]) * (1/r)
    theta_hat = np.array([r_[2]*r_[0], r_[2]*r_[1], -rho**2]) * (1/(r*rho))
    phi_hat = np.array([-r_[1], r_[0], 0] ) * (1/rho)
    
    vr = np.dot(v_, r_hat)
    vtheta = np.dot(v_, theta_hat)
    vphi = np.dot(v_, phi_hat)     
    return (vr, vtheta, vphi)

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
    
mod_r = []    
for i in posr:
    mod_r.append(np.linalg.norm(i))

TE = []
for i,j in zip(KE, PE):
    TE.append(i+j)

plt.figure(4)
TE, = plt.plot(t, TE, linestyle='--', color='black')
PE, = plt.plot(t, PE, color='red')
KE, = plt.plot(t, KE, color='blue')
plt.legend([TE, PE, KE], ['TE', 'PE', 'KE'])
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")


""" Trajectories """

plt.figure(5)
xpos, = plt.plot(t, x, color='red')
ypos, = plt.plot(t, y, color='blue')
zpos, = plt.plot(t, z, color='green')
r, = plt.plot(t, mod_r, linestyle='--', color='black')
plt.legend([xpos, ypos, zpos, r], ['xpos', 'ypos', 'zpos', 'r'])
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")

""" Velocities in spherical coordinates """

vrtp = []
for i in soln:
    vrtp.append(rec_spherical(i))

vradial = []
vtheta = []
vphi = []
vtotal = []

for i in vrtp:
    vradial.append(i[0])
    vtheta.append(i[1])
    vphi.append(i[2])
    vtotal.append((i[0]**2+i[1]**2+i[2]**2)**0.5)
    
plt.figure(6)
radial_v, = plt.plot(t, vradial, color='red')
theta_v, = plt.plot(t, vtheta, color='blue')
phi_v, = plt.plot(t, vphi, color='green')
v_total, = plt.plot(t, vtotal, linestyle='--', color='black')
plt.legend([xpos, ypos, zpos, r], ['radial_v', 'theta_v', 'phi_v', 'total speed'])
plt.xlabel("Time (s)")
plt.ylabel("Velocity (ms-1)")