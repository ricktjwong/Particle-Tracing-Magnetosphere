#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 14:35:03 2017

@author: ricktjwong
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as spi
from mpl_toolkits.mplot3d import Axes3D

#plt.rcParams["figure.figsize"] = (4,3)
#plt.rcParams["font.family"] = "Times New Roman"
#plt.rcParams['lines.linewidth'] = 0.8

""" Global quantities """

mu = 4*np.pi*1e-7
q = 1.6 * 1e-19
mp = 1.67 * 1e-27
m = np.array([0.0, 0.0, 5e22]) 
r0 = np.array([0.0, 0.0, 0.0])
RE = 6.4e6

def b_field(m, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(m, r_diff))/r_mag**5
    term2 = m/r_mag**3
    B = mu*(term1 - term2)/(4*np.pi)
    uniform = [0,0,5e-9]    # option for constant B field
    return B

def deriv(x,t):
    xx, xy, xz = x[0], x[1], x[2]   # Initial conditions position
    vx, vy, vz = x[3], x[4], x[5]   # Initial conditions velocity
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    B = b_field(m, x, r0)
    if np.linalg.norm(x) > RE:
        a = q * (np.cross(v,B) + E) / mp
        print 1
    else:
        a = np.array([0,0,0])
        v = np.array([0,0,0])
        print 2
    return (v[0], v[1], v[2], a[0], a[1], a[2])
    
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
    
xinit = [-5*RE, 0.0, 0.5*RE, 400e3, 0.0, 0.0]
x0 = np.array([xinit[0], xinit[1], xinit[2]])
v0 = np.array([xinit[3], xinit[4], xinit[5]])
E = np.array([1e-2, 0, 0])
B0 = b_field(m,[xinit[0],xinit[1],xinit[2]],r0)
binit = np.linalg.norm(B0)
r = mp*xinit[3]/(q*binit)   # Larmar radius
T = 2*np.pi*r/xinit[3]      # Gyroperiod particle drift
print r
print T
print B0

t = np.linspace(0,14,2000)

soln = spi.odeint(deriv,xinit,t)    # Solve ODE

x, y, z = soln[:,0], soln[:,1], soln[:,2]
vx, vy, vz = soln[:,3], soln[:,4], soln[:,5]

plt.figure(1)
plt.plot(x,y)
#plt.xlim([-4e7, 0])
#plt.ylim([0,1])
plt.ticklabel_format(useOffset=False)
plt.xlabel("position, x")
plt.ylabel("position, y")

plt.figure(2)
plt.plot(x,z)
plt.xlabel("position, x")
plt.ylabel("position, z")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.plot(x,y,z)
plt.show()

posr = np.vstack((x,y,z)).T

def CheckEnter(r):
    Enter = False
    for i in r:
        dist = np.linalg.norm(i)
        if dist < RE :
            Enter = True
        elif dist > RE :
            Enter = False
        if Enter == True:
            print "Enter"
#        elif Enter == False:
#            print "did not enter"
            
#Track = CheckEnter(posr)

""" Energy Calculations """

v_xyz = np.vstack((vx,vy,vz)).T
v_ = []
for i in v_xyz:
    v_.append(np.linalg.norm(i))
v_ = np.array(v_)
KE = 0.5*mp*v_**2

mod_r = []    
for i in posr:
    mod_r.append(np.linalg.norm(i))
    
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

""" finding tangential vel (along B field lines) """
#
#v_para_ = []
#for i,j in zip(v_xyz, posr):
#    B_hat = (b_field(m,j,r0) / np.linalg.norm(j))
#    v_para_.append(np.dot(i , B_hat)*B_hat)
#v_para_mod =[]
#for i in v_para_:
#    v_para_mod.append(np.linalg.norm(i))
#    
#v_parax=[]
#v_paray=[]
#v_paraz=[]    
#
#for i in v_para_:
#    v_parax.append(i[0])
#    v_paray.append(i[1])
#    v_paraz.append(i[2])
#
#v_perp_ = []
#for i,j in zip(v_xyz, posr):
#    B_hat = (b_field(m,j,r0) / np.linalg.norm(j))
#    v_perp_.append(np.cross(i, B_hat))
#v_perp_mod =[]
#for i in v_perp_:
#    v_perp_mod.append(np.linalg.norm(i))
#    
#plt.figure()
#plt.plot(t , v_para_mod , 'r')
#plt.title("Parallel Velocity Vs Time")
#plt.xlabel('time')
#plt.ylabel('parallel velocity')
#
#plt.figure()
#plt.plot(t , v_perp_mod , 'r')
#plt.title("Perpendicular Velocity Vs Time ")
#plt.xlabel('time')
#plt.ylabel('perp velocity')
#
#plt.figure()
#plt.plot(x , v_para_mod , 'r')
#plt.title("Parallel Velocity Vs Distance")
#plt.xlabel('x-distance/m')
#plt.ylabel('parallel velocity')
#
#plt.figure()
#plt.plot(mod_r, v_para_mod ,' b')
#plt.title("Parallel Velocity Vs AbsDistance")
#plt.xlabel('AbsDistance/m')
#plt.ylabel('parallel velocity')
#
#plt.figure()
#fig4, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax3 = ax1.twinx()
#
#ax1.plot(x,v_parax,'r')
#ax2.plot(x,v_paray,'g')
#ax3.plot(x,v_paraz,'b')
#
#ax1.set_xlabel('x-distance/m')
#ax1.set_ylabel('v_parax', color='red')
#ax1.tick_params('y', colors='red')
#
#ax2.set_ylabel('v_paray', color='g')
#ax2.tick_params('y', colors='g')
#
#ax3.set_ylabel('v_paraz', color='b')
#ax3.tick_params('y', colors='b')
#
#plt.show()
#          
#    