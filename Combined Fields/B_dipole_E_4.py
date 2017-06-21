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
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['lines.linewidth'] = 0.8

""" Global quantities """

mu = 4*np.pi*1e-7
q = 1.6 * 1e-19
mp = 1.67 * 1e-27
m = np.array([0.0, 0.0, -5e22]) 
r0 = np.array([0.0, 0.0, 0.0])
RE = 6.4e6
t = np.linspace(0,180,4500)
omega = np.array([0.0, 0.0, 7.2921150e-5])

def b_field(m, r, r0):
    r_diff = r - r0
    r_mag = np.linalg.norm(r_diff)
    term1 = 3.*r_diff*(np.dot(m, r_diff))/r_mag**5
    term2 = m/r_mag**3
    B = mu*(term1 - term2)/(4*np.pi)
    return B

def e_field(r):
    v_convec = np.array([400e3, 0, 0])
    B_convec = np.array([0, 0, -5e-9])
    E_convec = -np.cross(v_convec, B_convec)
    B_corot = b_field(m, r, r0)
    E_corot = -np.cross(np.cross(omega, r),  B_corot)
    total_E = E_corot + E_convec
    return total_E

def deriv(x,t):
    xx, xy, xz = x[0], x[1], x[2]   # Initial conditions position
    vx, vy, vz = x[3], x[4], x[5]   # Initial conditions velocity
    x=np.array([xx, xy, xz])
    v=np.array([vx, vy, vz])
    B = b_field(m, x, r0)
    E = e_field(x)
    if np.linalg.norm(x) > RE:
        a = q * (np.cross(v,B) + E) / mp
    elif np.linalg.norm(x) < RE:
        a = np.array([0,0,0])
        v = np.array([0,0,0])
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
    
def CheckEnter(trajectory):
    Enter = False
    last_posvel = trajectory[-1]
    xx, xy, xz = last_posvel[0], last_posvel[1], last_posvel[2]
    r_ = np.array([xx,xy,xz])
    r_mod = np.linalg.norm(r_)
    if r_mod <= RE:
        Enter = True
    elif r_mod > RE:
        Enter = False
    return Enter
    
def CheckEnter2(trajectory):
    Enter = False
    step = 50
    sampled_traj = []
    true_list = []
    posr = []
    for i in trajectory[0::step]:
        sampled_traj.append(i)
    for i in sampled_traj:
        xx , xy , xz = i[0] , i[1] , i[2]
        r_ = np.array([xx,xy,xz])
        posr.append(r_)
    for i in posr:
        r_mod = np.linalg.norm(i)
        print r_mod
        if r_mod <= RE:
            true_list.append(True)
            print true_list
    if len(true_list) > 1:
        Enter = True
    else:
        Enter = False
    return Enter
    
def search(pos_vel, step, count, sweep_param):
    soln_set = []
    nonsoln_set = []
    delta = step  # increment of initial x to sweep
    n = 1           # counter
    while n <= count:
        traj = spi.odeint(deriv,pos_vel,t)
        if CheckEnter(traj) == True:
            soln_set.append(pos_vel)
            print traj[-1]
            print "Enter!"
        elif CheckEnter(traj) == False:
            mod_r=[]
            x,y,z = traj[:,0], traj[:,1], traj[:,2]
            posr = np.vstack((x,y,z)).T
            for i in posr:
                mod_r.append(np.linalg.norm(i))
            print np.min(mod_r)
            #print traj[-1]
            print "Not Enter!"
            nonsoln_set.append(pos_vel)
        
        if sweep_param == "posx":
            pos_vel = pos_vel + np.array([delta,0.0,0.0,0.0,0.0,0.0])
        elif sweep_param == "posy":
            pos_vel = pos_vel + np.array([0.0,delta,0.0,0.0,0.0,0.0])
        elif sweep_param == "posz":
            pos_vel = pos_vel + np.array([0.0,0.0,delta,0.0,0.0,0.0])
        elif sweep_param == "velx":
            pos_vel = pos_vel + np.array([0.0,0.0,0.0,delta,0.0,0.0])
        elif sweep_param == "vely":
            pos_vel = pos_vel + np.array([0.0,0.0,0.0,0.0,delta,0.0])
        elif sweep_param == "velz":
            pos_vel = pos_vel + np.array([0.0,0.0,0.0,0.0,0.0,delta])  
        n += 1
    
    return (soln_set, nonsoln_set)

def search2(pos_vel):
    soln_set = []
    nonsoln_set = []
    for i in pos_vel:
        traj = spi.odeint(deriv,i,t)
        if CheckEnter(traj) == True:
            soln_set.append(i)
            print traj[-1]
            print "Enter!"
        elif CheckEnter(traj) == False:
            mod_r=[]
            x,y,z = traj[:,0], traj[:,1], traj[:,2]
            posr = np.vstack((x,y,z)).T
            for j in posr:
                mod_r.append(np.linalg.norm(j))
            print np.min(mod_r)
            #print traj[-1]
            print "Not Enter!"
            nonsoln_set.append(i)
                
    return (soln_set, nonsoln_set)
    

def rng_posvel(lower_x , upper_x , lower_y , upper_y , lower_z , upper_z, ExpectV , n):
    v_sq = ExpectV**2
    Pos_Vel = []
    v_ = []
    
    x_pos = np.random.randint(lower_x, upper_x,size = (n,1))
    y_pos = np.random.randint(lower_y,upper_y,size=(n,1))
    z_pos = np.random.randint(lower_z,upper_z,size=(n,1))
    xyz_pos = np.hstack((x_pos,y_pos,z_pos))
    print xyz_pos
    
    xyz_vel = (v_sq)*np.random.dirichlet(np.ones(3),size = n)
    for i in xyz_vel:
        vx,vy,vz = np.sqrt(i[0]), (random.choice([+1, -1]))*np.sqrt(i[1]), (random.choice([+1, -1]))*np.sqrt(i[2])
        v_.append( [vx, vy, vz] )
    print v_
    
    Pos_Vel = np.hstack((xyz_pos, v_))
    
    return Pos_Vel

#soln_pos, nonsoln_pos = search(xinit, RE, 10, 'pos')
#soln_vel, nonsoln_vel = search(xinit, 100e3, 10, 'velx')

def saveToFile(valuesToWrite, sweep_param):
    text_file = open('Data/'+sweep_param+".txt", "w")
    headers=['Entered','Not entered']
    for i,j in zip(headers,valuesToWrite):
        text_file.write(str(i)+": %s \n\n" % str(j))
    text_file.close()

#random_posvel_array = rng_posvel(-10*RE, -2*RE, -3*RE , 3*RE , -3*RE , 3*RE, 400.0e3 , 10)
#soln_posvel , nonsoln_posvel = search2(random_posvel_array)
#print soln_posvel
#print nonsoln_posvel    
    
#saveToFile([soln_pos, nonsoln_pos], 'posx')
#saveToFile([soln_vel, nonsoln_vel], 'velx')

xinit = [-5*RE, 0.0, 0.5*RE, 400e3, 400e3, 100e3]

x0 = np.array([xinit[0], xinit[1], xinit[2]])
v0 = np.array([xinit[3], xinit[4], xinit[5]])
B0 = b_field(m,[xinit[0],xinit[1],xinit[2]],r0)
binit = np.linalg.norm(B0)
r = mp*xinit[3]/(q*binit)   # Larmar radius
T = 2*np.pi*r/xinit[3]      # Gyroperiod particle drift
print r
print T
print B0

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
    PE.append(q*np.dot((x0-i), e_field(i)))

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
