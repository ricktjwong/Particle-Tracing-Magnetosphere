#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 01:49:33 2017

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from VTK_to_Numpy import vtk_subs
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import odeint
from scipy.constants import electron_mass, proton_mass, elementary_charge

me = electron_mass
mp = proton_mass
qe = elementary_charge

RP = 6.4e6

plt.close('all')


#%% Load Data

# Load grid spacing
x,y,z = vtk_subs.import_space('source/x13_rho1-2000.vti')
dx = x[1]-x[0]
dy = y[1]-y[0]
dz = z[1]-z[0]

# Load variables
rho = vtk_subs.import_scalar('source/x13_rho1-2000.vti','rho1')
E = vtk_subs.import_vector('source/x13_Evec-2000.vti','Evec')
B = vtk_subs.import_vector('source/x13_Bvec-2000.vti','Bvec')

fig = plt.figure()
ax = fig.gca(projection='3d')

x, y, z = np.meshgrid(np.linspace(-60, 60, 120),
                      np.linspace(-60, 60, 120),
                      np.linspace(-60, 60, 120))

Bx = B[:,:,:,0]
By = B[:,:,:,1]
Bz = B[:,:,:,2]
#Bx1 = Bx[::10,::10,::10]
#By1 = By[::10,::10,::10]
#Bz1 = Bz[::10,::10,::10]

print np.shape(x)
print np.shape(Bx1)

ax.quiver(x[::10,::10,::10], y[::10,::10,::10], z[::10,::10,::10], Bx[::10,::10,::10], By[::10,::10,::10], Bz[::10,::10,::10], length=12000)
#ax.quiver(x, y, z, Bx, By, Bz, length=0.2)
#ax.quiver(x, y, z, Bx1, By1, Bz1, length=10)


#ax.set_xlim3d(-60,60)
#ax.set_ylim3d(-60,60)
#ax.set_zlim3d(-60,60)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.show()