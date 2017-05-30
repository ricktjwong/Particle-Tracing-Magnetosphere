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

x, y, z = np.meshgrid(np.linspace(-60, 60, 120),
                      np.linspace(-60, 60, 120),
                      np.linspace(-60, 60, 120))

Bx, By, Bz = B[:,:,:,0], B[:,:,:,1], B[:,:,:,2]
Ex, Ey, Ez = E[:,:,:,0], E[:,:,:,1], E[:,:,:,2]

def plotQuiver_B(scale):    
    ax.quiver(x[::scale,::scale,::scale], 
              y[::scale,::scale,::scale],
              z[::scale,::scale,::scale], 
              Bx[::scale,::scale,::scale],
              By[::scale,::scale,::scale],
              Bz[::scale,::scale,::scale], length=5, arrow_length_ratio=.5, normalize=True)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
def plotQuiver_E(scale):    
    ax.quiver(x[::scale,::scale,::scale], 
              y[::scale,::scale,::scale],
              z[::scale,::scale,::scale], 
              Ex[::scale,::scale,::scale],
              Ey[::scale,::scale,::scale],
              Ez[::scale,::scale,::scale], length=10, arrow_length_ratio=.5, normalize=True)
        
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

fig = plt.figure(1)
ax = fig.gca(projection='3d')

plotQuiver_B(15)

#ax.set_xlim3d(-60,60)
#ax.set_ylim3d(-60,60)
#ax.set_zlim3d(-60,60)

fig2 = plt.figure(2)
ax = fig2.gca(projection='3d')
plotQuiver_E(20)