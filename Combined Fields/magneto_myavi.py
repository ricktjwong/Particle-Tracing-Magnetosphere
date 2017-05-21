#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 12:50:16 2017

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from VTK_to_Numpy import vtk_subs
import mayavi.mlab as mlab

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

def plotQuiver_B():    
    obj = mlab.quiver3d(x, y, z, Bx, By, Bz,
                   line_width=0.5, 
                   scale_factor=4.0, 
                   scale_mode='scalar', 
                   mode='arrow', 
                   mask_points=300)
    return obj
    
    
def plotQuiver_E():    
    obj = mlab.quiver3d(x, y, z, Ex, Ey, Ez,
                   line_width=0.5, 
                   scale_factor=4.0, 
                   scale_mode='scalar', 
                   mode='arrow', 
                   mask_points=300)
    return obj

#B_field = plotQuiver_B()
E_field = plotQuiver_E()
#mlab.show(B_field)
mlab.show(E_field)
