#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:06:40 2017

@author: ricktjwong
"""

from matplotlib import rcParams
import numpy as np
from matplotlib import pyplot as plt
import mayavi.mlab as mlab

rcParams.update({'figure.autolayout': True})

def dipole(m, r, r0):
    """Calculate a field in point r created by a dipole moment m located in r0.
    Spatial components are the outermost axis of r and returned B.
    """
    # we use np.subtract to allow r and r0 to be a python lists, not only np.array
    R = np.subtract(np.transpose(r), r0).T
    
    # assume that the spatial components of r are the outermost axis
    norm_R = np.sqrt(np.einsum("i...,i...", R, R))
    
    # calculate the dot product only for the outermost axis,
    # that is the spatial components
    m_dot_R = np.tensordot(m, R, axes=1)

    # tensordot with axes=0 does a general outer product - we want no sum
    B = 3 * m_dot_R * R / norm_R**5 - np.tensordot(m, 1 / norm_R**3, axes=0)
    
    # include the physical constant
    B *= 1e-7

    return B

def plotQuiver_B():    
    obj = mlab.quiver3d(x, y, z, Bx, By, Bz,
                   line_width=0.1, 
                   scale_factor=0.1, 
                   scale_mode='scalar', 
                   mode='arrow', 
                   mask_points=100)
    return obj

X, Y, Z = np.linspace(-1, 1), np.linspace(-1, 1), np.linspace(-1, 1)
x, y ,z = np.meshgrid(X, Y, Z)

# Magnetic dipole orientated in the z direction
Bx, By, Bz = dipole(m=[0.0, 0.0, 0.1], r=np.meshgrid(X, Y, Z), r0=[0.0,0.0,0.0])

B_field = plotQuiver_B()
mlab.show(B_field)