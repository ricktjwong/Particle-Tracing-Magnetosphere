#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:06:40 2017

@author: ricktjwong
"""

from matplotlib import rcParams
import numpy as np
from matplotlib import pyplot as plt

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

X = np.linspace(-1, 1)
Y = np.linspace(-1, 1)
Z = np.linspace(-1, 1)

# Magnetic dipole orientated in the z direction
#Bx, By, Bz = dipole(m=[0.0, 0.0, 0.1], r=np.meshgrid(X, Y, Z), r0=[0.0,0.0,0.0])
#print np.shape(Bz)
print dipole(m=[0.0, 1.0, 3.0], r=[0.0, 1.0, 2.0], r0=[0.0, 1.0, 0.0])

#Bx, By = dipole(m=[0.0, 0.1], r=np.meshgrid(X, Y), r0=[0.0,0.0])
#
#plt.figure(figsize=(5, 5))
#plt.streamplot(X, Y, Bx, By)
#plt.margins(0, 0)


# Works in 3 dimensions as well:

#X1 = np.linspace(-1, 1, 6)
#print X1
#Y1 = np.linspace(-1, 1, 6)
#Z1 = np.linspace(-1, 1, 6)
#Bx, By, Bz = dipole(m=[0.0, 0.0, 3.0], r=np.meshgrid(X1,Y1,Z1), r0=[0, 0, 0])
#print np.shape(Bx)
#print Bx
