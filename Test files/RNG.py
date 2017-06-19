# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:11:39 2017

@author: Danielsrq
"""

import numpy as np
import random

## where lower_r upper_r are range of coordinates ExpectV is expectation of final V and n is size


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
    
rng1 = rng_posvel(-20, -10, -10, 10, -10, 10, 5, 5)        
