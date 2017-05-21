#% Import Modules

import numpy as np
import matplotlib.pyplot as plt
from VTK_to_Numpy import vtk_subs
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

#%%
''' coordinate plane fields ; index = 60 is approx at x=y=z=0 ''' 
## field on xy plane ; z constant
BZx = B[:,:,60,0]
BZy = B[:,:,60,1]
BZz = B[:,:,60,2]
speedxy = np.sqrt(BZx*BZx + BZy*BZy)
strengthxy = speedxy*2000 / speedxy.max()

## field on xz plane ; y constant
BYx = B[:,60,:,0]
BYy = B[:,60,:,1]
BYz = B[:,60,:,2]
speedxz = np.sqrt(BYx*BYx + BYz*BYz)
strengthxz = speedxz*2000 / speedxz.max()

## field on yz plane ; x constant
BXx = B[60,:,:,0]
BXy = B[60,:,:,1]
BXz = B[60,:,:,2]
speedyz = np.sqrt(BXy*BXy + BXz*BXz)
strengthyz = speedyz*5 / speedyz.max()

## create grid (arbitrary spacing)
grid = np.linspace(-59,60,120)
X,Y = np.meshgrid(grid,grid) ## axes names are dummy names

## streamplots
# xy slice at z=0
plt.figure(1)
plt.streamplot(X,Y,BZx,BZy,linewidth = strengthxy)
plt.xlabel("x / relative units")
plt.ylabel("y / relative units")
plt.show()

# xz slice at y=0
plt.figure(2)
plt.streamplot(X,Y,BYx,BYz,linewidth = strengthxz)
plt.xlabel("x / relative units")
plt.ylabel("z / relative units")
plt.show()

# yz slice at x=0 
plt.figure(3)
plt.streamplot(X,Y,BXy,BXz, linewidth = strengthyz)
plt.xlabel("y / relative units")
plt.ylabel("z / relative units")
plt.show()