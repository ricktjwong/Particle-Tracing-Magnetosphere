#% Import Modules

import numpy as np
import matplotlib.pyplot as plt
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

r_IB = 6*dx # Inner boundary

# Plot

iy = int(0.5*y.shape[0])

fig, ax = plt.subplots(figsize=(11,5))
p = ax.pcolormesh(x/RP,z/RP,rho[:,iy,:].T)
plt.colorbar(p,ax=ax)
ax.set(aspect='equal',xlabel='x')
plt.show()

#%% Particle Streamtracing

d = np.array([dx,dy,dz])
    
def lorentz(eta,t,m,q,x,y,z,E,B,d):
    
    pi = eta[:3]
    vi = eta[3:]
    
    f_err = [0,0,0,0,0,0]
    
    dh = -0.5*d
    
    ri = np.sqrt(np.sum(pi**2))
    if(pi[0]<=x[0] or pi[1]<=y[0] or pi[2]<=z[0] or 
       pi[0]>x[-1] or pi[1]>y[-1] or pi[2]>z[-1] or 
       ri<6*6.4e6):
        return f_err
    
    Bx = interp(pi+[0,1,1]*dh,x,y,z,B[:,:,:,0],d)
    By = interp(pi+[1,0,1]*dh,x,y,z,B[:,:,:,1],d)
    Bz = interp(pi+[1,1,0]*dh,x,y,z,B[:,:,:,2],d)
    
    Ex = interp(pi+[1,0,0]*dh,x,y,z,E[:,:,:,0],d)
    Ey = interp(pi+[0,1,0]*dh,x,y,z,E[:,:,:,1],d)
    Ez = interp(pi+[0,0,1]*dh,x,y,z,E[:,:,:,2],d)
    
    f = [vi[0], vi[1], vi[2],
         q/m*(Ex + vi[1]*Bz-vi[2]*By),
         q/m*(Ey + vi[2]*Bx-vi[0]*Bz),
         q/m*(Ez + vi[0]*By-vi[1]*Bx)]
    return f

def interp(pi,x,y,z,f,d):
    
    ix = int((pi[0]-x[0])/d[0])
    iy = int((pi[1]-y[0])/d[1])
    iz = int((pi[2]-z[0])/d[2])
    
    xi =  (x[ix+1]-pi[0])/d[0]
    yi =  (y[iy+1]-pi[1])/d[1]
    zi =  (z[iz+1]-pi[2])/d[2]
    
    fi = trilinear(xi,yi,zi,f[ix:ix+2,iy:iy+2,iz:iz+2])
    
    return fi    

def trilinear(x,y,z,f):
    xm = 1-x
    ym = 1-y
    zm = 1-z
    
    if(x<=0 or x>1):
        print(x)
    
    c00 = f[0,0,0]*xm + f[1,0,0]*x
    c01 = f[0,0,1]*xm + f[1,0,1]*x
    c10 = f[0,1,0]*xm + f[1,1,0]*x
    c11 = f[0,1,1]*xm + f[1,1,1]*x
    
    c0 = c00*ym+c10*y
    c1 = c01*ym+c11*y
    
    c = c0*zm+c1*z
    
    return c

#%% Define initial parameters and run

m = mp
q = qe

eta0 = np.array([-18*RP,-5*RP,2*RP,4e3,0,0])

t = np.linspace(0,1,1000)*600
eta = odeint(lorentz,eta0,t,args=(m,q,x,y,z,E,B,d))

# Filter out 'bad' values
el_x = np.logical_and(eta[:,0]>x.min(),eta[:,0]<x.max())
el_y = np.logical_and(eta[:,1]>y.min(),eta[:,1]<y.max())
el_z = np.logical_and(eta[:,2]>z.min(),eta[:,2]<z.max())
el_r = np.sqrt(np.sum(eta[:,:3]**2,axis=1))>r_IB
el = np.logical_and(el_x,el_y)
el = np.logical_and(el,el_z)
el = np.logical_and(el,el_r)

t = t[el]
eta = eta[el]

# Plot
fig, ax = plt.subplots(1,2,sharex=True,figsize=(11,5))
r = np.sqrt(np.sum(eta[:,:3]**2,axis=1))
v = np.sqrt(np.sum(eta[:,3:]**2,axis=1))
ax[0].plot(t,eta[:,:3]/RP,t,r/RP)
ax[1].plot(t,eta[:,3:]/1e3,t,v/1e3)
for axi in ax:
    axi.legend(['x','y','z','mag'])
plt.show()

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(eta[:,0]/RP,eta[:,1]/RP,eta[:,2]/RP)
ax.set(xlim=[x.min()/RP,x.max()/RP],ylim=[y.min()/RP,y.max()/RP],zlim=[z.min()/RP,z.max()/RP],
      xlabel='x',ylabel='y',zlabel='z',
      aspect='equal')
plt.show()



