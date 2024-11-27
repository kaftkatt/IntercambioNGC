# put this file in the input directory for the barotrophic experiment,
# if you are working in python. This creates the bathymetry and the wind
# forcings and can also be used to implement an initial ssh

import numpy as np
import array
import sys

def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=600  # ocean depth in meters
nx=20    # number of gridpoints in x-direction
ny=100    # number of gridpoints in y-direction
nz=20     # number of gridpoints in z-direction
nt=10  #for 3 years with timestep 60 s

# Flat bottom at z=-Ho
ho=-Ho*np.ones((nx,ny),dtype=dtype)
#dmax=-10 # wind stress maximum
#x=(np.linspace(0,nx-1,nx)-0.5)/(nx-2) # x-coordinate, centered at grid u-point
#y=(np.linspace(0,ny-1,ny)-0.5)/(ny-2) # y-coordinate, centered at grid v-point
#X,Y=np.meshgrid(x,y,indexing='ij')
#ho[0:100,:]=-45*np.sin(2*np.pi*X[10:110]) # generate file for -np.cos(y) profile between 0-1200km
#ho[120:,:]=-45*np.sin(np.pi*X[0:120])-32.5
#ho[80:120,:]=-17*np.sin(4*np.pi*X[70:110])-40
ho[0,:]=0
ho[-1,:]=0
ho[:,-1]=0
f = open('bathy2.bin','wb')
write_with_byte_inversion(f, ho)
f.close()

# Initial salinity 
#s=np.ones((nx,ny,nz), dtype=dtype)
#s[:,:,:]=20
#s[:120,:,30:]=35
#s[120:,:,45:]=35
#f = open('salinit2.bin','wb')
#write_with_byte_inversion(f, s)
#f.close()

Q=np.zeros((nx,ny,nt), dtype=dtype)
Q[:,70:90,1]=-80

f = open('Qnet.bin','wb')
write_with_byte_inversion(f, Q)
f.close()

#ev=np.zeros((nx,ny,nt), dtype=dtype)
#ev[0:20,:,1]=10**(-5)

#f = open('evap2.bin','wb')
#write_with_byte_inversion(f, ev)
#f.close()

# Wind-stress
#tauMax=0.1 # wind stress maximum
#x=(np.linspace(0,nx-1,nx)-0.5)/(nx-2) # x-coordinate, centered at grid u-point
#y=(np.linspace(0,ny-1,ny)-0.5)/(ny-2) # y-coordinate, centered at grid v-point
#X,Y=np.meshgrid(x,y,indexing='ij')
#tau=-tauMax*np.sin(4*np.pi*Y) # generate file for -np.cos(y) profile between 0-1200km
#tau[:,0:60]=tauMax*np.sin(4*np.pi*Y[:,0:60])
#tau=tau*0
#tau[0:12,55:65]=-0.2

#f = open('windx.bin','wb')
#write_with_byte_inversion(f, np.array(tau, dtype=dtype))
#f.close()
#tau=tauMax*np.sin(np.pi*Y/120) # generate file for sin(y) profile between 0-1200km
#f=open('windx_siny.bin','wb')
#write_with_byte_inversion(f, np.array(tau, dtype=dtype))
#f.close()
