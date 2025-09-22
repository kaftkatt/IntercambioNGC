import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

project='SpallGOTnoRBCS'

def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=600  # ocean depth in meters
nx=40    # number of gridpoints in x-direction
ny=120    # number of gridpoints in y-direction
nz=10     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s

dx=2.5
dy=5
dz=60

x=np.arange(0,nx*dx,dx)
y=np.arange(0,ny*dy,dy)
z=np.arange(0,nz*dz,dz)

lat=(22.92+31.70)/2

omega=7.2921*(10**(-5)) 

f=2*omega*sin(np.deg2rad(lat))

def gaussian(x, mu, sig):
    return (
        1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

def gauss2d(mat, sigma, center):
    gsize = np.shape(mat)
    [R,C] = np.meshgrid(np.arange(0,gsize[0],1),np.arange(0,gsize[1],1))
    mat = gaussC(R,C, sigma, center)
    return mat

def gaussC(x, y, sigma, center):
    xc = center[0]
    yc = center[1]
    exponent = (np.power((x-xc),2) + np.power((y-yc),2))/(2*sigma)
    val= np.exp(-exponent)
    return val

#Defining heatflux

X,Xno=np.meshgrid(x,y)
out=gauss2d(X,100,[80,10])
Q=np.zeros((nx,ny,nt), dtype=dtype)
for i in range(12):
    Q[:,:,i]=out*-300*((gaussian(np.arange(12), 5.5, 2)-np.min(gaussian(np.arange(12), 5.5, 2)))*5)[11-i]

#Defining Relaxation mask
#Input temperature and relaxation values
tempRelax=15*np.ones((nz,ny,nx),dtype=dtype)

#Values of temperature at open boundary 
tempBo=15*np.ones((nx,nz,nt),dtype=dtype)

#Bathymetry definition
ho=-Ho*np.ones((nx,ny),dtype=dtype)
ho[0,:]=0
ho[-1,:]=0
ho[:,-1]=0

'''
#Plotting it 
params = {'font.size': 8,
          'figure.figsize': (7, 4),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig3 = plt.figure(constrained_layout=True)
gs = fig3.add_gridspec(2, 5)
ax = fig3.add_subplot(gs[0, 0:2])

ax.plot(y,Q[10,:,0] )
ax.set(xlabel='Distance [km]',ylabel='Heatflux [W/m$^2$]')

ax = fig3.add_subplot(gs[1, 0:2])

ax.plot(tempBo[10,:,0],-z )
ax.set(xlabel='Temperature [°C]',ylabel='Depth [m]')


ax = fig3.add_subplot(gs[:, 2])
cax=ax.pcolormesh(x,y,Q[:,:,1].T)
ax.axis('equal')
cbar= plt.colorbar(cax)
cbar.set_label('Heatflux [W/m$^2$]')
ax.set(xlabel='Distance [km]',ylabel='Distance [km]')


ax = fig3.add_subplot(gs[:, 3])

cax=ax.pcolormesh(x,y,ho.T)
ax.axis('equal')
cbar= plt.colorbar(cax)
cbar.set_label('Depth [m]')
ax.set(xlabel='Distance [km]')
ax.yaxis.set_tick_params(labelleft=False)

ax = fig3.add_subplot(gs[:, 4])

cax=ax.pcolormesh(x,y,T_mask[0,:,:].T)
ax.axis('equal')
cbar= plt.colorbar(cax)
cbar.set_label('Restoring temperature area ')
ax.set(xlabel='Distance [km]')
ax.yaxis.set_tick_params(labelleft=False)



fig3.patch.set_alpha(0.0)
plt.savefig('../../IntercambioNGC/Figures/' + str(project) + '/InputVars.png', bbox_inches='tight')
'''

#Saving the files
f = open('Qnet.bin','wb')
write_with_byte_inversion(f, Q)
f.close()

f = open('bathy.bin','wb')
write_with_byte_inversion(f, ho)
f.close()

f = open('tempOBCS.bin','wb')
write_with_byte_inversion(f, tempBo)
f.close()


f = open('temperature.bin','wb')
write_with_byte_inversion(f, tempRelax)
f.close()
