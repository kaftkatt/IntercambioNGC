import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

sponge='True'

project='highRes'
location='/HOME/users/amelia.th/Model/IntercambioNGC/'
run='/run8/'
def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=1280  # ocean depth in meters
nx=1    # number of gridpoints in x-direction
ny=240    # number of total gridpoints in y-direction

nz=160     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s

dy=1
dx=1
dz=8

y=np.arange(0,ny*dy,dy)
x=np.arange(0,nx*dx,dx)
z=np.arange(0,nz*dz,dz)

lat=(22.92+31.70)/2

omega=7.2921*(10**(-5))

f=2*omega*sin(np.deg2rad(lat))

def gaussian(x, mu, sig):
    return (
         np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

#Defining heatflux

temp=15*np.ones((nx,ny,nz),dtype=dtype)
temp[:,int(ny*0.5):,:]=temp[:,int(ny*0.5):,:]+1

#Q created to make the plot something
Q=np.zeros((nx,ny,nt),dtype=dtype)

#Values of temperature at open boundary
tempBo=15*np.ones((nx,nz,nt),dtype=dtype)

#Values of salinity at open boundary
salBo=np.zeros((nx,nz,nt),dtype=dtype)

#Bathymetry definition
dep=(gaussian(np.arange(ny), ny/2, 10)*800-Ho)
ho=np.tile(dep,(nx,1))
ho[:,-1]=0


#Saving the files
f = open('temp.bin','wb')
write_with_byte_inversion(f, temp)
f.close()

f = open('bathy.bin','wb')
write_with_byte_inversion(f, ho)
f.close()

f = open('tempOBCS.bin','wb')
write_with_byte_inversion(f, tempBo)
f.close()

f = open('salOBCS.bin','wb')
write_with_byte_inversion(f, salBo)
f.close()

#### Plotting #####

params = {'font.size': 8,
          'figure.figsize': (8, 6),
         'font.family':'sans'}
pl.rcParams.update(params)

fig3 = plt.figure()
gs = fig3.add_gridspec(2, 4)
ax = fig3.add_subplot(gs[0, 0:2])

ax.plot(y,temp[0,:,6],'k' )
ax.set(xlabel='Distance [km]',ylabel='Inital Temperature [C°]')
ax.minorticks_on()

if sponge=='True':
        ax.axvline(y[20],linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[20]+10,-600, 'Sponge\nend', va='center',color='brown',fontsize='small')


ax.grid(which='major',alpha=0.7)

ax.text(0, 1.02, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[1, 0:2])

ax.plot(tempBo[0,:,0],-z )
ax.set(xlabel='Temperature Open boundary [°C]',ylabel='Depth [m]')
ax.text(0, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax.transAxes)


ax = fig3.add_subplot(gs[:, 2])

cax=ax.pcolormesh(z,y,temp[0,:,:])

cbar= plt.colorbar(cax)
cbar.set_label('Inital Temperature [°C]')
ax.set(xlabel='Depth [m]',ylabel='Distance [km]')
if sponge=='True':
	ax.axvline(y[20],linestyle='dashed',c='brown',linewidth=1)
	ax.text(y[20]+2,2, 'Sponge\nend', va='center',color='brown',fontsize='small')

ax.text(-0.3, 1.02, '(c)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[:, 3])

cax=ax.pcolormesh(-ho.T)

if sponge=='True':
	ax.axhline(y[20],linestyle='dashed',c='white',linewidth=1)
	ax.text(0,y[20]+10, 'Sponge\nend', va='center',color='white',fontsize='small')


cbar= plt.colorbar(cax)
cbar.set_label('Depth [m]')
ax.set(ylabel='Distance [km]',xlabel='Distance [km]')
ax.text(0, 1.02, '(d)', fontweight='bold', color='k', 
        transform=ax.transAxes)

plt.tight_layout()
plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png',dpi=400)
plt.show()
