import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

sponge='True'

project='highRes'
location='/HOME/users/amelia.th/Model/IntercambioNGC/'
run='/run10/'
def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=1280  # ocean depth in meters
nx=1    # number of gridpoints in x-direction
ny=240    # number of total gridpoints in y-direction
nyS=20

nz=160     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s

dx=1
dz=8

inc=0.12
ninc=30
start=1 + inc
maxi=ninc*inc+start
dy=1000*np.append(np.ones(nyS)*maxi,np.append(np.arange(maxi,start,-inc),np.append(np.ones(ny-nyS-ninc*2),np.arange(start,maxi,inc))))

y=np.cumsum(dy)/1000
x=np.arange(0,nx*dx,dx)
z=np.arange(0,nz*dz,dz)

umbral=ny*0.4
lat=(22.92+31.70)/2

omega=7.2921*(10**(-5))

f=2*omega*sin(np.deg2rad(lat))

def gaussian(x, mu, sig):
    return (
         np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

#Defining heatflux
out=gaussian(y, 0.8*y[-1], 10)
Q=np.zeros((nx,ny,nt), dtype=dtype)
for i in range(12):
    Q[:,:,i]=-750*out*(gaussian(np.arange(12), 5.5, 2)-np.min(gaussian(np.arange(12), 5.5, 2)))[11-i]

temp=15*np.ones((nx,ny,nz),dtype=dtype)
temp[:,:int(umbral),110:]=temp[:,:int(umbral),110:]-1

#Q created to make the plot something
#Q=np.zeros((nx,ny,nt),dtype=dtype)

#Values of temperature at open boundary
tempBo=14*np.ones((nx,nz,nt),dtype=dtype)

#Values of salinity at open boundary
salBo=np.zeros((nx,nz,nt),dtype=dtype)

#Bathymetry definition
dep=(gaussian(np.arange(ny), umbral, 10)*800-Ho)
ho=np.tile(dep,(nx,1))
ho[:,-1]=0


#Saving the files
f = open('Qnet.bin','wb')
write_with_byte_inversion(f, Q)
f.close()

f = open('bathy.bin','wb')
write_with_byte_inversion(f, ho)
f.close()

f = open('temp.bin','wb')
write_with_byte_inversion(f, temp)
f.close()

f = open('tempOBCS.bin','wb')
write_with_byte_inversion(f, tempBo)
f.close()

f = open('salOBCS.bin','wb')
write_with_byte_inversion(f, salBo)
f.close()

f = open('dy.bin','wb')
write_with_byte_inversion(f, dy)
f.close()
#### Plotting #####

params = {'font.size': 8,
          'figure.figsize': (8, 6),
         'font.family':'sans'}
pl.rcParams.update(params)

fig3 = plt.figure()
gs = fig3.add_gridspec(2, 4)
ax = fig3.add_subplot(gs[0, 0:2])

ax.plot(y,Q[0,:,6],'k',label='Heatflux [W/m$^2$]')
ax2=ax.twinx()
ax2.plot(y,ho[0,:],'red',label='Depth [m]')
ax3=ax.twinx()
ax3.plot(y,dy,'k',linestyle='dashed',label='Resolution [m]')
ax.set(xlabel='Distance [km]',ylabel='Heatflux [W/m$^2$]')
ax.minorticks_on()
ax.text(0, 0.9, f'Mean heatflux {np.mean(Q[:,nyS:,:]):.1f} W/m$^2$', color='k',
        transform=ax.transAxes)

fig3.legend(bbox_to_anchor=[0.6, 1.05],ncol=2)

if sponge=='True':
        ax.axvline(y[20],linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[20]+10,-600, 'Sponge\nend', va='center',color='brown',fontsize='small')


ax.grid(which='major',alpha=0.7)

ax.text(0, 1.02, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[1, 0:2])

ax.pcolormesh(y,-z,temp[0,:,:].T)
ax.plot(y,ho[0,:])
ax.set(xlabel='Temperature south of the sill [°C]',ylabel='Depth [m]')
ax.text(0, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax.transAxes)


ax = fig3.add_subplot(gs[:, 2])

cax=ax.pcolormesh(np.arange(nt),y,Q[0,:,:])

cbar= plt.colorbar(cax)
cbar.set_label('Heatflux [W/m$^2$]')
ax.set(xlabel='Time [months]',ylabel='Distance [km]')
#if sponge=='True':
#	ax.axvline(y[20],linestyle='dashed',c='brown',linewidth=1)
#	ax.text(y[20]+2,2, 'Sponge\nend', va='center',color='brown',fontsize='small')

ax.text(-0.3, 1.02, '(c)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[:, 3])

cax=ax.pcolormesh(-ho.T)

#if sponge=='True':
#	ax.axhline(y[20],linestyle='dashed',c='white',linewidth=1)
#	ax.text(0,y[20]+10, 'Sponge\nend', va='center',color='white',fontsize='small')


cbar= plt.colorbar(cax)
cbar.set_label('Depth [m]')
ax.set(ylabel='Distance [km]',xlabel='Distance [km]')
ax.text(0, 1.02, '(d)', fontweight='bold', color='k', 
        transform=ax.transAxes)

plt.tight_layout()
plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png',dpi=400,bbox_inches='tight')
plt.show()
