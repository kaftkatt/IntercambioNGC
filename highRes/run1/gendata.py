import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

sponge='True'

project='highRes'
location='/HOME/users/amelia.th/Model/IntercambioNGC/'
run='/run1/'
def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=1280  # ocean depth in meters
nx=1    # number of gridpoints in x-direction
ny=240    # number of total gridpoints in y-direction
nyA=210    # number of active gridpoints in y-direction
nyS=20    # number of longer dy gridpoints in y-direction
nySend=10
nz=160     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s

dx=1
dyS=5
dy=1
dz=8

x=np.arange(0,nx*dx,dx)
y=np.append(np.append(np.arange(0,nyS*dyS,dyS),np.arange(nyS*dyS,nyA*dy+(nyS*dyS),dy)),np.arange(nyA*dy+(nyS*dyS),nySend*dyS+(nyA*dy+(nyS*dyS)),dyS))
z=np.arange(0,nz*dz,dz)

lat=(22.92+31.70)/2

omega=7.2921*(10**(-5))

f=2*omega*sin(np.deg2rad(lat))

def gaussian(x, mu, sig):
    return (
         np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

#Bathymetry definition
dep=(gaussian(np.arange(nyA), nyA/2, 10)*800-Ho)
ho=np.ones((nx,ny))*(-Ho)
ho[:,nyS:-nySend]=np.tile(dep,(nx,1))
ho[:,-1]=0


#Values of temperature at open boundary
tempBo=15*np.ones((nx,nz,nt),dtype=dtype)

#Values of salinity at open boundary
salBo=np.zeros((nx,nz,nt),dtype=dtype)

#Defining heatflux
out=gaussian(y, 0.75*y[-1], 10)
Q=np.zeros((nx,ny,nt), dtype=dtype)
for i in range(12):
    Q[:,:,i]=-750*out*(gaussian(np.arange(12), 5.5, 2)-np.min(gaussian(np.arange(12), 5.5, 2)))[11-i]

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

f = open('salOBCS.bin','wb')
write_with_byte_inversion(f, salBo)
f.close()

#### Plotting #####
#### Plotting #####

params = {'font.size': 8,
          'figure.figsize': (8, 6),
         'font.family':'sans'}
pl.rcParams.update(params)

fig3 = plt.figure()
gs = fig3.add_gridspec(2, 3)
ax = fig3.add_subplot(gs[0, 0:2])

ax.plot(y,Q[0,:,6],'k' )
ax.plot(y,ho[0,:],'r')
ax.set(xlabel='Distance [km]',ylabel='Heatflux [W/m$^2$]')
ax.minorticks_on()

if sponge=='True':
        ax.axvline(y[20],linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[20]+10,-600, 'Sponge\nend', va='center',color='brown',fontsize='small')
        ax.axvline(y[-nyS],linestyle='dashed',c='black',linewidth=1)
        ax.text(y[-nyS]+10,-600, 'High res\nend', va='center',color='black',fontsize='small')


ax.grid(which='major',alpha=0.7)

ax.text(0, 1.02, '(a)', fontweight='bold', color='k',
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[1, 0:2])

ax.plot(tempBo[0,:,0],-z )
ax.set(xlabel='Temperature [°C]',ylabel='Depth [m]')
ax.text(0, 1.02, '(b)', fontweight='bold', color='k',
        transform=ax.transAxes)


ax = fig3.add_subplot(gs[:, 2])

cax=ax.pcolormesh(y,range(nt),Q[0,:,:].T)
ax.text(y[-1]+100,11,f'Average over active basin:\n{np.sum(np.mean(np.mean(Q[:,nyS:,:],axis=0),axis=1))/len(np.arange(100,y[-1],dy)):.3} W/m$^2$')
cbar= plt.colorbar(cax)
cbar.set_label('Heatflux [W/m$^2$]')
ax.set(xlabel='Distance [km]',ylabel='Time [Months]')
if sponge=='True':
        ax.axvline(y[nyS],linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[nyS]+2,2, 'Sponge\nend', va='center',color='brown',fontsize='small')
        ax.axvline(y[-nyS],linestyle='dashed',c='black',linewidth=1)
        ax.text(y[-nyS]+2,2, 'High res\nend', va='center',color='black',fontsize='small')

ax.text(-0.3, 1.02, '(c)', fontweight='bold', color='k',
        transform=ax.transAxes)


plt.tight_layout()
plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png',dpi=400)
plt.show()
