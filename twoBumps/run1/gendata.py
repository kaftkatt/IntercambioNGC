import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

spongeN='True'
sponge='True'

project='twoBumps'
location='/HOME/users/amelia.th/Model/IntercambioNGC/'
run='/run1/'
def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=1280  # ocean depth in meters
nx=40    # number of gridpoints in x-direction
ny=500    # number of total gridpoints in y-direction
nyS=20

nz=80     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s
dx=0.5
dz=16

OB=15

x=np.arange(0,nx*dx,dx)
z=np.arange(0,nz*dz,dz)

nyIn=80

dyIn=np.zeros(nyIn)
dyIn[0]=500
for i in np.arange(1,nyIn,1):
    dyIn[i] = dyIn[i-1] + dyIn[i-1]*0.035

ny2=250
dyIn2=np.zeros(ny2)

dyIn2[:80]=np.flip(dyIn)
dyIn2[80:208]=500
dyIn2[208:250]=dyIn[:42]

dy=np.append(dyIn2,np.flip(dyIn2))
y=np.cumsum(dy)

lat=(22.92+31.70)/2

omega=7.2921*(10**(-5))

f=2*omega*sin(np.deg2rad(lat))

def gaussian(x, mu, sig):
    return (
         np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )


#Values of salinity at open boundary
salBo=np.zeros((nx,nz,nt),dtype=dtype)

#Bathymetry definition
yB=np.cumsum(dyIn2)
depIN=(gaussian(yB,yB[int(ny2*0.6)], 10000)*880-Ho)
dep=np.append(depIN,np.flip(depIN))
ho=np.tile(dep,(nx,1))
ho[0,:]=0
ho[-1,:]=0

#Values of temperature at open boundary NORTH
indzSill=len(z[z<-np.max(ho[ho!=0])])
tempN=np.zeros(nz)
tempN[:indzSill]=np.arange(15,13,(13-15)/indzSill)
tempN[indzSill-1:]=np.arange(tempN[indzSill-1],11,(10.99-tempN[indzSill-1])/(len(tempN[indzSill-1:])))
tempBoN=np.tile(np.tile(tempN,(nx,1)).T,(nt,1,1)).T

#Values of temperature at open boundary SOUTH
indz600=len(z[z<600])
tempS=np.zeros(nz)
tempS[:indzSill]=np.arange(16,8.5,(8.5-16)/indzSill)
tempS[indzSill-1:indz600]=np.arange(tempS[indzSill-1],7,(7-tempS[indzSill-1])/(len(tempS[indzSill-1:indz600])))
tempS[indz600-1:]=np.arange(tempS[indz600-1],5,(5.02-tempS[indz600-1])/(len(tempS[indz600:])))
tempBoS=np.tile(np.tile(tempS,(nx,1)).T,(nt,1,1)).T

#Values of temperature in basin 

temp=14*np.ones((nx,ny,nz),dtype=dtype)
temp[:,:OB,:]=tempS
temp[:,-OB:,:]=tempN

#Saving the files
f = open('bathy.bin','wb')
write_with_byte_inversion(f, ho)
f.close()

f = open('tempOBCSn.bin','wb')
write_with_byte_inversion(f, tempBoN)
f.close()

f = open('tempOBCSs.bin','wb')
write_with_byte_inversion(f, tempBoS)
f.close()

f = open('salOBCS.bin','wb')
write_with_byte_inversion(f, salBo)
f.close()

f = open('dy.bin','wb')
write_with_byte_inversion(f, dy)
f.close()

f = open('temp.bin','wb')
write_with_byte_inversion(f, temp)
f.close()

#### Plotting #####


params = {'font.size': 8,
          'figure.figsize': (8, 6),
         'font.family':'sans'}
pl.rcParams.update(params)

fig3 = plt.figure()
gs = fig3.add_gridspec(2, 3)
ax = fig3.add_subplot(gs[0, 0:2])

ax.plot(y/1000,ho[1,:],'green',label='Depth [m]')

ax1=ax.twiny()
ax1.plot(tempN,-z,'blue',label='Northern boundary')
ax1.plot(tempS,-z,'red',label='Southern boundary')
ax1.set(xlabel='Temperature [degC]')

ax2=ax.twinx()
ax2.plot(y/1000,dy,'k',label='Resolution in y', linestyle='dashed')
ax.minorticks_on()

ax.set(xlabel='Distance [km]',ylabel='Depth [m]')
ax.minorticks_on()

fig3.legend(bbox_to_anchor=[0.6, 0.65])

if sponge=='True':
        ax.axvline(y[nyS]/1000,linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[nyS]/1000+5,-200, 'Sponge\nend', va='center',color='brown',fontsize='small')
if spongeN=='True':
        ax.axvline(y[-nyS]/1000,linestyle='dashed',c='brown',linewidth=1)
        ax.text(y[-nyS]/1000+5,-200, 'Sponge\nend', va='center',color='brown',fontsize='small')

ax.grid(which='major',alpha=0.7)

ax.text(0, 1.02, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[1, 0])

caxT=ax.pcolormesh(y[-OB:]/1000,-z,temp[5,-OB:,:].T)
ax.set(xlabel='Distance [km]',ylabel='Depth [m]')
ax.text(0, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax.transAxes)

cbar= plt.colorbar(caxT)
cbar.set_label('Temperature Northern Boundary [°C]')

ax = fig3.add_subplot(gs[1, 1])

caxT=ax.pcolormesh(y[:OB]/1000,-z,temp[5,:OB,:].T)
ax.set(xlabel='Distance [km]',ylabel='Depth [m]')
ax.text(0, 1.02, '(c)', fontweight='bold', color='k', 
        transform=ax.transAxes)

cbar= plt.colorbar(caxT)
cbar.set_label('Temperature southern Boundary[°C]')

ax = fig3.add_subplot(gs[:, 2])

cax=ax.pcolormesh(x,y/1000,-ho.T)

cbar= plt.colorbar(cax)
cbar.set_label('Depth [m]')
ax.set(ylabel='Distance [km]',xlabel='Distance [km]')
ax.text(0, 1.02, '(d)', fontweight='bold', color='k', 
        transform=ax.transAxes)

plt.tight_layout()
plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png',dpi=400,bbox_inches='tight')
plt.show()
