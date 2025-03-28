import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl

project='2DsillOB'
location='/HOME/users/amelia.th/Model/IntercambioNGC/'
run='/run4/'
def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

dtype=np.float32

Ho=1200  # ocean depth in meters
nx=1    # number of gridpoints in x-direction
ny=240    # number of gridpoints in y-direction
nz=80     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s

dx=1
dy=1
dz=15

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
out=gauss2d(X,100,[0.75*ny,nx/2])
Q=np.zeros((nx,ny,nt), dtype=dtype)
for i in range(12):
    Q[:,:,i]=out*-600*((gaussian(np.arange(12), 5.5, 2)-np.min(gaussian(np.arange(12), 5.5, 2)))*5)[11-i]


#Values of temperature at open boundary 
tempBo=15*np.ones((nx,nz,nt),dtype=dtype)

#Values of salinity at open boundary
salBo=np.zeros((nx,nz,nt),dtype=dtype)

#Values of velocities at open boundary, random numbers between average(abs) in 
#outputs from the previous model, run2
Vvel=(np.random.random(size=(nx,nz,nt))*10**(-4)).astype(dtype)
Uvel=(np.random.random(size=(nx,nz,nt))*10**(-1)).astype(dtype)
Wvel=(np.random.random(size=(nx,nz,nt))*10**(-5)).astype(dtype)

#Bathymetry definition
dep=(gaussian(np.arange(ny), ny/2, 10)*20000-Ho)
ho=np.tile(dep,(nx,1))
ho[:,-1]=0

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

f = open('VvelOBCS.bin','wb')
write_with_byte_inversion(f, Vvel)
f.close()
f = open('UvelOBCS.bin','wb')
write_with_byte_inversion(f, Uvel)
f.close()
f = open('WvelOBCS.bin','wb')
write_with_byte_inversion(f, Wvel)
f.close()


fig3 = plt.figure(constrained_layout=True)
gs = fig3.add_gridspec(2, 5)
ax = fig3.add_subplot(gs[0, 0:2])


#### Plotting #####

params = {'font.size': 10,
          'figure.figsize': (8, 6),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

ax.plot(range(nt),Q[0,int(ny*0.75),:],'k' )
ax.set(xlabel='Time [Months]',ylabel='Heatflux [W/m$^2$]')
ax.minorticks_on()

ax.grid(which='major',alpha=0.7)

ax.text(0, 1.02, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[1, 0:2])

ax.plot(tempBo[0,:,0],-z )
ax.set(xlabel='Temperature [°C]',ylabel='Depth [m]')
ax.text(0, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax.transAxes)


ax = fig3.add_subplot(gs[:, 2])

cax=ax.pcolormesh(Q[:,:,6].T)

cbar= plt.colorbar(cax)
cbar.set_label('Heatflux [W/m$^2$]')
ax.set(xlabel='Distance [km]',ylabel='Distance [km]')

ax.axhline(y[20],linestyle='dashed',c='brown',linewidth=1)
ax.text(0,y[20]+2, 'Sponge\nend', va='center',color='brown',fontsize='small')

ax.axhline(y[-20],linestyle='dashed',c='brown',linewidth=1)
ax.text(0,y[-20]+2, 'Sponge\nend', va='center',color='brown',fontsize='small')

ax.text(-0.3, 1.02, '(c)', fontweight='bold', color='k', 
        transform=ax.transAxes)

ax = fig3.add_subplot(gs[:, 4])

cax=ax.pcolormesh(-ho.T)

ax.axhline(y[20],linestyle='dashed',c='white',linewidth=1)
ax.text(0,y[20]+2, 'Sponge\nend', va='center',color='white',fontsize='small')

ax.axhline(y[-20],linestyle='dashed',c='white',linewidth=1)
ax.text(0,y[-20]+2, 'Sponge\nend', va='center',color='white',fontsize='small')


cbar= plt.colorbar(cax)
cbar.set_label('Depth [m]')
ax.set(xlabel='Distance [km]')
ax.yaxis.set_tick_params(labelleft=False)
ax.text(0, 1.02, '(d)', fontweight='bold', color='k', 
        transform=ax.transAxes)

#plt.savefig('test.png')
plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png', bbox_inches='tight')

