import numpy as np
import array
import sys
import matplotlib.pyplot as plt
from math import sin
import pylab as pl
from matplotlib.gridspec import GridSpec
import cmocean
from scipy.interpolate import make_smoothing_spline

spongeN='True'
sponge='True'

figshow='TRUE'
saveBathy='TRUE'
figsave='FALSE'
baty='real' #options are: 'real', '2sills', '1sill','flat'
closedBC='yes' #options are: 'no' or 'yes', this is for Y, x are always closed


def write_with_byte_inversion(f, tab):
    tab2write = array.array('f', tab.T.flatten())
    if sys.byteorder == 'little':
        tab2write.byteswap()
    f.write(tab2write)

def read_with_byte_inversion(file_path,rows,cols):
    # 1. Calculate the total expected number of float elements
    total_elements = rows * cols
    
    # 2. Open file in binary read mode ('rb')
    with open(file_path, "rb") as f:
        # Create an empty float array and fill it from the file
        tab_read = array.array('f')
        tab_read.fromfile(f, total_elements)
    
    # 3. Swap the bytes back if the current machine is little-endian
    if sys.byteorder == 'little':
        tab_read.byteswap()
    
    # 4. Convert the native Python array directly to a 1D NumPy array
    # 'f' corresponds to 32-bit floating point numbers (np.float32)
    data_1d = np.asarray(tab_read, dtype=np.float32)
    
    # 5. Reshape using Fortran order ('F') to account for the original .T.flatten()
    original_matrix = data_1d.reshape((rows, cols), order='F')
    
    return original_matrix

dtype=np.float32

Ho=80*19  # ocean depth in meters
nx=40    # number of gridpoints in x-direction
ny=910    # number of total gridpoints in y-direction
nyS=25

nz=80     # number of gridpoints in z-direction
nt=12 #nt=10  #for 3 years with timestep 60 s
dx=0.5
dz=19


x=np.arange(0,nx*dx,dx)
z=np.arange(0,nz*dz,dz)

nyIn=nyS

dyIn=np.zeros(nyIn)
dyIn[0]=500
for i in np.arange(1,nyIn,1):
    dyIn[i] = dyIn[i-1] + dyIn[i-1]*0.1

ny2=455
dyIn2=np.zeros(ny2)

dyIn2[:nyIn]=np.flip(dyIn)
dyIn2[nyIn:]=500

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
salBo=np.ones((nx,nz,nt),dtype=dtype)*34

#Bathymetry definition

if baty == 'real':
    ho = read_with_byte_inversion("bathyReal.bin",nx, ny)
elif baty == '2sills':
    depIN=(gaussian(y[:int(ny*0.5)],y[int(ny*0.3)], 8000)*880-Ho)
    dep=np.append(depIN,np.flip(depIN))
    ho=np.tile(dep,(nx,1))
elif baty == '1sill':
    dep=(gaussian(y,y[int(ny*0.5)], 8000)*880-Ho)
    ho=np.tile(dep,(nx,1))
elif baty == 'flat': 
    ho=np.ones((nx,ny))*-Ho

if closedBC == 'yes':
    ho[:,0]=0
    ho[:,-1]=0  

ho[0,:]=0
ho[-1,:]=0



#Values of temperature profile
indzSill=len(z[z<400])
indz600=len(z[z<600])
tempS=np.zeros(nz)
tempS[:indzSill]=np.arange(16,8,(8-16)/indzSill)
tempS[indzSill-1:indz600]=np.arange(tempS[indzSill-1],7,(7-tempS[indzSill-1])/(len(tempS[indzSill-1:indz600])))
tempS[indz600-1:]=np.arange(tempS[indz600-1],5,(5.02-tempS[indz600-1])/(len(tempS[indz600:])))
splnew = make_smoothing_spline(z, tempS, lam=1e6)
tempNew=splnew(z)

#Temperature at Boundary
tempBoS=np.tile(np.tile(tempNew,(nx,1)).T,(nt,1,1)).T

#Values of temperature in basin 

temp=tempNew*np.ones((nx,ny,nz),dtype=dtype)

#Values of passive tracer
trac=np.zeros((nx,ny,nz),dtype=dtype)
trac[:,:int(ny2*0.55),indzSill-5:]=1

# Values of heatflux
Q=np.zeros((nx,ny,nt))
Qin=-80*(gaussian(np.arange(12), 5.5, 2))
Q=np.tile(Qin,(nx,ny,1))

#Relaxation mask for rbcs 
scurve = np.zeros_like(y)
t=(y[:nyS]/y[nyS])
curve=1 - (3*t**4-2*t**6)
scurve[:5]=1
scurve[-5:]=1
scurve[5:nyS+5] = curve
scurve[-nyS-5:-5] = np.flip(curve)

relaxmask = scurve[None, :, None]*np.ones([nx, ny, nz])

# Create row indices (y-axis)
depmask = np.linspace(0,nz*dz,nz)[:, None]*-1   # shape (h, 1)

# Broadcast comparison
mask = depmask <= ho[20,:]  # True = below (masked), False = above

#Saving the files
if saveBathy=='TRUE':
	f = open('bathy.bin','wb')
	write_with_byte_inversion(f, ho)
	f.close()

f = open('relax_mask.bin','wb')
write_with_byte_inversion(f, relaxmask)
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

f = open('trac.bin','wb')
write_with_byte_inversion(f, trac)
f.close()

f = open('Qnet.bin','wb')
write_with_byte_inversion(f, Q)
f.close()

#### Plotting #####
params = {'font.size': 15,
          'figure.figsize': (10, 5),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

levels=np.arange(6,16,1)
fig = plt.figure()
gs = GridSpec(nrows=4, ncols=10,hspace=0.2,wspace=0.001)

ax = fig.add_subplot(gs[0, :8])
ax.plot(y/1000,Q[20,:,6],color='red',linewidth=2,label=f'Mean over active domain {np.mean(Q):.1f}')
ax.set_xlim((np.min(y/1000),np.max(y/1000)))
ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False) 
ax.set_ylabel(r'Q [W/m$^{2}$]',color='red')
ax.tick_params(axis='y', colors='red')
ax.legend()

ax2=ax.twinx()
ax2.plot(y/1000,dy/1000,color='k',linestyle='dotted',linewidth=2)
ax2.plot(y/1000,relaxmask[20,:,10],color='purple',linewidth=1)
ax2.set(ylabel=r'$\Delta$y [km]')
ax.minorticks_on()
 

ax.grid(which='minor',linestyle='--', alpha=0.5)
ax.grid(which='major',alpha=0.7)
ax2.spines['left'].set_color('red')

ax = fig.add_subplot(gs[1:, :])
cax=ax.pcolormesh(y/1000,-z,np.ma.masked_array(temp[20].T,mask=mask),cmap=cmocean.cm.thermal)
caxC=ax.contour(y/1000,-z,np.ma.masked_array(temp[20].T,mask=mask),levels=levels,colors='white')
ax.clabel(caxC, inline=1, fontsize=8)
cbar=plt.colorbar(cax)
cbar.set_label('Temperature [degC]')

#ax.plot(y/1000,deptry,'red')
ax.set(xlabel='Distance [km]',ylabel='Depth [m]')

ax.axvline(y[nyS]/1000,c='black',linestyle='dashed')
ax.text(0.06, 0.15, 'End of \nSponge', color='black', 
    transform=ax.transAxes,fontsize=8)
ax.axvline(y[-nyS]/1000,c='black',linestyle='dashed')
ax.text(0.8, 0.15, 'End of \nSponge', color='black', 
    transform=ax.transAxes,fontsize=8)

if figshow=='TRUE':
	plt.show()
if figsave=='TRUE':
	plt.savefig(str(location)+str(project)+'/Figures'+str(run)+'InputVars.png',dpi=400,bbox_inches='tight')
