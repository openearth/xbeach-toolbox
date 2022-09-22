import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
from matplotlib import cm
import sys
import os


## import xbeach tools
sys.path.append(os.path.abspath(os.path.join('..' )))

from xbTools import xgrid,rotate_grid,xb_run_script_win, ygrid,grid_world2local, seaward_extend, XBeachModelSetup, offshore_depth, lateral_extend
plt.style.use(os.path.join('..','xbTools','xb.mplstyle'))


###############################################################################
###  input                                                                  ###
###############################################################################

zs0 = 5
Hm0 = 9
Tp  = 15

## rotation
theta = 10

###############################################################################
###  load data                                                              ###
###############################################################################



## load data
bathy = np.loadtxt('clean//bathy.dep')


## set bathy grid
nx = 124
ny = 72
dx = 5
dy = 20

x = np.linspace(0,(nx-1)*dx,nx)
y = np.linspace(0,(ny-1)*dy,ny)

X, Y = np.meshgrid(x,y)

## plot
plt.figure()
plt.pcolor(x,y,bathy)
plt.colorbar()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('bathy')

fig     = plt.figure()
ax      = Axes3D(fig)
surf    = ax.plot_surface(X, Y, bathy, cmap=cm.coolwarm,  linewidth=0, antialiased=False)
plt.xlabel('x [m]')
plt.ylabel('y [m]')

###############################################################################
###  rotate                                                                 ###
###############################################################################

[X_world,Y_world] = rotate_grid(X,Y,np.deg2rad(theta))

plt.figure()
plt.pcolor(X_world,Y_world,bathy)
plt.colorbar()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('bathy')

## determine rotation
xl, yl, alpha= grid_world2local(X_world,Y_world)

print('Rotation: ',np.rad2deg(alpha))

###############################################################################
###  x grid                                                                 ###
###############################################################################

xgr,zgr = xgrid(x, bathy[20,:],dxmin=2,Tm=Tp,wl=zs0)


plt.figure()
plt.subplot(2,1,1)
plt.plot(x,bathy[20,:],'-o')
plt.plot(xgr,zgr,'.-')
plt.legend(['Bathy','xgr'])
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.subplot(2,1,2)
plt.plot(xgr[0:-1],np.diff(xgr))
plt.xlabel('x [m]')
plt.ylabel('z [m]')

###############################################################################
###  y grid                                                                 ###
###############################################################################

ygr = ygrid(y)

plt.figure()
plt.plot(y[:-1],np.diff(y),'-o')
plt.plot(ygr[:-1],np.diff(ygr),'.-')
plt.legend(['y','ygr'])
plt.xlabel('y [m]')
plt.ylabel('dy [m]')

###############################################################################
###  interpolate                                                            ###
###############################################################################

f = interpolate.interp2d(x, y, bathy, kind='linear')

zgr = f(xgr,ygr)

plt.figure()
plt.pcolor(xgr,ygr,zgr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('xb bathy')


xgr, ygr = np.meshgrid(xgr,ygr)

###############################################################################
###  seaward extend                                                         ###
###############################################################################
d_start, slope, Hm0_shoal = offshore_depth(Hm0=Hm0, Tp=Tp, depth_offshore_profile=abs(bathy[0,0])+zs0, depth_boundary_conditions=20)

xgr, ygr, zgr = seaward_extend(xgr,ygr,zgr,slope=slope,depth=-d_start)

plt.figure()
plt.pcolor(xgr,ygr,zgr)

plt.figure()
plt.plot(xgr[:,:].T,zgr[:,:].T)
plt.xlabel('x [m]')
plt.ylabel('z [m]')

###############################################################################
###  lateral extend                                                         ###
###############################################################################


xgr,ygr,zgr = lateral_extend(xgr,ygr,zgr,n=5)

plt.figure()
plt.pcolor(xgr,ygr,zgr)

###############################################################################
###  create model setup                                                     ###
###############################################################################





xb_setup = XBeachModelSetup('Test som 1')

print(xb_setup)

xb_setup.set_grid(xgr,ygr,zgr,alfa=theta)

#xb_setup.set_waves('parametric',{'Hm0':2,'Tp':5,'gammajsp':3.3, 's' : 10000, 'mainang':270,'fnyq':1})
xb_setup.set_waves('jonstable',{'Hm0':[1.5, 2, 1.5],'Tp':[4, 5, 4],'gammajsp':[3.3, 3.3, 3.3], 's' : [20,20,20], 'mainang':[270,280, 290],'duration':[3600, 3600, 3600],'dtbc':[1,1,1]})

xb_setup.set_params({'Wavemodel':'surfbeat',
                     'morphology':1,
                     'befriccoef':0.01,
                     'zs0':zs0,
                     'single_dir':1,
                     'tstop':3600,
                     'nglobalvar':['zb','zs','H'],
                     'npointvar':['zs','H'],
                     'nmeanvar':['zb'],
                     'npoints':['1 0', '6 0', '10 0', '12 0']})


sim_path = os.path.join('xb-2D')
if not os.path.exists(sim_path):
    os.mkdir(sim_path)
xb_setup.write_model(os.path.join(sim_path))

xb_run_script_win(xb_setup, N=1, maindir='', xbeach_exe='d:\\XBeach\\XBeachXFINAL\\xbeach.exe')
