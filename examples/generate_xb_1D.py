import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
from matplotlib import cm
import sys
import os


## import xbeach tools
sys.path.append(os.path.abspath(os.path.join('..')))

from xbTools import xgrid, ygrid, seaward_extend, XBeachModelSetup, offshore_depth
plt.style.use(os.path.join('..','xbTools','xb.mplstyle'))
###############################################################################
###  load data                                                              ###
###############################################################################

## load data
bathy = np.loadtxt('clean//bathy.dep')
bathy = bathy[0,:]

## set bathy grid
nx = 124
ny = 72
dx = 5
dy = 20

x = np.linspace(0,(nx-1)*dx,nx)

## plot
plt.figure()
plt.plot(x,bathy)
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.title('bathy')




###############################################################################
###  x grid                                                                 ###
###############################################################################

xgr,zgr = xgrid(x, bathy,dxmin=2)


plt.figure()
plt.plot(x,bathy,'-o')
plt.plot(xgr,zgr,'.-')
plt.legend(['Bathy','xgr'])
plt.xlabel('x [m]')
plt.ylabel('z [m]')



###############################################################################
###  interpolate                                                            ###
###############################################################################


zgr = np.interp(xgr, x, bathy)

plt.figure()
plt.plot(x,bathy,'-o')
plt.plot(xgr,zgr,'.-')
plt.xlabel('x [m]')
plt.ylabel('x [m]')
plt.title('xb bathy')


###############################################################################
###  seaward extend                                                         ###
###############################################################################

d_start, slope, Hm0_shoal = offshore_depth(Hm0=9, Tp=15, depth_offshore_profile=abs(bathy[0]), depth_boundary_conditions=20)

xgr, ygr, zgr = seaward_extend(xgr,[0],zgr,slope=slope,depth=d_start*-1)



plt.figure()
plt.plot(xgr.T,zgr[:,:].T)
plt.xlabel('x [m]')
plt.ylabel('z [m]')

###############################################################################
###  create model setup                                                     ###
###############################################################################



xb_setup = XBeachModelSetup('Test som 2')

print(xb_setup)

xb_setup.set_grid(xgr,None,zgr)

xb_setup.set_waves('params',{'Hm0':2,'Tp':5,'gammajsp':3.3, 's' : 10000, 'mainang':270,'fnyq':1})
#xb_setup.set_waves('jonstable',{'Hm0':[1.5, 2, 1.5],'Tp':[4, 5, 4],'gammajsp':[3.3, 3.3, 3.3], 's' : [20,20,20], 'mainang':[270,280, 290],'duration':[3600, 3600, 3600],'dtbc':[1,1,1]})

xb_setup.set_params({'Wavemodel':'surfbeat',
                     'morphology':0,
                     'befriccoef':0.01,
                     'tstop':3600,
                     'nglobalvar':['zb','zs','H'],
                     'npointvar':['zs','zb'],
                     'nmeanvar':['zb'],
                     'npoints':['1 0', '6 0', '10 0', '12 0']})



sim_path = os.path.join('xb-1D')
if not os.path.exists(sim_path):
    os.mkdir(sim_path)
xb_setup.write_model(sim_path)


