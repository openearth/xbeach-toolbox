import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
from matplotlib import cm
import sys
import os


## import xbeach tools
sys.path.append(os.path.join('..','scripts' ))
from xbeachtools import xgrid, ygrid, seaward_extend, XBeachModelSetup

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

fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(X, Y, bathy, cmap=cm.coolwarm,  linewidth=0, antialiased=False)


###############################################################################
###  x grid                                                                 ###
###############################################################################

xgr,zgr = xgrid(x, bathy[20,:],dxmin=2)


plt.figure()
plt.plot(x,bathy[20,:],'-o')
plt.plot(xgr,zgr,'.-')


###############################################################################
###  y grid                                                                 ###
###############################################################################

ygr = ygrid(y)

plt.figure()
plt.plot(y,'-o')
plt.plot(ygr,'.-')

###############################################################################
###  interpolate                                                            ###
###############################################################################

f = interpolate.interp2d(x, y, bathy, kind='linear')

zgr = f(xgr,ygr)

plt.figure()
plt.pcolor(xgr,ygr,zgr)

xgr, ygr = np.meshgrid(xgr,ygr)
###############################################################################
###  seaward extend                                                         ###
###############################################################################


xgr, ygr, zgr = seaward_extend(xgr,ygr,zgr,slope=1/20,depth=-20)

plt.figure()
plt.pcolor(xgr,ygr,zgr)

plt.figure()
plt.plot(xgr[:,:].T,zgr[:,:].T)

###############################################################################
###  create model setup                                                     ###
###############################################################################

xb_setup = XBeachModelSetup('Test som 1')

print(xb_setup)

xb_setup.set_grid(xgr,ygr,zgr)

xb_setup.set_params({'Wavemodel':'surfbeat',
                     'befriccoef':0.01})

xb_setup.write_model('test')
