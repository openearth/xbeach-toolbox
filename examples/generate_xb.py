import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
from matplotlib import cm

from scripts.xbeachtools import xgrid

## load data
bathy = np.loadtxt('clean//bathy.dep')

nx = 124
ny = 72
dx = 5
dy = 20

x = np.linspace(0,(nx-1)*dx,nx)
y = np.linspace(0,(ny-1)*dy,ny)

X, Y = np.meshgrid(x,y,)

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