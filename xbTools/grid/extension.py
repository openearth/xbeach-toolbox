# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:10:02 2023

@author: Menno de Ridder, Cas van Bemmelen
collection that allows for the extension of basic XBeach grids
module contains:
    - lateral_extend
    - seaward_extend
"""
import numpy as np

def lateral_extend(x,y,z,n=5):
    '''
    Extend the model domain at both lateral sides with n number of cells

    Parameters
    ----------
    x : array
        x-grid.
    y : TYPE
        y-grid.
    z : array
        bathymetry.
    n : int, optional
        extension at both lateral sides. The default is 5.

    Returns
    -------
    xnew : array
        x-grid.
    ynew : array
        y-grid.
    znew : array
        bathymetry

    '''
    #assert x.ndim<2,'x must be a matrix'
    #assert y.ndim<2,'y must be a matrix'
    #assert z.ndim<2,'z must be a matrix'
    #assert z.shape==x.shape==y.shape,'shape of input matrix is not the same'
    
    dy1 = y[1,0]-y[0,0]
    dy2 = y[-1,0]-y[-2,0]
    
    xnew = np.zeros((x.shape[0]+n*2, x.shape[1]))
    ynew = np.zeros((y.shape[0]+n*2, y.shape[1]))
    znew = np.zeros((z.shape[0]+n*2, z.shape[1]))
    
    xnew[n:-n] = x
    ynew[n:-n] = y
    znew[n:-n] = z
    
    for i in range(n):
        ## update x
        xnew[i,:]       = x[0,:]
        xnew[-(i+1),:]  = x[-1,:]
        ## update z
        znew[i,:]       = z[0,:]
        znew[-(i+1),:]  = z[-1,:]       
        ## update y
        ynew[i,:]       = y[0,:]-dy1*n+dy1*i
        ynew[-(i+1),:]  = y[-1,:]+dy2*n-dy2*i
    return xnew, ynew, znew

def seaward_extend(x,y,z,slope=1/20,depth=-20):
    '''
    Compute the seaward extend of the bathymery based on an artificial  slope and required offshore depth

    Parameters
    ----------
    x : array
        x coordinates of the grid.
    y : array
        y coordinates of the grid.
    z : array
        array with the bathymetry. positive upwards
    slope : float, optional
        artificial slope applied at the offshore to boundary. The default is 1/20.
    depth : float, optional
        Required offshore depth at the boundary. The default is -20.

    Returns
    -------
    xgr : array
        x grid.
    ygr : array
        y grid.
    zgr : array
        bathymetry.

    '''
    if len(z.shape)==1:
        z = np.reshape(z,(1,len(z)))
        x = np.reshape(x,(1,len(x)))
        y = np.reshape(y,(1,len(y)))

    ## maximum bed level at offshore boundary. 
    z0max = np.max(z[:,0])
    
    ## dx at offshore boundary. It assumes a constant dx at the boundary!
    dx_grid = x[0,1]-x[0,0]
    
    ## maximum distance
    distance    = (z0max - depth)/slope
    ## prevent very small grid sizes!
    distance    = np.ceil(distance/dx_grid) * dx_grid
    ## dummy array
    x_dummy     = np.arange(x[0, 0]-distance, x[0, 0], dx_grid)
    
    x_extend    = np.ones((x.shape[0], len(x_dummy) ))
    x_extend    = x_extend * x_dummy
    z_extend    = np.ones_like(x_extend)
    y_extend    = np.ones_like(x_extend)
    z_extend    = z_extend * depth
    
    for ii, z0 in enumerate(z[:,0]):
        if z0 < depth:
            continue
        ## required dx and dz
        dz = z0 - depth
        dx = 1/slope * dz
        
        xnew = np.arange(x[ii,0]-dx,x[ii,0],dx_grid)
        
        xp = np.array([x[ii,0]-dx,x[ii,0]])
        zp = np.array([depth,z0])
        
        znew = np.interp(xnew, xp, zp)
        
        N = len(znew)
        
        z_extend[ii,-N:]    = znew
        y_extend[ii,:]      = y[ii,0]
    
    xgr = np.concatenate((x_extend,x),1)
    zgr = np.concatenate((z_extend,z),1)
    ygr = np.concatenate((y_extend,y),1)
    
    return xgr, ygr, zgr