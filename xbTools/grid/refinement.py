# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:10:02 2023

@author: Menno de Ridder, Cas van Bemmelen
collection that contains refinement functions for regular meshes 
module contains:
    - grid_refine_grid
"""
# general import
import numpy as np

# toolbox specific import
from ..general.geometry import rotate_grid

def grid_refine_grid(xgr,ygr,xfactor = 2, yfactor = 1):
    '''
    grid_refine_grid(xgr,ygr,zgr,xfactor = 2, yfactor = 1, ne_layer=None)
    refines the grid with the factor xfactor in xdirection and yfactor in y direction
    works only on rectilinear grids
    returns refined grid where the grid has kept its coordinates
    
    Author: Marlies van der Lugt
    Revision 0 
    # todo: add better docstring and syntax
    '''    
    #rotation of the grid
    alpha = np.arctan2(ygr[0,-1]-ygr[0,0],xgr[0,-1]-xgr[0,0])
    
    #the origin of the grid
    x0 = xgr[0,0]
    y0 = ygr[0,0]
    
    #to local coordinates
    xl,yl = rotate_grid(xgr-x0,ygr-y0,alpha) 

    #refine the xgrid with the factor
    xx = xl[0,:]    
    dx = np.diff(xx)
 
    xr = (xx[:-1,None] + np.linspace(0,dx,xfactor,endpoint=False).T).ravel() # xx is input array
    xr = np.append(xr,xx[-1])
    
    #refine the ygrid with the factor    
    yy = yl[:,0]
    dy = np.diff(yy)
    yr = (yy[:-1,None] + np.linspace(0,dy,yfactor,endpoint=False).T).ravel() # yy is input array
    yr = np.append(yr,yy[-1])
    
    xr,yr = np.meshgrid(xr,yr)
    
    #back to world coordinates
    xnew,ynew = rotate_grid(xr,yr,-alpha)
    
    #add origin again
    xgr2 = xnew + x0
    ygr2 = ynew + y0
    
    return xgr2, ygr2 