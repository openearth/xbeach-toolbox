# -*- coding: utf-8 -*-
"""
--- Synopsis --- 
This script converts degrees to u,v vectors. Intensity of a parameter can be included to scale the vectors for quiver plots.

--- Remarks --- 
See also: 
To-Do: 
Dependencies: 

--- Version --- 
Created on Fri October 05 16:34:50 2018
@author: BEMC
Project: general
Script name: deg2uv.py 

--- Revision --- 
To-Do Subversioncontrol 
Status: Unverified 

Witteveen+Bos Consulting Engineers 
Van Twickelostraat 2 
P.O. box 233 
7400 AE Deventer 
The Netherlands 
	
	Give an array of directions will return the U and V
	components. intensity is optional.
	
"""
# 1. Import modules 
import numpy as np 

# 2. Define new functions 
def deg2uv(direction,intensity=np.nan):
    """
    Function converts degrees to vector, denoted as
    u and v
    
    Parameters
    ----------
    direction : numpy array
        direction in degrees in numpy array.
    intensity : numpy array
        length of the vector. The default is np.nan. and 1 will be used

    Returns
    -------
    u : np array
        magnitude in x direction
    v : np array
        magnitude in y direction

    """
    Dir_nan_index = np.isnan(direction)
    
    # Make nan 0 to perform sin cos operations
    Dir_0 = direction
    if np.any(Dir_nan_index):
        Dir_0[Dir_nan_index] = 0
        intensity[Dir_nan_index] = 0
    
    rad = -4.0*np.arctan(1.0)/180.

    if ~np.all(np.isnan(intensity)):
        u   = intensity*np.sin(Dir_0*rad)
        v   = -intensity*np.cos(Dir_0*rad)
    else:
        u   = np.sin(Dir_0*rad)
        v   = -np.cos(Dir_0*rad)
    
    # Replace all original nan with nan
    #u[Dir_nan_index] = np.nan
    #v[Dir_nan_index] = np.nan
    return u,v

def uv2deg(u,v,convention):
    """
    Function converts u,v vector to direction, denoted as
    direction
    
    Parameters
    ----------
    u : numpy array
        x-component (e.g. velocity).
    v : numpy array
        y-component (e.g. velocity)
    convention: nautical or cartesian

    Returns
    -------
    direction : np array
    """

    if convention == 'cartesian':
        direction = np.mod(np.arctan2(u,v)*180/np.pi,360);
    elif convention == 'nautical':
        direction = np.mod(np.arctan2(-u,-v)*180/np.pi,360)
    else:
        print('please provide a valid convention: nautical or cartesian')
            
    return direction