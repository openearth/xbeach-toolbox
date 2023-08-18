# -*- coding: utf-8 -*-
"""
Created on Wed May 5 10:04:00 2023

@author: Cas van Bemmelen
collection containing geometry functionalities
"""
# import python modules
import os
import numpy as np
import math
from matplotlib.path import Path
from shapely.geometry import Polygon, Point

def rotate(x, y, theta, origin = [0,0]):
    '''
    Rotates the coordinates (x,y) counterclockwise through an angle theta around origin

    Parameters
    ----------
    x : float
        x-coordinate.
    y : float
        y-coordinate.
    theta : floatd
        angle in radians.
    origin: list length 2
        list of floats [x0, y0] origin coordinate. Default [0, 0]

    Returns
    -------
    rotated arrays x,y

    '''
    x0, y0 = origin
    ny, nx = x.shape
    coords = np.vstack((x.flatten()-x0, y.flatten()-y0))

    rotMatrix = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    xrot, yrot = rotMatrix @ coords

    return xrot.reshape(ny, nx)+x0, yrot.reshape(ny, nx)+y0


def rotate_grid(xgr, ygr, theta):
    '''
    Rotate_grid(xgr,ygr,theta)
    rotates a grid xgr,ygr over the angle theta (in radians)

    Parameters
    ----------
    xgr : array
        x grid.
    ygr : array
        y grid.
    theta : floatd
        ange in radians.

    Returns
    -------
    uv : array
        rotated grid.
    vd : array
        rotated grid.

    '''
    
    ny, nx = xgr.shape
    coords = np.vstack((xgr.reshape(-1), ygr.reshape(-1)))
    rotMatrix = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    uv, vd = rotMatrix @ coords
    uv = uv.reshape([ny, nx])
    vd = vd.reshape([ny, nx])
    return uv, vd
    
def grid_world2local(xgr, ygr):
    '''
    Identifies the rotation angle of the x-axis of a 2D XBeach grid and returns grid in local coordinates
    (i.e. x only cross shore, y only alongshore)
    output: rotated grid x,y and grid angle alpha

    Parameters
    ----------
    xgr : array
        x-grid.
    ygr : array
        y-grid.

    Returns
    -------
    xl : array
        x-grid in local coordinates
    yl : array
        y-grid in local coordinates
    alpha : TYPE
        rotation.
    '''
    
    #rotation of the grid
    alpha = np.arctan2(ygr[0,-1]-ygr[0,0],xgr[0,-1]-xgr[0,0])
    
    #the origin of the grid
    x0 = xgr[0,0]
    y0 = ygr[0,0]
    
    #to local coordinates
    xl,yl = rotate_grid(xgr-x0,ygr-y0,alpha) 
    
    return xl, yl, alpha
    
def samples_to_local_grid(xs,ys,x0,y0,theta):
    '''
    # shift samples in local grid coordinates for simple interpolation and modification
    Parameters
    ----------
    xs: x-world coordinates of samples
    ys: y-world coordinates of samples
    x0: x-origin of grid in world coordinates
    y0: y-origin of grid in world coordinates
    theta: angle of x-dir (nautical
    Returns
    -------
    xl: x-coordinates of samples in local coordinates
    yl: y-coordintes of samples in local coordinates
    '''

    xs2 = xs-x0
    ys2 = ys-y0
    return rotate(xs2,ys2,-theta)

def in_polygon(x,y,poli):
    '''
    checks whether the coordinate (x,y) or list of coordinates (xi,yi) fall 
    within the polygon poli.
    
    Parameters
    ----------
    x : numpy array, either 1D or 2D
        x coordinates.
    y : numpy array, either 1D or 2D
        y coordintes.
    poli : shapely Polygon geometry
        polygon to test.

    Returns
    -------
    ip : numpy array, either 1D or 2D
        mask being 1 if in polyon, 0 outside polygon.

    Author: Marlies van der Lugt
    Revision 0 
    '''
    ny,nx = x.shape
    p = [Point(ix, iy) for ix, iy in zip(x.flatten(),y.flatten())]
    # pdb.set_trace()
    ip = np.array([poli.contains(p[i]) for i in range(len(p))]).reshape(ny,nx)
    return ip

# def rotate(origin, point, angle):
#     """
#     Rotate a point counterclockwise by a given angle around a given origin.
#
#     The angle should be given in radians.
#     """
#     ox, oy = origin
#     px, py = point
#
#     qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
#     qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
#     return qx, qy
   
def get_points_in_polygon(pol, x_coords, y_coords):
    """
    to-do, replace in_polygon by this function 
    returns True for points inside the polygon and False for points outside
    
    Note: If a point falls on the border of the polygon, a False is returned
    Parameters
    ----------
    pol : shapely.geometry.Polygon
        polygon to get points in.
    x_coords : np.array
        x coords
    y_coords : np.array
        y coords

    Returns
    -------
    ind_inside : np.array
        boolean array with True if point inside and False if outside.

    """
    # try to convert if not already in np.array, except raise error
    try:
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
    except:
        raise UserWarning('Please supply x and y coords as np.array, you supplied: '
                          f'{type(x_coords) and {type(y_coords)}}')
    
    # initiate with all not inside the polygon
    ind_inside = np.full_like(x_coords, False, dtype=bool)
    
    # first filter on rough extend
    x_pol, y_pol = pol.exterior.coords.xy
    
    ind_roughly_in = (x_coords <= np.max(x_pol)) & (x_coords >= np.min(x_pol)) & \
                     (y_coords <= np.max(y_pol)) & (y_coords >= np.min(y_pol))
    
    # exit early if nothing is inside the rough extend
    if not np.any(ind_roughly_in):
        return ind_inside
    
    # now check exactly if in     
    # con catenate the coordinates to points to determine if inside
    point_coords = list(zip(x_coords[ind_roughly_in], y_coords[ind_roughly_in]))

    # make a matplotlib path of the polygon to check if points fall in it
    pol_path = Path(list(pol.exterior.coords))
    
    print(f'start determining if {len(point_coords):.2e} points are inside polygon')    
    # get the indices of points that are inside and reshape this to original format
    ind_exactly_inside = pol_path.contains_points(point_coords)
    
    # now we need to combine two indice lists and set some values on True if inside
    ind_inside[ind_roughly_in] = ind_exactly_inside
    
    return ind_inside


def path_distance(polx, poly):
    '''
    computes the distance along a polyline
    ----------
    polx : TYPE array
        X COORDINATES.
    poly : TYPE array
        Y COORDINATES.

    Returns: TYPE array
        PATHDISTANCE
    -------
    python alternative to matlab function.

    '''
    dx = np.diff(polx)
    dy = np.diff(poly)

    dr = np.sqrt(dx ** 2 + dy ** 2)

    pathDistance = np.insert(np.cumsum(dr), 0, 0, axis=0)

    return pathDistance
    