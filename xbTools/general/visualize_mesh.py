# -*- coding: utf-8 -*-
"""
Created on Wed May 5 10:04:00 2023

@author: Cas van Bemmelen
collection that allows for easy mesh plotting and conversion to .shp file for GIS programs
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import geopandas
from shapely.geometry import Point, LineString, Polygon, mapping
import fiona
from fiona.crs import from_epsg

def plot_mesh(mesh_x,mesh_y,thinning=1,ax=None,**kwargs):
    '''
    Function to plot mesh based on meshgrid regular grid, option for thinning

    Parameters
    ----------
    mesh_x : array
        Array with x grid.
    mesh_y : array
        Array with y grid.
    thinning : integer to specify thinning of grid to plot. Default = 1
    ax : specify ax to which the mesh should be added. Default = None, which makes a new figure
    **kwargs : kwargs to feed into the plot function (e.g. color, linewidth)

    Returns
    -------
    none.

    '''
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_axes([0.15, 0.15, 0.85, 0.75])

    # thin on y axis and then on x axis if value > 1
    ax.plot(mesh_x[:,::thinning], mesh_y[:,::thinning], **kwargs)
    ax.plot(mesh_x[::thinning,:].T, mesh_y[::thinning,:].T, **kwargs)

def write_mesh_to_shp(mesh_x,mesh_y,file_name,modelnr,EPSG=3844):
    '''
    Function to write mesh to shapefile

    Parameters
    ----------
    mesh_x : array
        x-coordinates of the mesh
    mesh_y : array
        y-coordinates of the mesh
    file_name : string
        where to store the file
    modelnr : float
        which number to give to the lines
    EPSG : int
        optional CRS choice. The default is 32635.

    Returns
    -------
    None.

    '''
    # create dir if doesnt exists
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))
    
    # create empty array to store lines
    lines = []
    row_nums = []
    col_nums = []
    
    # loop over both directions in i and j indices and store line segment to empty lines array
    for i in range(np.shape(mesh_x)[1]):
        coords = [(mesh_x[0][i],mesh_y[0][i]),(mesh_x[-1][i],mesh_y[-1][i])]
        lines.append(LineString(coords))
        col_nums.append(f"{i:05d}")
        row_nums.append('')

    for j in range(np.shape(mesh_x)[0]):
        coords = [(mesh_x[j,0],mesh_y[j,0]),(mesh_x[j,-1],mesh_y[j,-1])]
        lines.append(LineString(coords))
        col_nums.append('')
        row_nums.append(f"{j:05d}")

    # schema: it is a simple dictionary with geometry and properties as keys
    schema = {'geometry': 'LineString','properties': {'modelnr': 'str',
                                                      'rownr' : 'str',
                                                      'colnr' : 'str'}}
    # for defining the geometry, you need Shapely
    # two simples geometries
    with fiona.open(file_name, 'w',crs=from_epsg(EPSG), driver='ESRI Shapefile', schema = schema) as layer:
        for line, row_num, col_num in zip(lines, row_nums, col_nums):
            # filling schema
            elem = {}
            # geometry with mapping function of shapely
            elem['geometry'] = mapping(line) 
            # attribute value (the same here)
            elem['properties'] = {'modelnr': modelnr,
                                  'rownr' : row_num,
                                  'colnr' : col_num}
            # writing element in the file
            layer.write(elem)

    print(f'{file_name} is saved')
    
def grid_parameter_to_shape(xgrid,ygrid,parameter_name,parameter_values,filename, EPSG = 3844):
    """
    saves a parameter defined at grid points to a shapefile

    Parameters
    ----------
    xgrid : array
        array with x-coordinates of grid
    ygrid : array
        array with y-coordinates of grid
    parameter_name : string
        name of parameter to save to shapefile
    parameter_values : array of size xgrid, ygrid
        parameter to save to shapefile
    filename : path
        where to save shapefile
    EPSG : integer
        EPSG-code of grid. The default is 4326.

    Returns
    -------
    None.

    """    
    
    # check if folder exists
    folder = os.path.dirname(filename)    
    if not os.path.exists(folder):
        os.makedirs(folder)  
    
    # flatten grid and parameter
    xgrid_flatten = np.ndarray.flatten(xgrid)
    ygrid_flatten = np.ndarray.flatten(ygrid)
    parameter_values_flatten = np.ndarray.flatten(parameter_values)
    
    # define points and add data
    points = [Point(x_val, y_val) for x_val, y_val in zip(xgrid.flatten(), ygrid.flatten())]
    shape_data = {'geometry':points,parameter_name:parameter_values_flatten}
    
    # make geodataframe and save
    gdf = geopandas.GeoDataFrame(shape_data, crs = "EPSG:"+str(EPSG))
    gdf.to_file(filename) 