# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:10:02 2023

@author: Menno de Ridder, Cas van Bemmelen
collection that allows for the creation of basic XBeach grids
module contains:
    - xgrid
    - ygrid
    - grid_transition
"""
import numpy as np
from ..general.wave_functions import dispersion

def xgrid(x,z,
          ppwl=20,
          dxmin=5,
          dxmax=np.inf,
          vardx=1,
          wl = 0,
          eps = 0.01,
          Tm = 8,
          xdry=None,
          zdry=None,
          dxdry = None,
          depthfac = 2,
          maxfac = 1.15,
          nonh = False):
    '''
    Compute spatially varying grid based on the local wave length 

    Parameters
    ----------
    x : array
        x points of the bathymetry.
    z : array
        bathymetry (positive upwards).
    ppwl : integer, optional
        Number of points per wave length. The default is 20.
    dxmin : float, optional
        minimum grid resolution. The default is 5.
    dxmax : float, optional
        maximum grid reslution. The default is np.inf.
    vardx : int, optional
        1=spatially varying grid; 0=equidistant  grid resolution. The default is 1.
    wl : float, optional
        Water level. The default is 0.
    eps : float, optional
        Minimum water depth. The default is 0.01.
    Tm : float, optional
        Mean wave period. The default is 8.
    xdry : float, optional
        bathymetry is considered dry for x-grids larger than xdry. The default is None.
    zdry : float, optional
        bathymetry is considered dry for z values larger than zdry. The default is None.
    dxdry : float, optional
        Resolution of dry cells. The default is None.
    depthfac : float, optional
         . The default is 2.
    maxfac : TYPE, optional
        DESCRIPTION. The default is 1.15.

    Returns
    -------
    xgr : array
        grid points.
    zgr : array
        depth points.

    #todo improve docstring and syntax of code
    '''

    ## set default values
    if dxdry is None:
        dxdry   = dxmin
    if zdry is None:
        zdry    = wl

    # make sure x is monitonically increasing:
    if x[0]>x[-1]:
        x = np.flipud(x)
        z = np.flipud(z)


    ## remove nan
    x = x[~np.isnan(z)]
    z = z[~np.isnan(z)]
    
    ## set values
    xstart = x[0]
    zstart = z[0]
    xend   = x[-1]
    zend   = z[-1]
    
    ## equidistant cross-shore grid
    if vardx == 0:
        
        nx = int((xend-xstart)/dxmin)
        
        ## constant dx
        xgr  = np.linspace(xstart,xend,nx+1)
        zgr  = np.interp(xgr, x, z)
    ## spatially varying grid
    else:
        ## water depth
        h = np.maximum(wl-z,eps)
        
        if h[0]>eps:
            # k       = dispersion(2*np.pi/Tm,h[-1])
            ## MATLAB
            k       = dispersion(2*np.pi/Tm,h[0])
            if nonh:
                Lshort  = 2*np.pi/k
                Lwave   = Lshort
            else:
                Lshort  = 2*np.pi/k
                Lwave   = 4 * Lshort               
        else:
            Lwave = 0
        
            ## MATLAB
            #k       = dispersion(2*np.pi/Tm,h[0])
            #Lshort  = 2*np.pi/k
            #Lwave   = 4 * Lshort
        
        ##
        i   = 0
        xgr = [xend]
        zgr = [zend]
        hgr = [h[-1]]
        dx  = []
        xlast   = xend
        while xlast > xstart:
            ## dry cells
            if xdry != None:
                drycell = xgr[i] > xdry
            else:
                drycell = zgr[i] > zdry ## error in Matlab?
            if drycell:
                localmin = dxdry
            else:
                localmin = dxmin  
            ##
            
            ##
            dxmax_cell = Lwave/ppwl
            dxmax_cell = np.minimum(dxmax_cell,dxmax)
            dx.append( np.maximum(depthfac*hgr[i], localmin) )
            
            if dxmax_cell > localmin:
                dx[i] = np.minimum(dx[i], dxmax_cell)
            else:
                dx[i] = localmin
            
            
            ## max factor
            if i>0:
                if dx[i] >=maxfac*dx[i-1]:
                    dx[i] = maxfac*dx[i-1]
                if dx[i] <=1/maxfac*dx[i-1]:
                    dx[i] = 1/maxfac*dx[i-1]
            ## update lists
            i = i + 1
            xgr.append(xgr[i-1] - dx[i-1])
            
            xtemp   = np.minimum(xgr[i],xend)
            hgr.append( np.interp(xtemp, x,h) )
            zgr.append( np.interp(xtemp, x,z) )
            xlast = xgr[i]
            
            if hgr[i]>eps:
                k       = dispersion(2*np.pi/Tm,hgr[i])
                if nonh:
                    Lshort  = 2*np.pi/k
                    Lwave   = Lshort
                else:
                    Lshort  = 2*np.pi/k
                    Lwave   = 4 * Lshort  
            else:
                Lwave = 0
        ##
        xgr = np.asarray(xgr)
        zgr = np.asarray(zgr)
        if (np.min(dx)<dxmin):
            print('Computed dxmax (= {} m) is smaller than the user defined dxmin (= {} m). Grid will be generated using constant dx = dxmin. Please change dxmin if this is not desired.'.format(dxmax,localmin) )
            dxmin = np.min(dx)

        ## chop off depth profile, if limit is exceeded
        if xlast>xstart:
            zgr[-1] = zstart
        ##  reverse order back to offshore --> onshore
        xgr = np.flip(xgr)
        zgr = np.flip(zgr)
        
        ## make sure horizontal reference is similar to the input profile
        xgr = xgr-xgr[-1]+xend
    return xgr, zgr


def ygrid(y,
           dymin = 5,
           dymax = 20,
           area_type='center',
           maxerror = 0.05,
           transition_distance = -0.1,
           area_size = 0.4):
    '''
    function to setup a basic ygrid array

    Parameters
    ----------
    y : TYPE
        DESCRIPTION.
    dymin : TYPE, optional
        DESCRIPTION. The default is 5.
    dymax : TYPE, optional
        DESCRIPTION. The default is 20.
    area_type : TYPE, optional
        DESCRIPTION. The default is 'center'.
    maxerror : TYPE, optional
        DESCRIPTION. The default is 0.05.
    transition_distance : TYPE, optional
        DESCRIPTION. The default is -0.1.
    area_size : TYPE, optional
        DESCRIPTION. The default is 0.4.

    Returns
    -------
    ygr : array
        array with y grid.
    '''
    
    retry = False
    if transition_distance < 0:
         transition_distance = np.abs(transition_distance)
         retry = True
         print('Enable optimization of transition distance')
     
     
    if len(y)==1:
         print('1D model')
         ygr = np.linspace(0,1,1) * dymin
    else:
         if dymin==dymax:
             ygr = np.arange(np.min(y),np.max(y)+dymax,dymax)
             print('Create equidistant alongshore grid')
             print('Grid size {}'.format(dymin))
         else:
             err = np.inf
             
             print('Area type {}'.format(area_type))
             
             while err > maxerror:
                 dy = np.max(y)-np.min(y)
                
                 if transition_distance < 1:
                     transition_distance = transition_distance * dy
                    
                 print('Transition {}'.format(transition_distance))
                 if area_type == 'center':
                     if area_size < 1:
                         area_size = area_size * dy
                     ygr = np.arange(np.mean(y)-area_size/2, np.mean(y)+area_size/2+dymin, dymin)
                 else:
                     ygr = np.mean(y) + np.array([-1, 1]) * dymin/2
                 ## grid transition
                 ff, nf, gridf, err = grid_transition(dymin, dymax, transition_distance)
                 
                 tmp = np.concatenate((ygr[0] - np.flip(gridf), ygr))
                 ygr = np.concatenate((tmp, ygr[-1] + gridf))
                    
                 if retry:
                     if err > maxerror and retry:
                         transition_distance = 1.1 * transition_distance
                     else:
                         break
             if err > maxerror:
                 print('Relative error in alongshore grid transition')
            
             tmp = np.concatenate((np.flip(np.arange( ygr[0] - dymax, np.min(y) - dymax, -1*dymax)), ygr ))
             ygr = np.concatenate((tmp, np.arange( ygr[-1] + dymax, np.max(y) + dymax, dymax) ))
            
    return ygr 

def grid_transition(cell1, cell2, distance):
    
    """
    function for grid_transition from one cell to the other

    Returns:
        ff, nf, gridf, error 
    
    #todo improve docstring and syntax of code
    """  

    precision = 1e-10
    maxloop = 1e2
    maxfac = 1.15
    
    nf = 0
    ff = 1
    gridf = []
    error = 0
    
    #assert(distance>0,'error, zero or negative distance')
    
    cells = np.array([cell1, cell2])
    
    nmin = np.max([1, np.floor(distance/np.max(cells))] )
    nmax = np.max([1, np.ceil(distance/np.min(cells))] )
    
    n = np.arange(nmin,nmax+1,1)
    
    i = 0
    f = []
    
    for ni in n:
        ni = int(ni)
        # start with equidistant grid
        fj      = np.array([1., 1.])
        Lj      = np.ones(3)*cell1*ni
        stage   = 1
        
        j = 0
        while np.abs(Lj[2]- distance ) > precision:
            if stage ==1:
                if cell1 > cell2:
                    fj[0] = 0.9 * fj[0]
                else:
                    fj[1] = 1.1 * fj[1]
            if stage ==2:
                if Lj[2] > distance:
                    fj[1] = np.mean(fj)
                elif Lj[2] < distance:
                    fj[0] = np.mean(fj)
            
            ##
            Lj[0] = cell1 * np.sum( np.power(fj[0], np.arange(1,ni+1,1))  )
            Lj[1] = cell1 * np.sum(np.power(fj[1],np.arange(1,ni+1,1)) )
            Lj[2] = cell1 * np.sum(np.power(np.mean(fj),np.arange(1,ni+1,1)) )
            
            if (Lj[0] > distance and Lj[1] < distance) or ( Lj[0] < distance and Lj[1] > distance ):
                stage = 2
                
            if fj[0] < precision or np.isinf(fj).any() or np.isnan(fj).any() or (np.diff(fj) < precision and j > maxloop):
                fj = np.array([ np.inf, np.inf])
                break
            j = j+1
        
        ##
        f.append( np.mean(fj) )
        
        if f[i] > maxfac:
            f[i] = np.nan
        
        i =i + 1
    ##
    errors =  np.abs(cell1 * np.power(np.asarray(f),(n+1)) - cell2 )
    
    i = np.where(errors==np.nanmin(errors))[0]
    
    if len(i)>0:
        nf = int(n[i])
        ff = np.asarray(f)[i]
        error = errors[i]/cell2
        gridf = np.cumsum(cell1 * np.power(ff,np.arange(1,nf+1,1)) )
    
    return ff, nf, gridf, error