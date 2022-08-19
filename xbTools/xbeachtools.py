import numpy as np
import os
from datetime import datetime
import json
import matplotlib.pyplot as plt
from shapely.geometry import Point

from xbTools.wave_functions import dispersion, wavecelerity


def xb_run_script_win(xb, N, maindir, xbeach_exe):
    '''
    Create batch script to run simulations.

    Parameters
    ----------
    xb : LIST or CLASS
        list with simulation paths or XBeachModelSetup class
    N : int
        Number of batch scripts.
    maindir : TYPE
        path where run script is created.
    xbeach_exe : TYPE
        path of XBeach executable.

    Returns
    -------
    None.

    '''
    ## check xb, and create path_sims list
    if isinstance(xb,list):
        if isinstance(xb[0],XBeachModelSetup):
            path_sims = []
            for item in xb:
                path_sims.append(item.model_path)
        elif isinstance(xb[0],str):
            path_sims = [xb]
        else:
            print('Unvalid path')
    else:
        path_sims = [xb.model_path]
    
    ## number of batch scripts
    Nmax = int(np.ceil(len(path_sims)/N))
    
    ## string
    string      = ''
    count       = 0
    run_number  = 0
    for ii, path_sim in enumerate(path_sims):
        string = string + 'cd {} \ncall {}\n'.format(path_sim,xbeach_exe)
        if count==Nmax:
            with open(os.path.join(maindir,'run{}.bat'.format(run_number)), 'w') as f:
                f.write(string)
            count = 0
            run_number = run_number + 1
            string = ''
        count = count +1
    if count<=Nmax:
        print(os.path.join(maindir,'run{}.bat'.format(run_number)))
        with open(os.path.join(maindir,'run{}.bat'.format(run_number)), 'w') as f:
            f.write(string)
    
    

def _celerity_ratio_equals_09(Tp,d_start):
    '''
    function to find water depth for which the n ration equal 0.9

    Parameters
    ----------
    Tp : float
        peak period.
    d_start : float
        Water depth.

    Returns
    -------
    d : float
        depth.

    '''
    d_dummy     = d_start
    count2      = 1
    n           = 1
    while n>0.9:
        cg, n       = wavecelerity(Tp, d_dummy)
        d_dummy     = d_dummy + 0.05;
        count2      = count2+1;
        if count2>500:
            print('No n value found!')
            break 
    d = d_dummy
    return d

def offshore_depth(Hm0, Tp, depth_offshore_profile, depth_boundary_conditions):
    '''
    
    compute required ofsshore water depth to correctly force the waves
    Parameters
    ----------
    Hm0 : float
        Wave height.
    Tp : float
        Peak period.
    depth_offshore_profile : float
        Offshore depth of the profile.
    depth_boundary_conditions : float
        Depth of the boundary conditions.

    Returns
    -------
    d_start : float
        Required offshore water depth.
    slope : float
        Artificial slope.
    Hm0_shoal : float
        Wave height at the boundary.

    '''
    cg_bc, dummy            = wavecelerity(Tp, depth_boundary_conditions)
    cg_profile, n_profile   = wavecelerity(Tp, depth_offshore_profile)
    
    Hm0_shoal           = Hm0 * np.sqrt(cg_bc/cg_profile)
    
    if Hm0_shoal/depth_offshore_profile < 0.3 and n_profile<0.9:
        slope   = None
        d_start = depth_offshore_profile
        print('No extension required')
        print('Hm0,shoal = {}'.format(Hm0_shoal))
        print('d start = {}'.format(depth_offshore_profile))
        print('Hm0,shoal/d = {}'.format(Hm0_shoal/depth_offshore_profile))
        print('n = {}'.format(n_profile))
    else:
        ## compute d_start en Hm0,shoal iterative
        d_start             = depth_offshore_profile
        d_start_previous    = 2 * depth_offshore_profile
        count               = 1
        d_n                 = _celerity_ratio_equals_09(Tp,d_start)
        while np.abs(d_start-d_start_previous)>0.05:
            ## update depth
            d_start_previous = d_start
            ## compute required depth
            d_start         = np.max([3.33333*Hm0_shoal, d_n])
            ## compute hm0 shoal
            cg, n_startdepth        = wavecelerity(Tp, d_start)
            Hm0_shoal               = Hm0 * np.sqrt(cg_bc/cg)
            ## update count
            count =+ 1
            if count>50:
                print('no convergence')
                break
        if Hm0_shoal/depth_offshore_profile>0.3 and n_profile>0.9:
            slope = 0.02
            print('Artificial slope of 1:50')
        else:
            slope = 0.1
            print('Artificial slope of 1:10')

        print('Hm0,shoal = {}'.format(Hm0_shoal))
        print('d start = {}'.format(d_start))
        print('Hm0,shoal/d profile = {}'.format(Hm0_shoal/depth_offshore_profile))
        print('Hm0,shoal/d slope = {}'.format(Hm0_shoal/d_start))
        print('n profile = {}'.format(n_profile))
        print('n slope = {}'.format(n_startdepth))
    return d_start, slope, Hm0_shoal






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
    assert(x.ndim<2,'x must be a matrix')
    assert(y.ndim<2,'y must be a matrix')
    assert(z.ndim<2,'z must be a matrix')
    assert(z.shape==x.shape==y.shape,'shape of input matrix is not the same')
    
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
    distance    = (z0max - depth) * 1/slope
    ## dummy array
    x_dummy     = np.arange(x[0,0]-distance,x[0,0],dx_grid)
    
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

    '''

    ## set default values
    dxdry   = dxmin
    zdry    = wl

    
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
        zgr  = np.interp(xgr,x,z)
    ## spatially varying grid
    else:
        ## water depth
        h = np.maximum(wl-z,eps)
        
        if h[-1]>eps:
            k       = dispersion(2*np.pi/Tm,h[-1])
            ## MATLAB
            # k       = dispersion(2*np.pi/Tm,h[0])
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
        dx  = np.asarray(dx)
        if (dx<=dxmin).all():
            print('Computed dxmax (= {} m) is smaller than the user defined dxmin (= {} m). Grid will be generated using constant dx = dxmin. Please change dxmin if this is not desired.'.format(dxmax,localmin) )
       
        ## chop off depth profile, if limit is exceeded
        if xlast>xstart:
            zgr[-1] = zstart
        ##  reverse order back to offshore --> onshore
        xgr = np.flip(xgr)
        zgr = np.flip(zgr)
        
        ## make sure horizontal reference is similar to the input profile
        xgr = xgr-xgr[-1]+xend
    return xgr, zgr


def grid_transition(cell1, cell2, distance):
    
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


    
def ygrid(y,
           dymin = 5,
           dymax = 20,
           area_type='center',
           maxerror = 0.05,
           transition_distance = -0.1,
           area_size = 0.4):
    '''
    

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
    ygr : TYPE
        DESCRIPTION.

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
                
 
def rotate(x, y, theta):
    '''
    rotates the coordinates (x,y) through an angle theta
    if x,y are matrices, these are first flattened.
    output: rotated arrays x,y

    author: Marlies van der Lugt
    revision: v0

    '''
    rotMatrix = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    coords = np.vstack((x.flatten(), y.flatten()))
    return rotMatrix @ coords


def rotate_grid(xgr, ygr, theta):
    '''
    rotate_grid(xgr,ygr,theta)
    rotates a grid xgr,ygr over the angle theta (in radians)

    author: Marlies van der Lugt
    revision: v0
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
    identifies the rotation angle of the x-axis of a 2D XBeach grid and returns grid in local coordinates
    (i.e. x only cross shore, y only alongshore)

    output: rotated grid x,y and grid angle alpha

    author: Marlies van der Lugt
    revision: v0

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
    # rephrase samples in local grid coordinates for simple interpolation and modification
    Input variables
    :param xs: x-world coordinates of samples
    :param ys: y-world coordinates of samples
    :param x0: x-origin of grid in world coordinates
    :param y0: y-origin of grid in world coordinates
    :param theta: angle of x-dir (nautical
    :return:
    xl: x-coordinates of samples in local coordinates
    yl: y-coordintes of samples in local coordinates

    author: Marlies van der Lugt
    revision: v0
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
    
def grid_refine_grid(xgr,ygr,xfactor = 2, yfactor = 1):
    '''
    grid_refine_grid(xgr,ygr,zgr,xfactor = 2, yfactor = 1, ne_layer=None)
    refines the grid with the factor xfactor in xdirection and yfactor in y direction
    works only on rectilinear grids
    returns refined grid where the grid has kept its coordinates
    
    Author: Marlies van der Lugt
    Revision 0 
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
    
 
    
class XBeachModelSetup():
    '''
    XBeach model setup class
    ''' 
    
    
    def __init__(self,fname):
        self.fname      = fname
        ## by default set wbctype and wavemodel to None
        self.wbctype    = None
        self.wavemodel  = None
        
        self.model_path = None
        
    def __repr__(self):
        return self.fname
    
    
    def set_params(self,input_par_dict):
        
        ## set wavemodel. Default is Surfbeat
        if 'Wavemodel' not in input_par_dict:
            print('No wavemodel defined. Wavemodel is set to Surfbeat')
            self.wavemodel = 'surfbeat'
        else:
            self.wavemodel = input_par_dict['Wavemodel']
        
        ## set wbctype
        if 'wbctype' in input_par_dict:
            self.wbctype = input_par_dict['wbctype'] 


        ## load parameters and categories
        f           = open(os.path.join(os.path.dirname(__file__), 'par.json'),'r')
        par_dict    = json.loads(f.read())
        ## create input dict
        self.input_par = {}
        self.input_par['par'] = {}
        ## loop over input parameters 
        value_added = False
        for input_par in input_par_dict:
            ## loop over categories
            for par_category in par_dict:
                ## if input parameter is in category, add parameter
                if input_par in par_dict[par_category]:
                    ## create category if not exist
                    if not par_category in self.input_par:
                        self.input_par[par_category] = {}
                    ## add parameter and value                    
                    self.input_par[par_category][input_par] = input_par_dict[input_par]
                    value_added = True
            if not value_added:
                self.input_par['par'][input_par] = input_par_dict[input_par]
            
        
    
    def set_grid(self,xgr,ygr,zgr, posdwn=1, xori=0, yori=0, thetamin=-90, thetamax = 90, dtheta=10):
        
        ##
        assert(xgr.shape==zgr.shape,'Shape of xgr is not equal to shape of zgr')
        
        ## 1D model
        if ygr is None or ygr.shape[0]==1:
            self.ygr = None
            ## make 2d matrix
            if xgr.ndim==1:
                self.xgr = xgr[np.newaxis, ...] 
                self.zgr = zgr[np.newaxis, ...]  
            else:
                self.xgr = xgr
                self.zgr = zgr            
        ## 2D model
        else:
            self.ygr = ygr
            self.xgr = xgr
            self.zgr = zgr
        
        ##
        self.nx = self.xgr.shape[1] - 1
        self.ny = self.xgr.shape[0] - 1
        ##
        
        ## 1D
        if ygr is None or ygr.shape[0]==1:
            self.fast1D = True
            self.ny = 0
        else:
            self.fast1D = False 
        ##
        self.posdwn = posdwn
        self.xori   = xori
        self.yori   = yori
        self.yori   = yori
        self.thetamin   = thetamin
        self.thetamax   = thetamax
        self.dtheta     = dtheta
        self.vardx  = 1
    
    def set_waves(self,wbctype, input_struct):
        self.wbctype = wbctype
        ##
        if wbctype=='jonstable':
            required_par = ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc']
        elif wbctype=='jons':
            required_par = ['Hm0','Tp','mainang','gammajsp','s','fnyq']
        else:
            assert False, 'Wrong wbctype'
        
        self.waves_boundary  = {}
        for item in required_par:
            assert item in input_struct, '{} missing'.format(item)
            self.waves_boundary[item] =  input_struct[item]
            
    def set_vegetation(self):
        pass

        
    def set_tide(self):
        pass
        
    def load_model_setup(self,path):
        ## todo
        pass    

        
        
    def write_model(self, path, figure=True):
        self.model_path = path
        path_params = os.path.join(path,'params.txt')
        
        assert os.path.exists(path), '{} does not exist'.format(path)
        
        
        current_date    = datetime.today().strftime('%Y-%m-%d %HH:%mm')
        user            =  os.path.basename(os.path.expanduser('~'))
        
        tabnumber = 20
        
        
        ## waves boundary
        if self.wbctype=='jons':
            if 'Wave boundary condition parameters' in self.input_par:
                self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonswap.txt'
            else:
               self.input_par['Wave boundary condition parameters'] = {}
               self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonswap.txt'
            required_par = ['Hm0','Tp','mainang','gammajsp','s','fnyq']
            with open(os.path.join(path,'jonswap.txt'),'w') as f:
                for par in required_par:
                    f.write('{}\t= {}\n'.format(par,self.waves_boundary[par]).expandtabs(tabnumber))
                
        elif self.wbctype=='jonstable':
            if 'Wave boundary condition parameters' in self.input_par:
                self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonstable.txt'
            else:
               self.input_par['Wave boundary condition parameters'] = {}
               self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonstable.txt'                
            required_par = ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc']
            with open(os.path.join(path,'jonstable.txt'),'w') as f:
                for ii in range(len(self.waves_boundary['Hm0'])):
                    for par in required_par:
                        f.write('{} '.format(self.waves_boundary[par][ii]))
                    f.write('\n')
        
        
        ## create params
        with open(path_params,'w') as f:
            ## meta data
            f.write('%% XBeach model: {} \n'.format(self.fname))
            f.write('%% Params created on {} \n'.format(current_date))
            f.write('%% Params created by {} \n'.format(user))
            f.write('\n')

            ## general
            f.write('%% General \n')
            f.write('\n')
            if self.wavemodel!=None:
                f.write('wavemodel\t= {}\n'.format(self.wavemodel).expandtabs(tabnumber))
            if self.wbctype!=None:
                f.write('wbctype\t= {}\n'.format(self.wbctype).expandtabs(tabnumber))
            f.write('\n')
            
            ## grid
            f.write('%% Grid \n')
            f.write('\n')
            f.write('vardx\t= {}\n'.format(self.vardx).expandtabs(tabnumber))
            f.write('nx\t= {}\n'.format(self.nx).expandtabs(tabnumber))
            f.write('ny\t= {}\n'.format(self.ny).expandtabs(tabnumber))
            f.write('xori\t= {}\n'.format(self.xori).expandtabs(tabnumber))
            f.write('yori\t= {}\n'.format(self.yori).expandtabs(tabnumber))     
            f.write('xfile\t= x.grd\n'.expandtabs(tabnumber))
            if not self.fast1D:
                f.write('yfile\t= y.grd\n'.expandtabs(tabnumber))
            f.write('depfile\t= bed.dep\n'.expandtabs(tabnumber))
            f.write('thetamin\t= {}\n'.format(self.thetamin).expandtabs(tabnumber))
            f.write('thetamax\t= {}\n'.format(self.thetamax).expandtabs(tabnumber))
            f.write('dtheta\t= {}\n'.format(self.dtheta).expandtabs(tabnumber))
            f.write('\n')
            
            ## write input vars
            for par_category in self.input_par:
                ## skip category starting with _
                if par_category[0]=='_':
                    continue
                
                ## write meta
                f.write('%% {} \n'.format(par_category))
                f.write('\n')
                for par in self.input_par[par_category]:
                    f.write('{}\t= {}\n'.format(par,self.input_par[par_category][par]).expandtabs(tabnumber))
                f.write('\n')
            ## write output variables
            if '_Output' in self.input_par:
                f.write('%% Output variables \n')
                f.write('\n')
                for par in self.input_par['_Output']:
                    dummy = self.input_par['_Output'][par]
                    f.write('{}\t= {}\n'.format(par,len(dummy)).expandtabs(tabnumber))
                    assert type(dummy)==list, 'expected a list for {}'.format(par)
                    for item in dummy:
                        f.write('{}\n'.format(item))
                    f.write('\n')
    
        ## write grid x
        with open(os.path.join(path,'x.grd'),'w') as f:
            for ii in range(self.ny+1):
                for jj in range(self.nx+1):
                    f.write('{} '.format(self.xgr[ii,jj]))
                f.write('\n')
        if not self.fast1D:
            ## write grid y
            with open(os.path.join(path,'y.grd'),'w') as f:
                for ii in range(self.ny+1):
                    for jj in range(self.nx+1):
                        f.write('{} '.format(self.ygr[ii,jj]))
                    f.write('\n')                    
       ## write dep
        with open(os.path.join(path,'bed.dep'),'w') as f:
            for ii in range(self.ny+1):
                for jj in range(self.nx+1):
                    f.write('{} '.format(self.zgr[ii,jj]))
                f.write('\n')             
                
        ## write figures
        if figure:
            ## plot and write domain
            self._plotdomain(path)
            ## plot and write wave boundary
            if self.wbctype=='jonstable' or self.wbctype=='jons':
                self._plot_boundary(path)

    def _plot_boundary(self,save_path=None):
        '''
        Plot boundary conditions

        Parameters
        ----------
        save_path : TYPE, optional
            Path were figure is saved. The default is None.

        Returns
        -------
        None.

        '''
        if self.wbctype=='jonstable':
            plt.figure()
            plt.subplot(3,1,1)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['Hm0'],'-o')
            plt.ylabel('$H_{m0}$')
            plt.subplot(3,1,2)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['Tp'],'-o')
            plt.ylabel('$T_{p}$')
            plt.subplot(3,1,3)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['mainang'],'-o')
            plt.ylabel('$D$')
            plt.xlabel('Time')
            if save_path!=None:
                plt.savefig(os.path.join(save_path,'jonstable.png'))
        elif self.wbctype=='jons':
            print('wbctype=jons cannot be plotted')
        else:
            print('Not possible to plot wave boundary')
            
    def _plotdomain(self,save_path=None):
        '''
        Plot the domain. 

        Parameters
        ----------
        save_path : string, optional
            Path were figure is saved. The default is None.

        Returns
        -------
        None.

        '''
        plt.figure(figsize=[10,10])
        if self.fast1D==True:
            plt.subplot(2,1,1)
            plt.plot(np.squeeze(self.xgr),np.squeeze(self.zgr)*self.posdwn)
            plt.xlabel('x')
            plt.ylabel('z')
            plt.subplot(2,1,2)
            plt.plot(np.squeeze(self.xgr)[1:],np.diff(np.squeeze(self.xgr)))
            plt.xlabel('x')
            plt.ylabel('dx')
            plt.subplot(2,1,1)
        else:
            plt.pcolor(self.xgr,self.ygr,self.zgr*self.posdwn)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
        plt.grid('on')
        plt.title(self.fname)
        if save_path!=None:
            plt.savefig(os.path.join(save_path,'domain.png'))