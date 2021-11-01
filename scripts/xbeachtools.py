import numpy as np
import os
from datetime import datetime
import json

def dispersion(w,d,max_error=0.00001):
    '''
    Computes the wave number given a radial frequeny and water depth

    Parameters
    ----------
    w : float
        Radial frequency.
    d : float
        water depth.
    max_error : float, optional
        maximum error. The default is 0.00001.

    Returns
    -------
    k : float
        wave number.

    '''
    g   = 9.81
    ## initial guess
    k1  = 1000
    ## iterate untill convergence
    for i in range(100000):
        ## update k
        k       = k1
        ## next iteration
        k1      = w**2.0/g/np.tanh(k * d)
        ## compute error
        error   = abs(k1-k)
        if error < max_error:
            break   
    ## no convergence
    if i==99999:
        print ('Warning: no convergence')
    return k 



def seaward_extend(x,y,z,slope=1/20,depth=-20):
    
    ## maximum bed level. 
    z0max = np.max(z[:,0])
    
    ## dx
    dx_grid = x[0,1]-x[0,0]
    
    ## maximum distance
    distance    = (z0max - depth) * 1/slope
    x_dummy     = np.arange(x[0,0]-distance,x[0,0],dx_grid)
    
    x_extend    = np.ones((x.shape[0], len(x_dummy) ))
    x_extend    = x_extend * x_dummy
    z_extend    = np.ones_like(x_extend)
    
    y_extend    = np.ones_like(x_extend)
        
    z_extend = z_extend * depth
    
    
    for ii, z0 in enumerate(z[:,0]):
        if z0 < depth:
            continue
        
        
        dz = z0 - depth
        dx = 1/slope * dz
        
        xnew = np.arange(x[ii,0]-dx,x[ii,0],dx_grid)
        
        xp = np.array([x[ii,0]-dx,x[ii,0]])
        zp = np.array([depth,z0])
        
        znew = np.interp(xnew, xp, zp)
        
        N = len(znew)
        
        z_extend[ii,-N:] = znew
        y_extend[ii,:] = y[ii,0]
    
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
          maxfac = 1.15):
    '''
    Computes optimal xgrid 

    Parameters
    ----------
    x : array
        x points of the bathymetry.
    z : array
        bathymetry.
    ppwl : integer, optional
        Number of points per wave length. The default is 20.
    dxmin : TYPE, optional
        DESCRIPTION. The default is 5.
    dxmax : TYPE, optional
        DESCRIPTION. The default is np.inf.
    vardx : TYPE, optional
        DESCRIPTION. The default is 1.
    wl : TYPE, optional
        DESCRIPTION. The default is 0.
    eps : TYPE, optional
        DESCRIPTION. The default is 0.01.
    Tm : TYPE, optional
        DESCRIPTION. The default is 8.
    xdry : TYPE, optional
        DESCRIPTION. The default is None.
    zdry : TYPE, optional
        DESCRIPTION. The default is None.
    dxdry : TYPE, optional
        DESCRIPTION. The default is None.
    depthfac : TYPE, optional
        DESCRIPTION. The default is 2.
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
            Lshort  = 2*np.pi/k
            Lwave   = 4 * Lshort
        else:
            Lwave = 0
        
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
                dx[i] = np.minimum(dx[i], dxmax)
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
                
 
    
 
    
class XBeachModelSetup():
    
    
    def __init__(self,fname):
        self.fname = fname
        
    def __repr__(self):
        return self.fname
    
    
    def set_params(self,input_par_dict):
        
        if 'wavemodel' not in input_par_dict:
            print('no wavemodel defined')
        
        f   = open(os.path.join(os.path.dirname(__file__), 'par.json'),'r')
        par_dict = json.loads(f.read())
        
        
        ## {'grid: {'nx: 10}}
        self.input_par = {}
        
        ## loop over input parameters 
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
        
    
    def set_grid(self,xgr,ygr,zgr):
        ## check dimensions
        
        
        ##
        self.nx = xgr.shape[1]
        self.ny = xgr.shape[0]
        ##
        
        ## 1D
        if ygr is None or ygr.shape[0]==1:
            self.fast1D = True
            self.ny = 0
        else:
            self.fast1D = False 
        ##
        self.posdwn = -1
        self.vardx = 0
    
    def set_waves(self):
        pass
        
    def set_tide(self):
        pass
        
    

        
        
    def write_model(self, path):
        path_params = os.path.join(path,'params.txt')
        
        assert os.path.exists(path), '{} does not exist'.format(path)
        
        current_date = datetime.today().strftime('%Y-%m-%d')
        ## create params
        with open(path_params,'w') as f:
            ## meta data
            f.write('%% XBeach model: {} \n'.format(self.fname))
            f.write('%% Params created on {} \n'.format(current_date))
            f.write('\n')

            ## grid
            f.write('%% General \n')
            f.write('\n')
            
            ## grid
            f.write('%% Grid \n')
            f.write('\n')
            f.write('nx \t = {}\n'.format(self.nx))
            f.write('ny \t = {}\n'.format(self.ny))
            
            f.write('ny \t = {}\n'.format(self.ny))
            
            
            for par_category in self.input_par:
                ## write meta
                f.write('%% {} \n'.format(par_category))
                f.write('\n')
                for par in self.input_par[par_category]:
                    f.write('{} \t = {}\n'.format(par,self.input_par[par_category][par]))
            

        pass

    def _plot_boundary(self):
        pass

    def _plotdomain(self):
        pass