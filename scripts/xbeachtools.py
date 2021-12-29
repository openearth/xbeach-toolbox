import numpy as np
import os
from datetime import datetime
import json
import matplotlib.pyplot as plt

def dispersion(w, d, max_error=0.00001):
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
                
 
    
 
    
class XBeachModelSetup():
    '''
    XBeach model setup class
    ''' 
    
    
    def __init__(self,fname):
        self.fname      = fname
        self.wbctype    = None
        self.wavemodel  = None
        
    def __repr__(self):
        return self.fname
    
    
    def set_params(self,input_par_dict):
        
        ## set wavemodel
        if 'Wavemodel' not in input_par_dict:
            print('no wavemodel defined. Set to Surfbeat')
            self.wavemodel = 'surfbeat'
        else:
            self.wavemodel = input_par_dict['Wavemodel']
        

        ## load parameters and categories
        f           = open(os.path.join(os.path.dirname(__file__), 'par.json'),'r')
        par_dict    = json.loads(f.read())
        
        self.input_par = {}
        ## loop over categories
        for par_category in par_dict:
            ## loop over input parameters 
            for input_par in input_par_dict:
                ## if input parameter is in category, add parameter
                if input_par in par_dict[par_category]:
                    ## create category if not exist
                    if not par_category in self.input_par:
                        self.input_par[par_category] = {}
                    ## add parameter and value                    
                    self.input_par[par_category][input_par] = input_par_dict[input_par]
        
    
    def set_grid(self,xgr,ygr,zgr, posdwn=1, xori=0, yori=0, thetamin=-90, thetamax = 90, dtheta=10):
        
        ## 1D model
        if ygr is None or ygr.shape[0]==1:
            self.ygr = None
            ## reduce size 
            self.xgr = np.reshape(xgr,len(np.squeeze(xgr) ))
            self.zgr = np.reshape(zgr,len(np.squeeze(zgr)))
        ## 2D model
        else:
            self.ygr = ygr
            self.xgr = xgr
            self.zgr = zgr
        
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
            

        
    def set_tide(self):
        pass
        
    def load_model_setup(self,path):
        ## todo
        pass    

        
        
    def write_model(self, path, figure=True):
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
            f.write('yfile\t= y.grd\n'.expandtabs(tabnumber))
            f.write('depfile\t= bed.dep\n'.expandtabs(tabnumber))
            f.write('thetamin\t= {}\n'.format(self.thetamin).expandtabs(tabnumber))
            f.write('thetamax\t= {}\n'.format(self.thetamax).expandtabs(tabnumber))
            f.write('dtheta\t= {}\n'.format(self.dtheta).expandtabs(tabnumber))
            f.write('\n')
            
            ## 
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
            ## 
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
            for ii in range(self.ny):
                for jj in range(self.nx):
                    f.write('{} '.format(self.xgr[ii,jj]))
                f.write('\n')
        ## write grid x
        with open(os.path.join(path,'y.grd'),'w') as f:
            for ii in range(self.ny):
                for jj in range(self.nx):
                    f.write('{} '.format(self.ygr[ii,jj]))
                f.write('\n')                    
       ## write dep
        with open(os.path.join(path,'bed.dep'),'w') as f:
            for ii in range(self.ny):
                for jj in range(self.nx):
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
            plt.ylabel('$YD$')
            plt.xlabel('Time')
            if save_path!=None:
                plt.savefig(os.path.join(save_path,'jonstable.png'))
        elif self.wbctype=='jons':
            print('wbctype=jons cannot be plotted')
        else:
            print('Not possible to plot wave boundary')
            
    def _plotdomain(self,save_path=None):
        plt.figure()
        if self.fast1D==True:
            plt.plot(self.xgr,self.zgr)
            plt.xlabel('x')
            plt.ylabel('z')
        else:
            plt.pcolor(self.xgr,self.ygr,self.zgr)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
        plt.grid('on')
        plt.title(self.fname)
        if save_path!=None:
            plt.savefig(os.path.join(save_path,'domain.png'))