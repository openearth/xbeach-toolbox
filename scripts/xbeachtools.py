import numpy as np


def dispersion(w,d):
    '''
    Compute wave number for given w and d
    '''
    g   = 9.81
    k1  = 1000
    for i in range(100000):
        k       = k1
        k1      = w**2.0/g/np.tanh(k * d)
        error   = abs(k1-k)
        if error < 0.00001:
            break   
    if i==99999:
        print ('Warning: no convergence')
    return k 

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
    
def ygrid(y,
           dymin = 10,
           dymax = 20,
           area_type='center'):
     '''
     '''
     
    if len(y)==1:
        print('1D model')
        ygr = np.linspace(0,1,1) * dymin
    else:
        if dymin==dymax:
            ygr = np.arange(np.min(y),np.max(y),dxmax)
        else:
            if area_type == 'center':
                