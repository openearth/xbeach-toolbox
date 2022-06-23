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
    if len(y) == 1:
        print('1D model')
        ygr = np.linspace(0,1,1) * dymin
    else:
        if dymin==dymax:
            ygr = np.arange(np.min(y),np.max(y),dxmax)
    #else:
    #    if area_type == 'center':
    return ygr

def grid_expand_ydir(xgr,ygr,zgr,nyx,side = 'both'):
    '''
    grid_expand_ydir(xgr,ygr,zgr,nyx,side = 'both')
    adds extra grid cells (nyx ones) in the y-direction based on the y spacing of the
    last or first ycells. 
    side = 'both' or 'left' (y=0) or 'right' (y=-1)
    returns expanded grid where the original grid has kept its coordinates
    
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
    if (side=='both') | (side=='left'):
        #left side 
        dy = yl[1,0]-yl[0,0]
        yex = np.arange(nyx)*dy
        yex = yex - np.max(yex)-dy
        xex = xl[0,:]
        Xex0,Yex0 = np.meshgrid(xex,yex)
        xl = np.vstack([Xex0, xl])
        yl = np.vstack([Yex0, yl])
        zgr = np.vstack([np.tile(zgr[0,:],(nyx,1)), zgr])
    if (side=='both') | (side=='right'):        
        #end side 
        dy = yl[-1,0]-yl[-2,0]
        yex = np.arange(nyx)*dy
        yex = yex + np.max(yl) + dy
        xex = xl[-1,:]
        Xex1,Yex1 = np.meshgrid(xex,yex)
        xl = np.vstack([xl, Xex1])
        yl = np.vstack([yl, Yex1])
        zgr = np.vstack([zgr, np.tile(zgr[-1,:],(nyx,1))])
    #back to world coordinates
    xnew,ynew = rotate_grid(xl,yl,-alpha)
    
    #add origin again
    xgr2 = xnew + x0
    ygr2 = ynew + y0
    
    return xgr2, ygr2,zgr

def grid_expand_xdir(xgr,ygr,zgr,nx):
    '''
    grid_expand_xdir(xgr,ygr,zgr,nx)
    adds extra grid cells (nx ones) in the offshore x-direction based on the x spacing of the
    last or first xcells. 
    returns expanded grid where the original grid has kept its coordinates

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

    #left side 
    dx = xl[0,1]-xl[0,0]
    xex = np.arange(nx)*dx
    xex = xex - np.max(xex)-dx
    yex = yl[:,0]
    Xex0,Yex0 = np.meshgrid(xex,yex)
    
    xl = np.hstack([Xex0, xl])
    yl = np.hstack([Yex0, yl])
    zgr = np.hstack([np.tile(zgr[:,0],(nx,1)).T, zgr])
 
    #back to world coordinates
    xnew,ynew = rotate_grid(xl,yl,-alpha)
    
    #add origin again
    xgr2 = xnew + x0
    ygr2 = ynew + y0
    
    return xgr2, ygr2,zgr

def rotate_grid(xgr,ygr,theta):
    '''
    rotate_grid(xgr,ygr,theta)
    rotates a grid xgr,ygr over the angle theta (in radians)
    
    Author: Marlies van der Lugt
    Revision 0     
    '''
    ny,nx = xgr.shape
    coords = np.vstack((xgr.reshape(-1),ygr.reshape(-1)))
    rotMatrix = np.array([[np.cos(theta), np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    uv,vd = rotMatrix @ coords
    uv = uv.reshape([ny,nx])
    vd = vd.reshape([ny,nx])
    return uv, vd

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