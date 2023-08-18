# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:25:00 2023
@author: Cas van Bemmelen
collection grid checks
"""

from matplotlib import pyplot as plt
import numpy as np
import os 

def diff(A,direc = 0):
    '''
    Computes the difference in subsequent items
    A is either a vector or a 2D matrix
    if A is a matrix, 
        direc = 0 performs operation horizontally
        direc = 1 performs operation vertically
    Author: Cas van Bemmelen
    Date: 26-06-2023
    
    Parameters
    ----------
    A : TYPE
        DESCRIPTION.
    direc : TYPE, optional
        DESCRIPTION. The default is 0.
    
    Returns
    -------
    dif : TYPE
        DESCRIPTION.
    
    '''
    # Make sure A is a numpy array
    A = np.array(A)

    # Check if A is a vector
    if A.ndim == 1:
        
        # Prealocate
        dif = A * 0;
        
        # Loop over all colums and rows to substract subsequent elements
        for ii in range(1,len(A)):
            dif[ii] = A[ii] - A[ii-1]
        
        return dif
    
    # Check A is a 2d matrix        
    if A.ndim == 2:
        
        # Prealocate
        dif = A * 0;
        
        # Loop over all colums to substract subsequent elements
        if direc == 1:
            for ii in range (1,A.shape[1]):
                dif[:,ii] = A[:,ii] - A[:,ii-1]
        
        # Loop over all rows to substract subsequent elements
        elif direc == 0:
            for ii in range(1,A.shape[0]):
                dif[ii,:] = A[ii,:] - A[ii-1,:]
        
        return dif

def smoothness(x, y, outputloc = ''):
    '''
    Perform grid checks based on D3D-quickin procedure for 
    smoothness in m and n-direction and stores figures
    note that it is advisable to use local coordinate system!
    Author: Cas van Bemmelen
    Date: 26-06-2023
    
    Parameters
    ----------
    gridname : TYPE
        DESCRIPTION.
    outputloc : TYPE, optional
        DESCRIPTION. The default is ''.
    
    Returns
    -------
    None.
    
    '''
    
    plt.close('all')
    
    #N-Smoothness
    
    #First import grid  and remove origin of the grid (in some cases this is defined at [0,0] [lon,lat])
    
    # Determine dimensions grid
    m = int(x.shape[0])
    n = int(x.shape[1])
    
    # Preallocate
    SMO_emp = np.empty((m,n))
    SMO_emp[:,:] = np.nan
    
    # Compute smoothness (using created diff() function)
    SMO_emp[0:m,:]      = np.sqrt(diff(x)**2 + diff(y)**2)
    SMO_emp[1:m-1,:]    = np.divide(SMO_emp[0:m-2,:],SMO_emp[1:m-1,:])
    SMO_empdiv          = np.divide(np.ones((m,n)),SMO_emp)
    # Remove large values due to the division step that is performed above
    SMO_empdiv[SMO_empdiv>100] = float('NaN')
    # Finlize smoothness and set nan values at zero
    smooth_matrix = np.maximum(SMO_emp,SMO_empdiv)
    smooth_matrix[np.isnan(smooth_matrix)] = 0
    
    # Plot and store figure
    fig1 = plt.figure()
    pcol = plt.pcolor(x[1:-1,1:-1],y[1:-1,1:-1],smooth_matrix[1:-1,1:-1],cmap = 'jet',vmin = 1-0.15/4,vmax = 1.15+0.15/4);
    plt.colorbar(pcol)
    plt.axis('equal')
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.axis('scaled')
    plt.title('n-smoothness')
    fig1.savefig(os.path.join(outputloc, 'n_smoothness.png'), dpi=300)
    
    #M-Smoothness   
    xx = x.transpose()
    yy = y.transpose()
    
    # Determine dimensions grid    
    m = int(xx.shape[0])
    n = int(xx.shape[1])

    # Prealocate
    SMO_emp = np.empty((m,n))
    SMO_emp[:,:] = np.nan

    # Compute smoothness (using created diff() function)
    SMO_emp[0:m,:]      = np.sqrt(diff(xx)**2 + diff(yy)**2)
    SMO_emp[1:m-1,:]    = np.divide(SMO_emp[0:m-2,:],SMO_emp[1:m-1,:])
    SMO_empdiv          = np.divide(np.ones((m,n)),SMO_emp)
    # Remove large values due to the division step that is performed above
    SMO_empdiv[SMO_empdiv>100] = float('NaN')
    # Finlize smoothness and set nan values at zero
    smooth_matrix = np.maximum(SMO_emp,SMO_empdiv)
    smooth_matrix[np.isnan(smooth_matrix)] = 0
    
    # Plot and store figure
    fig2 = plt.figure()
    pcol = plt.pcolor(xx[1:-1,1:-1],yy[1:-1,1:-1],smooth_matrix[1:-1,1:-1],cmap = 'jet',vmin = 1-0.15/4,vmax = 1.15+0.15/4);
    plt.colorbar(pcol)
    plt.axis('equal')
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.axis('scaled')
    plt.title('m-smoothness')
    fig2.savefig(os.path.join(outputloc, 'm_smoothness.png'), dpi=300)

def aspect_ratio(x, y, outputloc = ''):    
    '''
    Perform grid checks based on D3D-quickin procedure for 
    aspect ratio and stores figures
    note that it is advisable to use local coordinate system!
    Author: Cas van Bemmelen
    Date: 26-06-2023
    
    Parameters
    ----------
    gridname : TYPE
        DESCRIPTION.
    outputloc : TYPE, optional
        DESCRIPTION. The default is ''.
    
    Returns
    -------
    None.
    
    '''  
    plt.close('all')

    # Determine dimensions grid
    m = int(x.shape[0])
    n = int(x.shape[1])
    
    # Prealocate
    asp_ratio = np.empty((m,n))
    
    # Loop over all celss to determine aspect ratio
    # Remark: could be optimized by using list comprehension
    for mm in range(1,m):
        for nn in range(1,n):
            asp_ratio[mm,nn] = (max(np.sqrt((x[mm,nn]-x[mm-1,nn])**2 + (y[mm,nn]-y[mm-1,nn])**2), np.sqrt(((x[mm,nn]-x[mm,nn-1])**2) + ((y[mm,nn]-y[mm,nn-1])**2))))/(min(np.sqrt(((x[mm,nn]-x[mm-1,nn])**2) + ((y[mm,nn]-y[mm-1,nn])**2)), np.sqrt(((x[mm,nn]-x[mm,nn-1])**2) + ((y[mm,nn]-y[mm,nn-1])**2)))) 
    
    # Plot and save figure
    fig1 = plt.figure()
    pcol = plt.pcolor(x[1:,1:],y[1:,1:],asp_ratio[1:,1:],cmap = 'jet',vmin = 0.5,vmax = 3.5);
    plt.colorbar(pcol)
    plt.axis('scaled')
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.title('aspect ratio')
    fig1.savefig(os.path.join(outputloc, 'aspectratio.png'), dpi=300)
    

def orthogonality(x, y, outputloc = ''):  
    '''
    Perform grid checks based on D3D-quickin procedure for 
    orthogonality and stores figures
    note that it is advisable to use local coordinate system!
    Author: Cas van Bemmelen
    Date: 26-06-2023
    
    Parameters
    ----------
    gridname : TYPE
        DESCRIPTION.
    outputloc : TYPE, optional
        DESCRIPTION. The default is ''.
    
    Returns
    -------
    None.
    
    '''  
    plt.close('all')
    
    # Determine dimensions grid
    m = int(x.shape[0])
    n = int(x.shape[1])
    
    # Prealocate
    dx1 = np.empty((m,n))
    dy1 = np.empty((m,n))
    dx2 = np.empty((m,n))
    dy2 = np.empty((m,n))
    
    # Compute orthogonality
    dx1[1:m-1,:] = x[2:m,:]-x[0:m-2,:]
    dy1[1:m-1,:] = y[2:m,:]-y[0:m-2,:]
    
    dx2[:,1:n-1] = x[:,2:n]-x[:,0:n-2]
    dy2[:,1:n-1] = y[:,2:n]-y[:,0:n-2]
    
    ds1sq = (dx1**2) + (dy1**2)
    ds2sq = (dx2**2) + (dy2**2)
    ds3sq = ((dx1+dx2)**2) + ((dy1+dy2)**2)
    ortho = np.absolute(np.divide((-ds3sq+ds1sq+ds2sq),(2*np.sqrt(ds1sq * ds2sq))))

    # Plot and save figure
    fig1 = plt.figure(figsize=(16,6))
    pcol = plt.pcolor(x,y,ortho,cmap = 'jet',vmin = 0,vmax = 0.05);
    plt.colorbar(pcol)
    plt.axis('scaled')
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.title('orthogonality')
    fig1.savefig(os.path.join(outputloc, 'orthogonality.png'), dpi=300)
