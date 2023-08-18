"""
Created on Wed May 5 10:04:00 2023

@author: Menno de Ridder
collection containing wave functions
"""
import numpy as np

def dispersion(w, d, max_error=0.00001,g=9.81):
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
    g : float, optional
        gravitational acceleration.

    Returns
    -------
    k : float
        wave number.

    '''

    ## initial guess
    k1  = 1000
    ## iterate untill convergence
    N = 100000
    for i in range(N):
        ## update k
        k       = k1
        ## next iteration
        k1      = w**2.0/g/np.tanh(k * d)
        ## compute error
        error   = abs(k1-k)
        if error < max_error:
            break   
    ## no convergence
    if i==N-1:
        print ('Warning: no convergence')
    return k 


def wavecelerity(Tp, d, g=9.81):
    '''
    Computes the group velocity and wave celerity ratio based on the wave period and water depth.

    Parameters
    ----------
    Tp : float
        Peak period.
    d : float
        Water depth.
    g : float, optional
        gravitational acceleration. The default is 9.81.

    Returns
    -------
    cg : float
        Group velocity.
    n : float
        celerity ratio.

    '''
    
    k       = dispersion(2*np.pi/Tp, d, g)
    n       = 0.5 * (1 + 2 * k * d * np.sinh(2 * k * d))
    c       = g * Tp/(2 * np.pi) * np.tanh(k * d)
    cg      = n * c
    return cg, n

def directional_spread_kuik(theta,ee):
    '''
    Determine mean wave direction and directional spreading, Kuik et al. 1988
    see SWAN user manual for definitions
    theta0 and dspr only for range [fmin,fmax]

    Parameters
    ----------
    theta : array
        Array with the directions.
    ee : TYPE
        Array with the energy density.

    Returns
    -------
    float
        Directional spreading in degrees.

    '''

    dtheta = theta[1] - theta[0]
    st = np.sin(theta)
    ct = np.cos(theta)
    m0Sin = np.dot(ee, st) * dtheta
    m0Cos = np.dot(ee, ct) * dtheta
    m0 = np.sum(ee)*len(theta)*dtheta
    dspr2 = 2 * (1 - np.sqrt((m0Sin / m0) ** 2 + (m0Cos / m0) ** 2))
    return np.rad2deg(np.sqrt(dspr2))  # directional spreading in degrees

def celerity_ratio_equals_09(Tp,d_start):
    '''
    Function to find water depth for which the n ration equal 0.9.

    Parameters
    ----------
    Tp : float
        Peak period.
    d_start : float
        Water depth.

    Returns
    -------
    d : float
        Depth.

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
    
    Compute required Offshore water depth to correctly force the waves
    
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
        d_n                 = celerity_ratio_equals_09(Tp,d_start)
        while np.abs(d_start-d_start_previous)>0.05:
            ## update depth
            d_start_previous = d_start
            ## compute required depth
            d_start         = np.max([3.33333*Hm0_shoal, d_n])
            ## compute hm0 shoal
            cg, n_startdepth        = wavecelerity(Tp, d_start)
            Hm0_shoal               = Hm0 * np.sqrt(cg_bc/cg)
            ## update count
            count =count+ 1
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