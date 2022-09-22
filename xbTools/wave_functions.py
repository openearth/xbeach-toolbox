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
        gravitational acceleration

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
    Computes the group velocity and wave celerity ratio based on the period and water depth

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
    determine mean wave direction and directional spreading, Kuik et al. 1988
    see SWAN user manual for definitions
    theta0 and dspr only for range [fmin,fmax]

    Parameters
    ----------
    theta : array
        DESCRIPTION.
    ee : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        directional spreading in degrees.

    '''

    dtheta = theta[1] - theta[0]
    st = np.sin(theta)
    ct = np.cos(theta)
    m0Sin = np.dot(ee, st) * dtheta
    m0Cos = np.dot(ee, ct) * dtheta
    m0 = np.sum(ee)*len(theta)*dtheta
    dspr2 = 2 * (1 - np.sqrt((m0Sin / m0) ** 2 + (m0Cos / m0) ** 2))
    return np.rad2deg(np.sqrt(dspr2))  # directional spreading in degrees

