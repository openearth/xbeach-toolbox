import numpy as np



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


def wavecelerity(Tp, d, g=9.81):
    '''
    

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
    
    k = dispersion(2*np.pi/Tp, d, g)
    n = 0.5 * (1 + 2 * k * d * np.sinh(2 * k * d))
    c = g * Tp/(2 * np.pi) * np.tanh(k * d)
    cg = n * c
    return cg, n


