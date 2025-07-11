"""
Created on Wed May 5 10:04:00 2023

@author: Menno de Ridder, Marlies van der Lugt
collection containing wave functions
"""
import numpy as np
import numpy as np
from scipy.special import gamma
from scipy import signal
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt

def dispersion(w, d, max_error=0.00001,g=9.81):
    '''
    Computes the wave number given a radial frequency and water depth

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

def disper(w, h, g=9.8):
    '''
    Computes the wave number given a radial frequeny and water depth using an approximate, 
    noniterative method to solve the linear dispersion relationship.
    absolute error in k*h < 5.0e-16 for all k*h
    
    Syntax:
    k = disper(w, h, [g])
    
    Input:
    w = 2*pi/T, were T is wave period (either float or array)
    h = water depth
    g = gravity constant
    
    Output:
    k = wave number
    
    Example
    k = disper(2*pi/5,5,g = 9.81);
    
    Copyright notice
    --------------------------------------------------------------------
    Copyright (C) 
    G. Klopman, Delft Hydraulics, 6 Dec 1994
    M. van der Lugt conversion to python, 11 Jan 2021
    
    '''    
    #make sure numpy array
    listType = type([1,2])
    Type = type(w)

    w = np.atleast_1d(w)
    
    #check to see if warning disappears
    wNul = w==0
    w[w==0] = np.nan

    
    w2 = w**2*h/g
    q = w2 / (1-np.exp(-w2**(5/4)))**(2/5)
    
    for j in np.arange(0,2):
        thq = np.tanh(q)
        thq2 = 1-thq**2
        aa = (1 - q*thq) * thq2
        
        #prevent warnings, we don't apply aa<0 anyway
        aa[aa<0] = 0
        
        bb = thq + q*thq2
        cc = q*thq - w2
        
        
        D = bb**2-4*aa*cc
        
        # initialize argument with the exception
        arg = -cc/bb
        
        # only execute operation on those entries where no division by 0 
        ix = np.abs(aa*cc)>=1e-8*bb**2 
        arg[ix] = (-bb[ix]+np.sqrt(D[ix]))/(2*aa[ix]) 

                
        q = q + arg

              
    k = np.sign(w)*q/h
    
    #set 0 back to 0
    k = np.where(wNul,0,k)

    #if input was a list, return also as list
    if Type==listType:
        k = list(k)
    elif len(k)==1:
        k = k[0]
        
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
    
    k       = dispersion(2*np.pi/Tp, d, max_error = 0.0001, g = g)
    n       = 0.5 * (1 + 2 * k * d / np.sinh(2 * k * d))
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


def jonswap2(f, Hm0, Tp, g=9.8, gamma=3.3, method='Yamaguchi', normalize=True, eps=1e-6, sa=0.07, sb=0.09, pdir=None, ms=None, directions=None):
    '''
        
        Computes jonswap spectrum based on a given Hm0 and Tp and frequency grid. Optionally creates 2D directional spectrum too
        
        Parameters
        ----------
        f   : array
            frequency axis.
        Hm0 : float
            Wave height.
        Tp : float
            Peak period.
        g : float
            gravitational acceleration.
        gamma : float
            jonswap form factor
        method : string
            switch between Yamaguchi and Goda approaches to compute alpha
        normalize : logical
            switch to enforce variance to be equivalent to Hm0
        eps : float
            small number to prevent division by zero
        sa : float
        sb: float
        pdir : float, optional
            peak wave direction
        ms : float
            directional spreading factor
        directions : array
            desired directional axis

        Returns
        -------
        jon : array
            returns either 1D or 2D jonswap spectrum

        '''
    
    # Pierson-Moskowitz
    if method=='Yamaguchi':
        alpha = 1/(0.06533*gamma**0.8015+0.13467)/16  # used in SWAN
    elif method=='Goda':
        alpha = 1/(0.23+0.03*gamma-0.185*(1.9+gamma)**(-1)/16)
    else:
        print('method unkown, choose either of Yamaguchi or Goda')
        return
    
    pm = alpha*Hm0**2*Tp**(-4)*f**(-5)*np.exp(-1.25*(Tp*f)**(-4))
    
    # apply JONSWAP shape
    jon = pm*gamma**(np.exp(-0.5*(Tp*f-1)**2/(_sigma(f, 1/Tp, sa, sb))**2))

    jon[np.isnan(jon)] = 0

    # optionally correct total energy of user-discretized spectrum to match Hm0, as above methods are only an approximation
    if normalize:
        corr = Hm0**2/(16*np.trapz(jon, f)+eps)
        jon = jon*corr

    #2D
    if pdir is None and ms is None and directions is None:
        return jon
    else:
        assert pdir is not None and ms is not None and directions is not None, 'for 2D spectrum prescribe all three of pdir, ms and directions'

        if len(np.atleast_1d(ms))==1:
            ms = np.array([ms]*len(f))

        jon2 = np.tile(jon, [len(directions), 1])
        for iff in range(len(f)):
            cdir = _directional_spreading(directions, pdir, ms[iff])
            jon2[:, iff] = cdir*jon2[:, iff]

        if normalize:
            corr = Hm0**2 / (16*np.trapz(np.trapz(jon2, directions), f)+eps)
            jon2 = jon2*corr

        return jon2


def _sigma(f,fpeak,sa,sb):
    '''
    auxiliary function for jonswap spectrum computation
    '''
    s = np.array([sa]*len(f))
    s[f > fpeak] = sb
    return s

def _directional_spreading(directions, pdir, ms, units='deg', plot=0, quiet=0):
    '''
    auxiliary function for directional spreading of 2D jonswap spectrum

    '''
    if units=='deg':
        directions = directions/180*np.pi
        pdir = pdir/180*np.pi

    A1 = 2**ms * (gamma(ms/2+1))**2 / (np.pi * gamma(ms+1)) # coefficient to make sure integral is 1

    cdir = 0*directions
    for id in range(len(directions)):
        acos = np.cos(directions[id]-pdir)
        if acos > 0:
            cdir[id] = A1 * np.maximum(acos**ms, 1e-10)

    if units=='deg':
        cdir = cdir*np.pi/180
        directions = directions*180/np.pi
        pdir = pdir*180/np.pi

    # check for ill-sampling
    int = np.trapz(cdir, directions)
    if not(np.abs(int-1)<1e-6) and not(quiet):
        print('integral not 1 for ms={}:{}'.format(ms, int))

    if plot:
        plt.figure()
        plt.plot(directions, cdir)
        plt.axvline(pdir)


    return cdir     


def spectrum_simple(x,y,
                    fresolution=0.01,
                    Nr=None,
                    detrend=True,
                    overlap=0.5,
                    windowFilter='Hamming',
                    tolerance = 1e-3,
                    strict = False,
                    correctvar = True):   
    """
    spectrum_simple(sf,zs) 
    wrapper function for fast fourier transform that takes care of windowing of the inputsignal to obtain a power spectrum at the desired resolution in frequency space.
    
    Parameters
    ----------
    x : FLOAT
        DATETIME64 TIME ARRAY OR SAMPLING FREQUENCY (WHEN ONLY ONE VALUE).
    y : FLOAT
        SURFACE ELEVATION [M].
    fresolution : FLOAT, optional
        DESIRED RESOLUTION IN THE FREQUENCY DOMAIN. The default is 0.01.
    detrend : BOOL, optional
        DETREND THE SIGNAL YES/NO. The default is True.
    overlap : FLOAT, optional
        OVERLAPPING PERCENTAGE OF THE WINDOWS. The default is 0.5.
    windowFilter : STRING, optional
        WINDOW TYPE. The default is 'Hamming'.
    tolerance : FLOAT, optional
        WARNING TOLLERANCE. The default is 1e-3.
    correctvar : BOOL, optional
        RESTORE TOTAL VARIANCE IN FREQ SPACE YES/NO. The default is True.

    Returns
    -------
    fx, vy
    fx = frequency axis [Hz]
    vy = power density [m2/Hz]
    
    Matlab to Python: Marlies van der Lugt 14-01-2020

    """
    import numpy as np
    from scipy.fft import fft, ifft
    #input checks differentiating for datetime64 and seconds
    if len(np.atleast_1d(x))==1:
        dt = 1/x # 1/frequency
    else:
        try:
            if ((np.min(np.diff(x))/np.timedelta64(1,'s')<0) or 
                ((np.max(np.diff(x))-np.min(np.diff(x)))/np.timedelta64(1,'s')>1e-8 )):   
                print('Input xx must be monotonic increasing')
                exit
            else:
                x = (x - x[0])/np.timedelta64(1,'s')
        except:
            if ((np.min(np.diff(x))<0) or 
            ((np.max(np.diff(x))-np.min(np.diff(x)))>1e-8 )):   
                print('Input xx must be monotonic increasing')
            else:
                x = x-x[0]
        
        # Time step in time/space axis
        dt= np.mean(np.diff(x))

    # Number of data points in the total signal
    N = len(y)   
    if Nr==None:
        # Number of data points required for desired fresolution      
        Nr = int(np.ceil(1/fresolution/dt))
        if strict:
            #TODO: make a nextpow2 function
            Nr = Nr
    
        # Check if there are sufficient data points to acquire the set frequency
        # resolution
        if Nr>N:
            # reset Nr
            Nr = N
            print('Warning: The required frequency resolution could not be achieved.')

    # Number of Welch repetitions
    Nw = int(np.ceil((N-Nr)/(overlap*Nr))+1)

    # Index input arrays for start and end of each Welch repetition
    indend = np.linspace(Nr,N,Nw).astype(int)
    indstart = (indend-Nr).astype(int)

    # Time and frequency array 
    T  = dt*Nr
    df = 1/T
    ffa = np.arange(0,np.round(Nr/2))
    ffb = -np.arange(np.round(Nr/2),0,-1)
    ff = df*np.append(ffa,ffb)
    fx = ff[0:int(Nr/2)]

    # Storage arrays
    vy = np.zeros([int(Nr/2)])

    # % Detrend input signal
    if  detrend:
        # pdb.set_trace()
        y = signal.detrend(y)
        
    varpre = np.var(y)
    
    # do Welch repetition
    for i in np.arange(0,Nw):
        d = y[indstart[i]:indend[i]]
        if (windowFilter == 'Hann'):
            d = d * 0.5*(1-np.cos(2*np.pi*np.arange(0,Nr)/(Nr-1)))
            varpost = np.var(d)            
        elif (windowFilter == 'Hamming'):
            d = d * (0.54-0.46*np.cos(2*np.pi*np.arange(0,Nr)/(Nr-1)))
            varpost = np.var(d)            
        elif (windowFilter == None):
            varpost = varpre

        # FFT
        Q = fft(d)
        # Power density
        V = 2*T/Nr**2*np.abs(Q)**2
        # Store in array
        vy = vy + 1/Nw*V[0:int(Nr/2)]


    # Restore the total variance
    if correctvar:
       vy = (varpre/np.trapz(vy,dx = df))*vy


    # input/output check
    hrmsin = 4*np.std(y)/np.sqrt(2)
    hrmsout =  np.sqrt(8*np.trapz(vy, dx = df))
    dif = np.abs(hrmsout-hrmsin)
    if (dif > tolerance):
        print('Warning: Difference in input and output wave height ({}) is greater than set tolerance ({})'.format(dif,tolerance))

    return fx, vy


def Guza_split_waves(t, zsi, umi, zb, boundopt='free', quietopt=False):
    """
    Guza_split_waves: Decomposes measured water level and velocity data into incident and reflected wave components.
    
    This function processes the total surface elevation and depth-averaged current measurements to separate them into
    incident and reflected components based on the chosen boundary condition option. It can operate in time space or
    Fourier space, depending on the specified boundary option, and can apply a high-frequency cutoff to filter noise.
    
    Parameters:
    - t: array-like
        Time vector in seconds.
    - zsi: array-like
        Total surface elevation vector in meters (relative to SWL).
    - umi: array-like
        Depth-averaged current vector in meters per second.
    - zb: float
        Bed level at the location of zs and um in meters (relative to SWL).
    - boundopt: str, optional
        Boundary option for wave separation. Options include:
        'boundin', 'boundout', 'boundall', 'free', 'boundupper', 'sqrt(gh)'.
        Default is 'free'.
    - quietopt: bool, optional
        If True, suppresses display messages. Default is False.
    
    Returns:
    - zsin: array-like
        Incident wave surface elevation vector in meters.
    - zsout: array-like
        Reflected wave surface elevation vector in meters.
    - uin: array-like
        Incident wave depth-averaged current vector in meters per second.
    - uout: array-like
        Reflected wave depth-averaged current vector in meters per second.
    
    Notes:
    - The function calculates the mean water level and current, adjusts the signals to zero-centered, and optionally
      detrends the data.
    - For the 'sqrt(gh)' option, the function computes the components directly in the time domain. For other options,
      it performs the decomposition in the Fourier domain, applying a high-frequency cutoff to isolate noise.
    - It includes optional plotting to visualize the separated components and the original measurements.
    
    Example usage:
    >>> t = np.linspace(0, 100, 1001)
    >>> zsi = np.sin(2 * np.pi * t / 10) + np.random.randn(len(t)) * 0.1
    >>> umi = np.cos(2 * np.pi * t / 10) + np.random.randn(len(t)) * 0.1
    >>> zb = 0
    >>> zsin, zsout, uin, uout = Guza_split_waves(t, zsi, umi, zb, boundopt='sqrt(gh)', quietopt=True)
    """
        
    # Constants
    g = 9.81

    # Mean water level and velocity
    zsm = np.mean(zsi)
    umm = np.mean(umi)
    
    if not quietopt:
        print(f"Mean water level found is : {zsm:.2f}m")
        print(f"Mean velocity found is    : {umm:.2f}m/s")
    
    # Average water depth
    h = zsm - zb

    # Adjust to zero-centered water level and velocity
    zs = zsi - zsm
    um = umi - umm

    # Detrend signals
    zsd = signal.detrend(zs)
    umd = signal.detrend(um)
    zsm += (zs - zsd)
    zs = zsd
    umm += (um - umd)
    um = umd

    if boundopt == 'sqrt(gh)':
        # In time space
        hh = zs + h
        c = np.sqrt(9.81 * hh)
        q = umd * hh
        ein = (zs * c + q) / (2 * c)
        eout = (zs * c - q) / (2 * c)
        zsin = ein + zsm
        zsout = eout + zsm
        uin = (np.sqrt(1. / hh**2) * c * ein) + umm
        uout = -(np.sqrt(1. / hh**2) * c * eout) + umm
    else:
        # Into Fourier space
        n = len(zs)
        Z = fft(zs, n)
        U = fft(um, n)
        df = 1 / (t[-1] - t[0])
        ff = df * np.concatenate((np.arange(0, np.ceil(n / 2)), np.arange(-np.floor(n / 2), 0)))
        w = 2 * np.pi * ff
        w[0] = 0
        k = disper(w, np.mean(h), g)  
        c = w / k
        c[0] = 0

        # hf cutoff
        filterfreq = 0.01
        minfreq = 0.03
        incl = ff >= minfreq
        ftemp = ff[incl]
        vartemp = np.real(Z)**2
        vartemp = vartemp[incl]
        window = max(1, round(filterfreq / df))
        vartemp = filterwindow(vartemp, window)  
        fp = ftemp[np.argmax(vartemp)]
        Tp = 1 / fp
        hfc = min(10 * fp, max(ff))
        
        if not quietopt:
            print(f"Peak period found is      : {Tp:.2f}s")
        
        # Cutoff high frequency c
        cg, n = wavecelerity(1 / hfc, h)
        c = np.maximum(c, cg/n)  
        
        # Properties of peak waves
        cgp, n = wavecelerity(Tp, h) 
        cp = cgp/n

        # Select frequencies for bound components
        freein = freeout = boundin = boundout = None
        if boundopt == 'boundin':
            freein = np.abs(ff) > 0.5 * fp
            boundin = np.abs(ff) <= 0.5 * fp
            freeout = np.ones_like(ff, dtype=bool)
            boundout = np.zeros_like(ff, dtype=bool)
        elif boundopt == 'boundout':
            freein = np.ones_like(ff, dtype=bool)
            boundin = np.zeros_like(ff, dtype=bool)
            freeout = np.abs(ff) > 0.5 * fp
            boundout = np.abs(ff) <= 0.5 * fp
        elif boundopt == 'boundall':
            freein = freeout = np.abs(ff) > 0.5 * fp
            boundin = boundout = np.abs(ff) <= 0.5 * fp
        elif boundopt == 'free':
            freein = freeout = np.ones_like(ff, dtype=bool)
            boundin = boundout = np.zeros_like(ff, dtype=bool)
        elif boundopt == 'boundupper':
            fac = 2
            freein = freeout = np.abs(ff) < fac * fp
            boundin = boundout = np.abs(ff) >= fac * fp
        
        # Find the velocity for all Fourier components
        cin = cout = np.zeros_like(ff)
        cin[freein] = c[freein]
        cout[freeout] = c[freeout]
        cout[boundout] = cgp
        if boundopt == 'boundupper':
            cin[boundin] = cout[boundout] = cp
        else:
            cin[boundin] = cout[boundout] = cgp

        # Maximize to the long wave celerity
        cin = np.minimum(cin, np.sqrt(9.81 * h))
        cout = np.minimum(cout, np.sqrt(9.81 * h))

        # Cut off hf noise
        mm = round(hfc / df)
        set = np.zeros(len(Z))
        set[1:mm] = 1
        set[-mm+1:] = 1

        ein = np.real(ifft(set * ((Z * cout + U * h) / (cin + cout))))
        eout = np.real(ifft(set * ((Z * cout - U * h) / (cin + cout))))
        zsin = ein + zsm
        zsout = eout + zsm
        uin = np.real(ifft(set * (np.sqrt(1. / h**2) * cin * fft(ein)))) + umm
        uout = -np.real(ifft(set * (np.sqrt(1. / h**2) * cout * fft(eout)))) + umm

    reflc = np.std(zsout)**2 / np.std(zsin)**2

    if not quietopt:
        print(f"Energy reflection found is: {reflc:.2f} [-]")
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(t, zsi, 'r', t, zsin, 'b', t, zsout, 'b--', t, zsin + zsout - zsm, 'g--', t, zsm, 'k--')
        plt.title('Water level')
        plt.legend(['zs measured', 'zsin', 'zsout', 'zsin+zsout', 'mean zs'])
        plt.subplot(2, 1, 2)
        plt.plot(t, umi, 'r', t, uin, 'b', t, uout, 'b--', t, uin + uout - umm, 'g--', t, umm, 'k--')
        plt.title('Velocity')
        plt.legend(['u measured', 'uin', 'uout', 'uin+uout', 'u mean'])
        plt.show()

    return zsin, zsout, uin, uout

def filterwindow(y, window, iterations=1, nanhandle=0, method='mean', direction='center'):
    """
    Applies a windowed filter to the input data 'y' using the specified method.
    
    Parameters:
    - y: array-like
        The input data to be filtered.
    - window: int
        The size of the window over which to apply the filter.
    - iterations: int, optional
        The number of times the filter is applied. Default is 1.
    - nanhandle: int, optional
        Handling of NaN values: 0 (ignore NaNs), 1 (handle NaNs). Default is 0.
    - method: str, optional
        The method of filtering: 'mean', 'min', 'max', 'var', 'std'. Default is 'mean'.
    - direction: str, optional
        The direction of the filter window: 'center', 'backwards', 'forwards'. Default is 'center'.
    
    Returns:
    - fy: numpy.ndarray
        The filtered data.
    """

    if iterations is None or iterations == 0:
        iterations = 1
    if nanhandle is None:
        nanhandle = 0
    if method is None:
        method = 'mean'
    if direction is None:
        direction = 'center'

    y = np.asarray(y)
    if y.ndim == 1:
        y = y[np.newaxis, :]  # Convert to row vector if necessary

    fy = np.zeros_like(y)
    iw = window // 2

    if direction == 'center':
        il = -iw
        ir = iw
    elif direction == 'backwards':
        il = -2 * iw
        ir = 0
    elif direction == 'forwards':
        il = 0
        ir = 2 * iw
    else:
        raise ValueError("Invalid direction value. Choose from 'center', 'backwards', or 'forwards'.")

    for _ in range(iterations):
        for i in range(y.shape[1]):
            window_slice = y[0, max(0, i + il):min(y.shape[1], i + ir)]
            if method == 'mean':
                if nanhandle == 0:
                    fy[0, i] = np.mean(window_slice)
                else:
                    fy[0, i] = np.nanmean(window_slice)
            elif method == 'min':
                if nanhandle == 0:
                    fy[0, i] = np.min(window_slice)
                else:
                    fy[0, i] = np.nanmin(window_slice)
            elif method == 'max':
                if nanhandle == 0:
                    fy[0, i] = np.max(window_slice)
                else:
                    fy[0, i] = np.nanmax(window_slice)
            elif method == 'var':
                if nanhandle == 0:
                    fy[0, i] = np.var(window_slice)
                else:
                    fy[0, i] = np.nanvar(window_slice)
            elif method == 'std':
                if nanhandle == 0:
                    fy[0, i] = np.std(window_slice)
                else:
                    fy[0, i] = np.nanstd(window_slice)
            else:
                raise ValueError("Invalid method value. Choose from 'mean', 'min', 'max', 'var', 'std'.")
        y = fy

    if y.shape[0] == 1:
        fy = fy[0, :]  # Convert back to original shape if necessary

    return fy