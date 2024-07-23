"""
Created on Wed May 5 10:04:00 2023

@author: Menno de Ridder, Marlies van der Lugt
collection containing wave functions
"""
import numpy as np
import numpy as np
from scipy.special import gamma
from scipy import signal
import matplotlib.pyplot as plt

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