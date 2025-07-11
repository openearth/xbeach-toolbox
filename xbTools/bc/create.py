import os
import numpy as np
import pandas as pd

def construct_bc_for_32Hr_storm(zsmax, Hs, Tp, dir, tide_amp, s=10, dt=1800, gamma=3.3, dtbc=1):
    '''
    function to generate jonstable and tide file for 32 hr storm superimposed on 2.5 tides
    immediately writes these to file in the folder 'fold'

    Parameters
    ----------
    zsmax: float
        maximum surge level
    Hs: float
        significant wave height
    Tp: float
        peak period
    dir: float
        wave direction
    tide_amp: float
        amplitude of the tide

    Returns
    ----------
    HR : dict
        dictionary containing timeseries of time varying boundary conditions
    tide : 2D array
        array containing timeseries of water level variations offshore: [t, zs]
    '''

    #gaussian shaped Hs buildup over 32 hours
    x = np.arange(-4, 4+0.05, 0.05)
    y = np.exp(-0.5 * x**2) / (np.sqrt(2*np.pi))
    y = y/np.max(y)
    ix1 = np.where(y>1/6)[0][1]
    ix2 = len(x)-ix1+1
    x = x[ix1:ix2] 
    y = y[ix1:ix2] 
    x = x-np.min(x); x = x/np.max(x); x = 32*x

    HR = {}
    HR['t'] = np.arange(0, 32.5, dt/3600)
    HR['Hm0'] = np.interp(HR['t'], x, Hs*y)

    #matching wave period in which steepness remains constant
    steep = 2*np.pi*Hs/9.8/Tp**2
    HR['Tp'] = np.sqrt(2*np.pi*HR['Hm0']/9.8/steep)

    # generate zs: 2.5 tides and evolution of maximum surge on top
    t2 = np.arange(0, 12.9, 0.1) 
    wl2 = tide_amp*np.sin(2*np.pi/12.8*t2)
    n_wl2 = len(wl2)
    t2 = np.arange(0, 32.1, 0.1)
    wl2 = np.append(np.append(wl2, wl2[1:]), wl2[1:int(n_wl2/2)+1])
    wl3 = np.interp(HR['t'], t2, wl2) + np.interp(HR['t'], x, zsmax*y)
    wl3 = wl3-np.min(wl3)
    wl3 = wl3/np.max(wl3)*zsmax
    # HR['zs'] = wl3
    tide = np.vstack([HR['t']*3600, wl3]).T

    # stationary parameters of jonstable
    HR['mainang'] = dir+0*HR['Hm0']
    HR['duration'] = dt+0*HR['Hm0']
    HR['gammajsp'] = gamma+0*HR['Hm0']
    HR['s'] = s+0*HR['Hm0']
    HR['dtbc'] = dtbc+0*HR['Hm0']

    return HR, tide
       
