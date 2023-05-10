# -*- coding: utf-8 -*-
"""
Created on Wed May 5 10:04:00 2023

@author: Cas van Bemmelen
collection containing executing run functionalities
"""
# import python modules
import os
import numpy as np

from xbTools.xbeachtools import XBeachModelSetup

def xb_run_script_win(xb, N, maindir, xbeach_exe):
    '''
    Create batch script to run simulations.

    Parameters
    ----------
    xb : LIST or CLASS
        list with simulation paths or XBeachModelSetup class
    N : int
        Number of batch scripts.
    maindir : TYPE
        path where run script is created.
    xbeach_exe : TYPE
        path of XBeach executable.

    Returns
    -------
    None.

    '''
    ## check xb, and create path_sims list
    if isinstance(xb,list):
        if isinstance(xb[0],XBeachModelSetup):
            path_sims = []
            for item in xb:
                path_sims.append(item.model_path)
        elif isinstance(xb[0],str):
            path_sims = [xb]
        else:
            print('Unvalid path')
    else:
        path_sims = [xb.model_path]
    
    ## number of batch scripts
    Nmax = int(np.ceil(len(path_sims)/N))
    
    ## string
    string      = ''
    count       = 0
    run_number  = 0
    for ii, path_sim in enumerate(path_sims):
        string = string + 'cd {} \ncall {}\n'.format(path_sim,xbeach_exe)
        if count==Nmax:
            with open(os.path.join(maindir,'run{}.bat'.format(run_number)), 'w') as f:
                f.write(string)
            count = 0
            run_number = run_number + 1
            string = ''
        count = count +1
    if count<=Nmax:
        print(os.path.join(maindir,'run{}.bat'.format(run_number)))
        with open(os.path.join(maindir,'run{}.bat'.format(run_number)), 'w') as f:
            f.write(string)
    