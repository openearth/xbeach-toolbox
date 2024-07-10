# -*- coding: utf-8 -*-
"""
Created on Wed May 5 10:04:00 2023

collection containing executing run functionalities
"""
# import python modules
import os
import numpy as np
import subprocess

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
    maindir : string
        path where run script is created.
    xbeach_exe : string
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
            path_sims = xb
        else:
            print('Unvalid path')
    else:
        path_sims = [xb.model_path]
    
    ## number of batch scripts
    # TODO: This will cause an error if a string is passed into this function. This would return the length of the string divided by the N
    # TODO: Don't think this makes any sense
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

def generate_batch_script(model_folder, exe_path, batch_file_name = "run_model.bat", include_cd = False, batch_file_folder = None):
    """
    Generates an batch script inside of the model directory, linking to the exe_path
    
    Parameters
    ----------
    model_folder: string
        Directory to the model that the batch script should be created for
    
    exe_path: string
        Directory to the executable that should be linked to in the batch file
    
    batch_file_name: string (Optional)
        Name of the batch script that should be generated. Defaults to "run_model.bat"
    
    include_cd: Boolean (Optional)
        Controls if a "cd {model directory}" statement is included in the batch file.

    batch_file_folder: string (Optional)
        Input a directory that should hold the batch file. This overrides including the file inside of the model folder
        NOTE: This causes include_cd to be True
    Returns
    -------
    None.
    
    """

    if batch_file_folder is None:
        # Make the batch file path
        batch_script_path = os.path.join(model_folder, batch_file_name)
    
    else:
        # Make the include_cd statement true, because for the batch file to be saved in a different directory than the model it will have to cd into
        # that directory
        include_cd = True
        batch_script_path = os.path.join(batch_file_folder)

    # init empty string to possibly hold the change directory statement
    cd_str = ""

    if include_cd:
        cd_str = "cd \"{}\"".format(model_folder)

    # Make the string that will be written to the file
    batch_str = cd_str + "call \"{}\"".format(exe_path)

    # Create the file and open it in write mode
    with open(batch_script_path, "w") as f:
        # write the string to the batch file
        f.write(batch_str)

def run_batch_script(batch_script_path, flag_print_Blog = False):
    # TODO: Need to test this some more. There are some troubles when the permissions 
    """
    Run a batch script given a path
    
    Parameters
    ----------
    batch_script_path : string
        The path to a batch script that the user wants to run. This will cause the script to be run inside of the python script

    Returns
    -------
    None.

    TODO: There's some trouble that when this is imported it tends to not work. When it is directly in the .ipynb file it works fine.
    TODO: Note sure what's going on with it
    '''
    """

    # Set the working directory to where the batch file is located
    working_directory = os.path.dirname(batch_script_path)

    try:
        # Execute the batch file, capturing both stdout and stderr
        result = subprocess.run(batch_script_path, check=True, shell=True, cwd=working_directory, capture_output=True, text=True)
        
        # Print success message and output
        print(f"Batch file '{batch_script_path}' executed successfully.")

        if flag_print_Blog:
            print("Output:")
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        # Print error message and captured stderr
        print(f"An error occurred while executing the batch file: {e}")
        print("Error output:")
        print(e.stderr)

if __name__ == "__main__":
    pass
