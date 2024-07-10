# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:10:02 2023

@author: Menno de Ridder, Cas van Bemmelen
main collection for the setup of XBeach models
module contains class to setup 1D and 2D XBeach models based on user input

"""
# general imports
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt

# import rotate grid from geometry
from .general.geometry import rotate_grid
from .general.deg2uv import deg2uv

# Dictionary functions
from xbTools.general.file_utils import (write_2d_arr_2_file, format_subsection_header, 
                                       format_header_line, get_json)


class XBeachModelSetup():
    '''
    XBeach model setup class
    ''' 
    
    def __init__(self,file_name):
        self.file_name      = file_name
        ## by default set wbctype and wavemodel to None
        self.wbctype    = None
        self.wavemodel  = None
        self.zs0type    = None
        
        ## set default values
        self.model_path = None
        self.friction_layer = None
        self.wavefriction_layer = None

        # WaveHello Added these 
        self.wavefriction = None
        self.friction = None
        self.nebed = None
        self.struct = None
        
    def __str__(self):
        """_summary_

        Returns:
            file_name (string): Name of the xBeach model file
        """        
        return self.file_name
    
    def set_params(self,input_par_dict, json_file_name = "par.json"):
        """
        Sets the wavemodel and any of the parameters in the par.json file for the model

        Args:
            input_par_dict (dict): dict of input param flags and the associated values
        """        
        ## set wavemodel. Default is Surfbeat
        if 'wavemodel' not in input_par_dict:
            print('No wavemodel defined. Wavemodel is set to Surfbeat')
            self.wavemodel = 'surfbeat'
        else:
            self.wavemodel = input_par_dict['wavemodel']
        
        ## set wbctype
        if 'wbctype' in input_par_dict:
            self.wbctype = input_par_dict['wbctype'] 
        
        # Get the json file that contains all of the parameters in ther proper subsection so that the self.input_par
        # can be put under the right headers
        folder_path = os.path.dirname(__file__)

        # Add store it in the class
        self.json_param_dict = get_json(folder_path, file_name = json_file_name)

        ## create input dict
        self.input_par = {}

        # Init 'par' key as a nested dict
        self.input_par['par'] = {}

        # Init flag to track if parameters aren't found in the JSON
        not_found_params = False

        ## loop over input parameters 
        for input_par in input_par_dict:
            value_added = False

            ## loop over categories from the par.json file
            for par_category in self.json_param_dict:
                ## if input parameter is in category, add parameter
                if input_par in self.json_param_dict[par_category]:
                    ## create category if not exist
                    if not par_category in self.input_par:
                        self.input_par[par_category] = {}
                        
                    ## add parameter and value                    
                    self.input_par[par_category][input_par] = input_par_dict[input_par]
                    value_added = True

            if not value_added:
                # Add all parameters that don't match a key in the par.json file to this new dict
                self.input_par['par'][input_par] = input_par_dict[input_par]

                not_found_params = True
        
        # Print out the params that weren't found. This will help catch typo errors
        if not_found_params:
            print(f"The following params were not found in the JSON: \n{self.input_par['par']}")

    def set_grid(self,xgr,ygr,zgr, posdwn=1, xori=0, yori=0,alfa=0, thetamin=-90, thetamax = 90, thetanaut = 0, dtheta=10, dtheta_s=10):
        """_summary_

        Args:
            xgr (array): x-grid
            ygr (array): y-grid
            zgr (array): z-grid
            posdwn (int, optional): _description_. Defaults to 1.
            xori (int, optional): _description_. Defaults to 0.
            yori (int, optional): _description_. Defaults to 0.
            alfa (int, optional): _description_. Defaults to 0.
            thetamin (int, optional): _description_. Defaults to -90.
            thetamax (int, optional): _description_. Defaults to 90.
            thetanaut (int, optional): _description_. Defaults to 0.
            dtheta (int, optional): _description_. Defaults to 10.
            dtheta_s (int, optional): _description_. Defaults to 10.
        """        
        ##
        assert xgr.shape==zgr.shape, 'Shape of xgr is not equal to shape of zgr'
        
        ## 1D model
        if ygr is None:
            self.ygr = None
            ## make 2d matrix
            if xgr.ndim==1:
                self.xgr = xgr[np.newaxis, ...] 
                self.zgr = zgr[np.newaxis, ...]  
            else:
                self.xgr = xgr
                self.zgr = zgr
        ## 2D model
        else:
            self.ygr = ygr
            self.xgr = xgr
            self.zgr = zgr
        
        ##
        self.nx = self.xgr.shape[1] - 1
        self.ny = self.xgr.shape[0] - 1
        
        ## 1D
        if ygr is None or ygr.shape[0]==1:
            self.fast1D = True
            self.ny = 0
        else:
            self.fast1D = False 
        
        # TODO: Store the grid data in a dictionary

        ## set values
        self.posdwn = posdwn
        self.xori   = xori
        self.yori   = yori
        self.thetamin   = thetamin
        self.thetamax   = thetamax
        self.thetanaut  = thetanaut
        self.dtheta     = dtheta
        self.vardx  = 1
        self.alfa = alfa
        self.dtheta_s = dtheta_s

    def set_nebed(self, nebed, struct=1):      
        '''
        function to set non erodible bed for the xbeach model

        Parameters
        ----------
        nebed : input of sandy layer thickness
        struct : optional. The default is 1.

        Returns
        -------
        None.

        '''
        self.nebed = nebed
        self.struct = struct

    def set_friction(self, friction, friction_layer = 1):      
        '''
        function to set friction layer for the xbeach model

        Parameters
        ----------
        friction : input of sandy layer thickness
        friction_layer : optional yes/no. The default is 1.

        Returns
        -------
        None.

        '''
        self.friction = friction
        self.friction_layer = friction_layer

        # TODO option to specify what kind of bed friction there is: Manning/Chezy or else

    def set_wavefriction(self, wavefriction, wavefriction_layer = 1):      
        '''
        function to set wave friction layer for the xbeach model

        Parameters
        ----------
        wavefriction : input of sandy layer thickness
        wavefriction_layer : optional yes/no. The default is 1.

        Returns
        -------
        None.

        '''
        self.wavefriction = wavefriction
        self.wavefriction_layer = wavefriction_layer

    @staticmethod
    def _get_wbctype_required_params(wbctype):
        """
        Returns the required inputs for the selcted wbctype.
        NOTE: This is only for the new boundary conditions
        
        Inputs:
            wbctype (string): One of the selected instat boundary condtion types

            Possible inputs for wbctype are: 

            swan: XBeach can read standard SWAN 2D variance density or energy density output files (+.sp2 files) as specified in the SWAN v40.51 manual. 
            
            vardens: 2D spectral information that is not in SWAN format can be provided using a formatted variance density spectrum file
            
            off: 
            
            jonstable:  Each line in the spectrum definition file contains a parametric definition of a spectrum, 
                         like in a regular JONSWAP definition file, plus the duration for which that spectrum is used during the simulation.
            
            reuse: makes XBeach reuse wave time series that were generated during a previous simulation. 
                   This can be a simulation using the same or a different model as long as the computational grids are identical
            
            ts_1: First-order time series of waves (keyword insat = ts_1). 
                  XBeach will calculate the bound long wave based on the theory of [LHS64]).
            
            ts_2: Second-order time series of waves (keyword insat = ts_2). 
                  The bound long wave is specified by the user via a long wave elevation.
            
            
            ts_nonh: Requires a boundary condition file that contains free surface elevations and velocities (both in u and v).
                     Possible input values are z (surf. elevation), t, u, v, w, dU, dV, q (discharge), dq

            parametric: 
        """

        # Using False to note that the required parameters haven't been included here yet
        # TODO: Update the ts_nonh required parameters to make it more general
        wbctype_options = {'swan'        : False,
                          'vardens'      : False,
                           'off'         : False,
                           'jonstable'   : ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc'],
                           'resuse'      : False, 
                           'ts_1'        : False, 
                           'ts_2'        : False, 
                           'ts_nonh'     : ["make_file", "file_name", "dimension", "variable_dict"],
                           'parametric': ['Hm0','Tp','mainang','gammajsp','s','fnyq']
                           }

        try:
            required_params = wbctype_options[wbctype]
        except KeyError as error:
            raise KeyError("Input: {} is not an option for wbctype conditions.\
                           Possible options are: {}".format(wbctype, wbctype_options.keys()))
        
        if required_params is False:
            raise ValueError(" Option is valid but, writing {} not implemented yet".format(wbctype))
        
        return required_params
    
    @staticmethod
    def _get_instat_required_params(instat_type):
        """
        Returns the required inputs for the input boundary condition

        NOTE: This is only for the old boundary condition type - instat

        Inputs:
            instat_type (string): One of the selected instat boundary condtion types

            Possible inputs for instat_type are:
            bichrom: 

            ts_1: First-order time series of waves (keyword insat = ts_1). 
                  XBeach will calculate the bound long wave based on the theory of [LHS64]).
            
            ts_2: Second-order time series of waves (keyword insat = ts_2). 
                  The bound long wave is specified by the user via a long wave elevation.
            
            jons : A JONSWAP wave spectrum is parametrically defined in a file that is referenced using the bcfile keyword. 
                  This file contains a single parameter per line in arbitrary order.
            
            swan: XBeach can read standard SWAN 2D variance density or energy density output files (+.sp2 files) as specified in the SWAN v40.51 manual. 
            
            vardens: 2D spectral information that is not in SWAN format can be provided using a formatted variance density spectrum file
            
            reuse: makes XBeach reuse wave time series that were generated during a previous simulation. 
                   This can be a simulation using the same or a different model as long as the computational grids are identical
            
            ts_nonh: Requires a boundary condition file that contains free surface elevations and velocities (both in u and v).
                     Possible input values are z (surf. elevation), t, u, v, w, dU, dV, q (discharge), dq
            off: 
            stat_table: Only in case of insat = stat_table the time-varying stationary wave boundary conditions are fully described
                       in an external file referenced by the bcfile keyword.
            jonstable:  Each line in the spectrum definition file contains a parametric definition of a spectrum, 
                         like in a regular JONSWAP definition file, plus the duration for which that spectrum is used during the simulation.
        """

        # Using False to note that the required parameters haven't been included here yet
        # TODO: Update the ts_nonh required parameters to make it more general
        instat_type_options = {
                                'bichrom'     : False,
                                'ts_1'        : False, 
                                'ts_2'        : False, 
                                'jons'        : False,
                                'swan'        : False,
                                'vardens'     : False,
                                'resuse'      : False, 
                                'ts_nonh'     : ["make_file", "file_name", "dimension", "variable_dict"],
                                'off'         : False,
                                'stat_table'  : ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc'],
                                'jonstable '  : ['Hm0','Tp','mainang','gammajsp','s','fnyq']
                                }

        try:
            required_params = instat_type_options[instat_type]
        except KeyError as error:
            raise KeyError("Input: {} is not an option for instat_type conditions.\
                           Possible options are: {}".format(instat_type, instat_type_options.keys()))
        
        if required_params is False:
            raise ValueError(" Option is valid but, writing {} not implemented yet".format(instat_type))
        
        return required_params

    def set_waves(self, wbctype, input_struct, instat_bc = False):
        """

        Gets the required inputs for the selected boundary condition and stores the values in a dict for 
        later output. 
   
        Args:
            wbctype (string)   : String input that selects one of the boundary condition options
            input_struct (dict): Dict contains the information
            instat_bc (bool): Flag for if the old boundary conditions should be used (Defaults to False)
        
        Calls:
            _get_wbctype_required_params: To get the required params for the wbctype boundary conditions
            _get_instat_required_params : To get the required params for the instat boundary conditions

        See the called functions for more information on the boundary conditions
        """        
        self.wbctype = wbctype
        
        # Get the thr requried input parameters

        if not instat_bc:
            # Use the wbctype conditions
            required_par = self._get_wbctype_required_params(wbctype)

        elif instat_bc:
            # Use the instat conditions
            required_par = self._get_instat_required_params(wbctype)
        else: 
            TypeError("instat_bc should be a boolean. Type is: {}".format(type(instat_bc)))

        # Store the required parameters for later usage
        self.required_wbc_params = required_par

        # Init dict to hold the waves_boundary conditions
        self.waves_boundary  = {}
                            
        # Loop over the required parameters...
        for item in required_par:
            # Check that the required param is in the input_struct
            if item not in input_struct:
                raise KeyError("Required parameter: {} isn't found in input_struct. Required parameters are: {}".format(item, required_par))
            
            # If param is in input string add it to waves_boundary dict for later output
            self.waves_boundary[item] =  input_struct[item]
            
    def set_vegetation(self):
        """_summary_
        """        
        pass

    def set_wind(self):
        """_summary_
        """        
        pass        

    def set_tide(self, zs0, optional_par_input = {}):
        """
        function to set offshore water level boundary conditions for the xbeach model

        Parameters
        ----------
        zs0 (array): table of water levels in format [t, zs1, (zs2), (zs3, zs4)]
        optional_par_input: dict with optional parameter settings: ['tideloc', 'tidetype', 'paulrevere']. If unspecified, tideloc is inferred from zs0.

        Returns
        -------
        None.
        """     

        if np.atleast_1d(zs0).size > 1:
            self.zs0type = 'list'
        else:
            self.zs0type = 'par'
        ## 

        self.tide_boundary = {}

        if self.zs0type == 'list':
            self.tide_boundary['zs0file'] = 'tide.txt'
            optional_par = ['paulrevere', 'tideloc', 'tidetype']

            for item in optional_par_input:
                if item in optional_par:
                    self.tide_boundary[item] = optional_par_input[item]
                else:
                    assert False, 'invalid tide parameter specified'

            self.tide_boundary['_tidelen'] = zs0.shape[0]

            # if not overruled by a specific par input, infer the tideloc from shape of zs0
            if not 'tideloc' in optional_par_input:
                tideloc_inferred = zs0.shape[1]-1
                if tideloc_inferred in [1, 2, 4]:    
                    self.tide_boundary['tideloc'] = tideloc_inferred
                else:
                    assert False, 'format of zs0 invalid: use either 1, 2 or 4 corners for zs0 specification'
               
        elif self.zs0type == 'par':
            optional_par = []
            self.tide_boundary['zs0'] = zs0

        else:
            assert False, 'Wrong zs0 format'

        self.zs0 =  zs0

        # TODO check that zs0file is at least as long as end time of simulation
    
    @staticmethod
    def write_zsinitfile(init_surface_elev, model_dir, file_name = "zsinitfile.txt", decimal_digits = 5):
        """
        Write the initial water surface elevation file  "zsinitfile"
        """
        
        # Construct the file path
        zs_init_file_path = os.path.join(model_dir, file_name)

        # Convert the data to be a least 2d, so that the writer function can work
        init_surface_elev = np.atleast_2d(init_surface_elev)

        # Write the file to the specified path
        write_2d_arr_2_file(init_surface_elev, zs_init_file_path, decimal_digits)
        

    def load_model_setup(self,path):
        """_summary_

        Args:
            path (_type_): _description_
        """        
        ## todo
        pass    

    def _write_boun_U_file(self, path, boun_U_dict, num_dec_dig):
        """
        Writes the boundary conditon file that is required for ts_nonh models
        
        Inputs:
            self:
            path (string)     : Path to where the file should be written
            boun_U_dict (dict): Dictionary that contains the data needed to make the boun_U.bcf
            file_name (string): Name of the boun_U file. This shouldn't need to change but just in case
        
        Below is an example of what the boun_U_dict should contain
        boun_U_dict = {"make_file": "Boolean",
                       "file_name": "boun_U.bcf",
                       "dimension": "scalar or vector",
                       "variable_dict": {
                           "t" : "time array",
                           "U" : "U-velocity data",
                           "V" : "V-Velocity data",
                           "W" : "W-Velocity data",
                           "dU": "dU data",
                           "dV": "dV data",
                           "q" : "q data",
                           "dq": "dq data"
                           "zs": "surface pertubation data" 
                       }         
            }

        Info:
            make_file (string)  : Used to determine if the boun_U should be created
            file_name (string)  : Name of the boun_U.bcf file
            variable_dict (dict): Dictionary of variable names that xBeach expects and the associated data.
                                  NOTE: Not all of the variables have to be provided the above is just an example 
                                  of what could happen.

        """
        
        variable_dict = boun_U_dict["variable_dict"]

        # Make a list of the allowed variables names
        allowed_variable_names = ['z','zs','t', 'u', 'v', 'w','du', 'dv', 'q', 'dq']

        # Get the variable names
        variable_names = list(boun_U_dict["variable_dict"].keys())

        # Convert them to lower case
        lowercase_variables = [item.lower() for item in variable_names]

        # Check that all the variables are in the name
        mask_list = [item in allowed_variable_names for item in lowercase_variables]

        # Check if there any input variables that aren't allowed
        if not all(mask_list):
            # Apply the mask to the inputs and get the input vars that aren't allowed
            not_allowed_vars = [item for item, mask in zip(variable_names, mask_list) if not mask]

            raise ValueError("Input variable(s): {} are not allowed. \
                             Allowed variables are: {}".format(not_allowed_vars, allowed_variable_names))
        
        # Now that the variables are checked generate the file
        file_path =  os.path.join(path, boun_U_dict["file_name"])

        # Loop over the data and stack it together
        data = np.column_stack([variable_dict[key] for key in variable_dict.keys()])

        string_variables = " ".join(variable_names)
        num_variables = len(variable_names)

        # Write the header and the data to a file
        with open(file_path, 'w') as file:
            file.write("{}\n".format(boun_U_dict["dimension"]))
            file.write("{}\n".format(num_variables))
            file.write("{}\n".format(string_variables))
            np.savetxt(file, data, delimiter=' ', fmt=f'%10.{num_dec_dig}f')

        print("Data written to {}".format(file_path))

    def _write_jonswap_file(self, path, tab_number, file_name = "jonswap.txt"):
        """
        Write the jonswap file

        Inputs:
            self: The current object
            path (string)     : Path that the file should be written into
            tab_number (int)  : Length of the tab used to seperate data in non-table files
            file_name (string): Name of the jonswap file
        """

        # TODO: Add other boundary conditions that use the jonswap.txt file
        write_table_dict = {
            'parametric': False,
            'jonstable' : True
        }

        # Check if the  jonswap should be a table or a single entry of values
        try:
            write_table = write_table_dict[self.wbctype]
        except KeyError as error:
            raise KeyError("Input wbctype: {} doesn't currently support writing a jonswap file. \
                           The conditons that are implemented are: {}".format(self.wbctype, write_table.keys()))

        # Store the name of the boundary condition file
        self.waves_boundary["bcfile"] = file_name

        # Open the jonswap file and write the data...
        with open(os.path.join(path, file_name), 'w') as f:
            # Writing if the data isn't a table
            if not write_table:
                for param in self.required_wbc_params:
                    f.write('{}\t= {}\n'.format(param, self.waves_boundary[param]).expandtabs(tab_number))
            
            # Writing if the data is a table
            elif write_table:
                #TODO: Think of a more elegant way to do this that is more future proof
                # Get number of entries under the first key
                num_vals = len(self.waves_boundary['Hm0'])
                
                # Loop over each row of the table...
                for i in range(num_vals):
                    # Write the param vals into each column
                    for param in self.required_wbc_params:
                        f.write('{} '.format(self.waves_boundary[param][i]))
                    
                    # Set up the next row
                    f.write('\n')
    
    def _write_wbc_file(self, path, num_dec_dig):
        """
        Wrapper for the functions that write the files that are needed for the wave boundary conditons
        """ 

        tab_number = 20
        
        # Modes that require the boun_U_dict file
        if self.wbctype == "ts_nonh":
            boun_U_dict = self.waves_boundary
            
            # Check if the file should be created...
            # The first index should be "make_file"
            # TODO: Update the required parameters to be a 
            if boun_U_dict["make_file"]:
                # If True make the file
                self._write_boun_U_file(path, boun_U_dict, num_dec_dig)

        # wbctypes that require the jonswap.txt file
        elif self.wbctype == "parametric" or self.wbctype == "jonstable":
            self._write_jonswap_file(path, tab_number, file_name = "jonswap.txt")

        elif self.wbctype == "off":
            # Don't need to do anything here
            pass
        else:
            raise Warning("Writing the files for wbctype: {} is not implemented".format(self.wbctype))
    
    @staticmethod
    def _write_params_metadata(file, current_date, user):
        """
        Write the Meta data of the params file
        """
        # Def a function to format the header lines 

        
        function_name = "_write_params_metadata"

        header = [
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",
        format_header_line("XBeach parameter settings input file"),
        format_header_line(""),
        format_header_line(f"date:     {current_date}"),
        format_header_line(f"Params created by {user}"),
        format_header_line(f"function: {function_name}"),
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n",
        ]
        
        # Make the list into a string and add the newline values
        header = "\n".join(header)

        # Write the metadata
        file.write(header)

    @staticmethod
    def _write_params_general_data(file, wavemodel, wbctype, tab_number):
        """
        Write the general data into the params file

        TODO: WaveHello - In the process of deprecating this function
        """

        ## general
        file.write('\n')
        if wavemodel is not None:
            file.write('wavemodel\t= {}\n'.format(wavemodel).expandtabs(tab_number))

        if wbctype is not None:
            file.write('wbctype\t= {}\n'.format(wbctype).expandtabs(tab_number))
        file.write('\n')

    def _write_params_grid_data(self, file, tab_number):
        """
        TODO: Check if this needs to be a internal module. 
        For the time being this is going to be the only one that stays inside of the class"""

        # TODO: Once the grid data is stored in a dictionary update this to follow the format of the other parameters

        ## grid
        file.write(f"{format_subsection_header("Grid parameters")}\n" )
        file.write('\n')
        file.write('vardx\t= {}\n'.format(self.vardx).expandtabs(tab_number))
        file.write('posdwn\t= {}\n'.format(self.posdwn).expandtabs(tab_number))
        file.write('nx\t= {}\n'.format(self.nx).expandtabs(tab_number))
        file.write('ny\t= {}\n'.format(self.ny).expandtabs(tab_number))
        file.write('xori\t= {}\n'.format(self.xori).expandtabs(tab_number))
        file.write('yori\t= {}\n'.format(self.yori).expandtabs(tab_number))
        file.write('alfa\t= {}\n'.format(self.alfa).expandtabs(tab_number)) 
        file.write('xfile\t= x.grd\n'.expandtabs(tab_number))

        # Check that the model has a y-grid (ie. that's its 2d)
        if not self.ygr is None:
            file.write('yfile\t= y.grd\n'.expandtabs(tab_number))

        file.write('depfile\t= bed.dep\n'.expandtabs(tab_number))
        file.write('thetamin\t= {}\n'.format(self.thetamin).expandtabs(tab_number))
        file.write('thetamax\t= {}\n'.format(self.thetamax).expandtabs(tab_number))
        file.write('thetanaut\t= {}\n'.format(self.thetanaut).expandtabs(tab_number))
        file.write('dtheta\t= {}\n'.format(self.dtheta).expandtabs(tab_number))
        file.write('dtheta_s\t= {}\n'.format(self.dtheta).expandtabs(tab_number))
        file.write('\n')

    def _write_params_tide_data(self, file, tab_number):
        """
        Write the tide data to the file
        """
        file.write(f"{format_subsection_header("Tide Boundary Conditions")}\n" )
        file.write('\n')     

        if self.zs0type == 'par':
            file.write('zs0\t= {}\n'.format())    

        elif self.zs0type == 'list':
            for item in self.tide_boundary:
                if item[0] != '_':
                    file.write('{}\t= {}\n'.format(item, self.tide_boundary[item]).expandtabs(tab_number))
       
        file.write('\n')

    def _write_params_input_vars(self, file, tab_number):
        """
        Writes the input parameters (self.input_par) to the params.txt file in the specified file

        Inputs:
            self: Current instance of class
            file: File that the param data should be written to
            tab_number: Number of tabs that should follow a parameter in the file

        """

        # Create a new ordered dictionary based on the order of json file
        # This way the order follows the JSON file
        # The JSON file can be put in the same order of the xBlog file making comparing the sections easier
        ordered_param_dict = {k: self.input_par[k] for k in self.json_param_dict if k in self.input_par}

        if "wbctype" not in ordered_param_dict["Wave boundary condition parameters"].keys():
            print(("\nWARNING: Include the wbctype in the set_params function. "
                   "In the process of moving away from printing the "
                   "wbctype at the top of the params.txt file"
                ))
        
        # Loop over the parameter categories in the input dict...
        for par_category in ordered_param_dict:
            ## skip category starting with _. The name of the category is set in the JSON file
            if par_category[0]=='_':
                continue
            
            ## write meta
            # Write the header of the category the parameter is in
            file.write(f"{format_subsection_header(par_category)}\n" )
            file.write('\n')

            # Go through all of the inputs in the category write those values
            for par in ordered_param_dict[par_category]:
                file.write('{}\t= {}\n'.format(par,ordered_param_dict[par_category][par]).expandtabs(tab_number))

            file.write('\n')
    
    def _write_params_output_vars(self, file, tab_number):
        """
        Write the variables that should be output by the xBeach model into the params file"""
        
        ## write output variables
        if '_Output' in self.input_par:
            file.write(f"{format_subsection_header("Output variables")}\n" )
            file.write('\n')
            for par in self.input_par['_Output']:
                dummy = self.input_par['_Output'][par]
                file.write('{}\t= {}\n'.format(par,len(dummy)).expandtabs(tab_number))

                if not isinstance(dummy, list):
                    raise TypeError("Expected a list for {}".format(par))

                for item in dummy:
                    file.write('{}\n'.format(item))
                file.write('\n')

    def _write_params_file(self, path_params, current_date, user, tab_number):

        """
        Write the params file. Acts as a wrapper for other functions that actually do the writing
        """

        ## Open and create the file
        with open(path_params,'w') as f:
            
            # TODO: Move these functions into a util file to make this class smaller

            # Write the metadata at the top of the file
            self._write_params_metadata(f, current_date, user)

            # Write the general data
            # 
            # NOTE: WaveHello - Commented out because the data is being written under the correct subheader now
            # self._write_params_general_data(f, self.wavemodel, self.wbctype, tab_number)
            
            self._write_params_grid_data(f, tab_number)

            ## tide 
            if self.zs0type is not None:
                self._write_params_tide_data(f, tab_number)

            ## write input vars
            self._write_params_input_vars(f, tab_number)

            # Write the output vars
            self._write_params_output_vars(f, tab_number)

    def write_model(self, path, figure=False, num_dec_dig = 5):
        """
        Wrapper for functions that write the input files for the model

        Args:
            path (string): Path to where the xBeach data is written to
            figure (bool, optional): _description_. Defaults to False.
        """        
        self.model_path = path
        path_params = os.path.join(path,'params.txt')
        
        # Raise an error if the folder isn't found
        if not os.path.exists(path):
            raise FileExistsError('{} does not exist'.format(path))
        
        current_date    = datetime.today().strftime('%Y-%m-%d %HH:%mm')
        user            =  os.path.basename(os.path.expanduser('~'))
        
        tab_number = 20
        
        # Write the wave boundary condition files
        self._write_wbc_file(path, num_dec_dig)
        
        # Write the params file
        self._write_params_file(path_params, current_date, user, tab_number)
    
        # Do a loop for the other files
        file_arr_dict = {"x.grd": self.xgr,
                         "y.grd": self.ygr, 
                         "bed.dep": self.zgr,
                         "ne_bed.dep": self.nebed,
                         "friction.dep": self.friction,
                         "wavefriction.dep": self.wavefriction
                         } 
        
        # Loop over the dict...
        for key, value in file_arr_dict.items():
            # Check that the value isn't none
            if value is not None:
                # Treat the arr as at least 2d, this allows the 1d models to be written with the 2d func
                value = np.atleast_2d(value)
                
                # Write the file
                dummy_path = os.path.join(path, key)
                write_2d_arr_2_file(value, dummy_path)        

        ## write tide boundary condition
        if self.zs0type == 'list':
            with open(os.path.join(path,'tide.txt'), 'w') as f:
                for ir in range(self.tide_boundary['_tidelen']):
                    for ic in range(self.tide_boundary['tideloc'] + 1):
                        f.write('{} '.format(self.zs0[ir, ic]))
                    f.write('\n')   

        ## write figures
        if figure:
            ## plot and write domain
            self._plotdomain(path)
            ## plot and write wave boundary
            if self.wbctype=='jonstable' or self.wbctype=='parametric':
                self._plot_boundary(path)

    def _plot_boundary(self,save_path=None):
        '''
        Plot boundary conditions

        Parameters
        ----------
        save_path : TYPE, optional
            Path were figure is saved. The default is None.

        Returns
        -------
        None.

        '''
        if self.wbctype=='jonstable':
            plt.figure()
            plt.subplot(3,1,1)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['Hm0'],'-o')
            plt.ylabel('$H_{m0}$')
            plt.subplot(3,1,2)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['Tp'],'-o')
            plt.ylabel('$T_{p}$')
            plt.subplot(3,1,3)
            plt.plot(np.cumsum(self.waves_boundary['duration']), self.waves_boundary['mainang'],'-o')
            plt.ylabel('$D$')
            plt.xlabel('Time')
            if save_path!=None:
                plt.savefig(os.path.join(save_path,'jonstable.png'))
        elif self.wbctype=='parametric':
            print('wbctype=parametric cannot be plotted')
        else:
            print('Not possible to plot wave boundary')
            
    def _make_theta_vectors(self):
        """function to generate thetavectors for plotting arrows in the domain plot

        Raises:
            UserWarning: _description_

        Returns:
            _type_: _description_
        """        
        # make theta vectors
        if self.thetanaut == 1: #use nautical conversion
            thetamin_uv = deg2uv(self.thetamin)
            thetamax_uv = deg2uv(self.thetamax)
        
        elif self.thetanaut == 0: #use cartesian conversion
            thetamin_naut = self.thetamin+90
            thetamax_naut = self.thetamax+90
            thetamin_uv = deg2uv(thetamin_naut)
            thetamax_uv = deg2uv(thetamax_naut)
        
        else: 
            raise UserWarning(f"you cannot use thetanaut = {self.thetanaut}, please keep it 0/1")
        
        return thetamin_uv, thetamax_uv

    def _plotdomain(self,save_path=None):
        '''
        Plot the domain. 

        Parameters
        ----------
        save_path : string, optional
            Path were figure is saved. The default is None.

        Returns
        -------
        None.

        '''
        fig1 = plt.figure(figsize=(8,8))
        thetamin_uv, thetamax_uv = self._make_theta_vectors()
        if self.fast1D==True:
            plt.subplot(2,1,1)
            plt.plot(np.squeeze(self.xgr),np.squeeze(self.zgr)*-self.posdwn)
            plt.ylabel('z')
            plt.subplot(2,1,2)
            plt.plot(np.squeeze(self.xgr)[1:],np.diff(np.squeeze(self.xgr)))
            plt.xlabel('x')
            plt.ylabel('dx')
        else:
            plt.subplot(2,2,1)
            plt.pcolor(self.xgr,self.ygr,self.zgr*-self.posdwn)
            plt.ylabel('y')
            plt.colorbar()
            plt.axis('equal')
            plt.title('Local coordinates - input')

            # add theta grid vectors
            ax = plt.subplot(2,2,2)

            ax.arrow(-thetamin_uv[0], -thetamin_uv[1], thetamin_uv[0], thetamin_uv[1], length_includes_head=True, head_width=0.08, head_length=0.2)
            ax.arrow(-thetamax_uv[0], -thetamax_uv[1], thetamax_uv[0], thetamax_uv[1], length_includes_head=True, head_width=0.08, head_length=0.2)
            
            ax.text(-thetamin_uv[0], -thetamin_uv[1], '$\Theta_{min}$ = '+f'{self.thetamin:.2f}', ha = 'left', va = 'center')
            ax.text(-thetamax_uv[0], -thetamax_uv[1], '$\Theta_{max}$ = '+f'{self.thetamax:.2f}', ha = 'left', va = 'center')
            
            ax.set_aspect('equal')
            ax.set_title(f'thetanaut = {self.thetanaut}')
            ax.axis('off')

            plt.subplot(2,2,3)
            [X_world,Y_world] = rotate_grid(self.xgr,self.ygr,np.deg2rad(self.alfa))
            plt.pcolor(X_world+self.xori,Y_world+self.yori,self.zgr*-self.posdwn)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
            plt.colorbar()
            plt.title('World coordinates - output')

        plt.suptitle(self.file_name)

        if self.struct == 1:
            if not self.fast1D == True:
                fig2 = plt.figure()

                plt.pcolor(self.xgr, self.ygr, self.nebed)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.colorbar()
                plt.title("ne_bed.dep (positive)")
                plt.axis("scaled")
                plt.grid("on")
                
        if self.friction_layer == 1:
            if not self.fast1D == True:
                fig3 = plt.figure()

                plt.pcolor(self.xgr, self.ygr, self.friction)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.colorbar()
                plt.title(
                    "friction.dep (positive)"
                )  # +self.input_par['Flow parameters']['bedfriction'])
                plt.axis("scaled")
                plt.grid("on")
                
        if self.wavefriction_layer == 1:
            if not self.fast1D == True:
                fig4 = plt.figure()

                plt.pcolor(self.xgr, self.ygr, self.wavefriction)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.colorbar()
                plt.title("wavefriction.dep (positive) - fw")
                plt.axis("scaled")
                plt.grid("on")

        if save_path!=None:
            fig1.savefig(os.path.join(save_path, "domain.png"), dpi=250)
            if self.struct == 1 and not self.fast1D == True:
                fig2.savefig(os.path.join(save_path, "ne_bed.png"), dpi=250)
            if self.friction_layer == 1 and not self.fast1D == True:
                fig3.savefig(os.path.join(save_path, "friction.png"), dpi=250)
            if self.wavefriction_layer == 1 and not self.fast1D == True:
                fig4.savefig(os.path.join(save_path, "wavefriction.png"), dpi=250)