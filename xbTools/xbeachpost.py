# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:10:02 2023

@author: Marlies van der Lugt
main collection for the postprocessing of XBeach models
module contains class for analysis of 2D XBeach models 
note that output and input must be available in directory
"""
import numpy as np
import pandas as pd
import os
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FuncFormatter, IndexLocator
from scipy.interpolate import interp2d
import xbTools as xb
import warnings
import re
# toolbox specific import
from .general.geometry import rotate_grid, path_distance

class XBeachModelAnalysis():
    '''
    XBeach model analysis class
    '''

    def __init__(self, fname, model_path):
        """_summary_

        Args:
            fname (string): File name - Doesn't seem to be used for anything
            model_path (string): Path to folder containing file
        """        
        self.fname = fname
        self.model_path = model_path


        # self.params = None
        self.grd = {}
        self.waves_boundary = {}
        self.dat = {}
        self.tide = {}
        self.wind = {}
        self.var = {}
        self.units = {}
        self.long_name = {}
        self.save_fig = False
        self.plot_localcoords = False
        self.plot_km_coords = False
        self.AOI = []
        self.globalstarttime = None

        self.vector_vars = ['u', 'v', 'ue', 've', 'Subg', 'Svbg', 'Susg', 'Svsg',
                          'u_mean', 'v_mean', 'ue_mean', 've_mean', 'Subg_mean', 'Svbg_mean', 'Susg_mean', 'Svsg_mean',
                          'u_max', 'v_max', 'ue_max', 've_max', 'Subg_max', 'Svbg_max', 'Susg_max', 'Svsg_max',
                          'u_min', 'v_min', 'ue_min', 've_min', 'Subg_min', 'Svbg_min', 'Susg_min', 'Svsg_min',
                          ]

        self.get_params()
        self.get_metadata()
        self.load_grid()
        self._cross_offset = 0
        self.load_output_coordinates()

    def get_metadata(self):
        '''
        function to print metadata from XBlog.txt
        '''
        with open(self.model_path+'\\XBlog.txt') as file:
        # setup meta dictionary
            meta_dictionary = {}
            # loop over all lines in file
            for line in file.readlines():
                if '=' in line: # check if = is in line
                    key, value = line.split('=') # plit on '='
                    meta_dictionary[key.replace(" ", "")] = value[:-1] 
                    # .replace(" ", "") to remove trailing and leading spaces in key
                    # [:-1] used to remove trailing /n
                if 'Finished reading input parameters' in line: # notice finished readline line, then break
                    break
        self.metadata = meta_dictionary # store to class

        self.get_mpi_boundaries()

    def get_mpi_boundaries(self):
        '''
        function to get mpi boundaries from XBlog.txt
        '''
        # first check if it is run with mpi
        if np.sum([1 if 'computational domains on processors' in x else 0 for x in self.metadata])==0:
            self.mpi_ix = []
            self.mpi_iy = []           
            return
        else:
            icdps = [True if 'computational domains on processors' in x else False for x in self.metadata]
            icdp = np.argwhere(icdps)[0][0]
            iadd = 0
            while not ('----') in line:       
                iadd+=1        
                data.append([float(x) for x in self.metadata[int(iadd+icdp)].split()])
                line = self.metadata[int(iadd+icdp+1)]
            if data!=[]:
                # do something with the data
                data = np.array(data)
                ixstop = data[:, 2]
                iystop = data[:, 4]
                self.mpi_ix = np.unique(ixstop)
                self.mpi_iy = np.unique(iystop)
            return


    def __repr__(self):
        return self.fname

    def set_save_fig(self, yesno):
        """_summary_

        Args:
            yesno (_type_): _description_
        """        
        assert type(yesno) is bool, 'input type must be bool'
        self.save_fig = yesno
        plotdir = os.path.join(self.model_path, 'fig')
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)

    def set_plot_localcoords(self, yesno):
        """_summary_

        Args:
            yesno (_type_): _description_
        """        
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_localcoords = yesno
        self.plot_km_coords = False
        self.load_output_coordinates()

    def set_plot_km_coords(self, yesno):
        """_summary_

        Args:
            yesno (_type_): _description_
        """        
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_km_coords = yesno
        self.load_output_coordinates()

    def set_aoi(self, AOI):
        """_summary_

        Args:
            AOI (_type_): _description_
        """        
        self.var = {}  # drop variables from memory because they need to be reloaded with appropriate AOI
        self.AOI = AOI

    def set_globalstarttime(self, tstart):
        """_summary_

        Args:
            tstart (_type_): _description_
        """        
        assert type(tstart) is str, 'tstart must be given as a string of the format 2021-10-11T13:00:00'
        self.globalstarttime = np.datetime64(tstart)

    def get_params(self):
        """
        reads parameters from params
        """

        #todo: better to read these from the xbeach.log instead of the params.txt

        assert os.path.exists(os.path.join(self.model_path, 'params.txt')), 'missing params.txt'

        f = open(os.path.join(self.model_path, 'params.txt'), 'r')
        dat = f.read().split('\n')

        # remove empty elements
        dat = [x for x in dat if not(len(x)==0)]

        # remove elements that are commented out with non alfanumeric character
        dat = [x for x in dat if (x.replace(' ', '')[0].isalnum())]    

        # read params from params file
        params = {}
        # TODO: This get's caught up when there's an equal sign in a comment. Should check if there's a 
                # comment sign on the line before trying to unpack values
        for d in dat:

            if '=' in d:
                x1, x2 = d.split('=')
                if x2.strip().isnumeric():
                    params[x1.strip()] = float(x2.strip())
                else:
                    params[x1.strip()] = x2.strip()

        # global variables
        ixlist = [i for i, var in enumerate(dat) if 'nglobalvar' in var]
        if len(ixlist) > 0:
            i0 = ixlist[0]
            params['globalvar'] = dat[i0+1:i0+int(params['nglobalvar']+1)]
        else:
            params['nglobalvar'] = 0
            params['globalvar'] = []

        # mean variables
        ixlist = [i for i, var in enumerate(dat) if 'nmeanvar' in var]
        if len(ixlist) > 0:
            i0 = ixlist[0]
            params['meanvar'] = dat[i0+1:i0+int(params['nmeanvar']+1)]
        else:
            params['nmeanvar'] = 0
            params['meanvar'] = []

        # point variables
        ixlist = [i for i, var in enumerate(dat) if 'npointvar' in var]
        if len(ixlist) > 0:
            i0 = ixlist[0]
            params['pointvar'] = dat[i0+1:i0+int(params['npointvar']+1)]
        else:
            params['npointvar'] = 0
            params['pointvar'] = []

        # output points
        if params['npointvar'] > 0:
            i0 = [i for i, var in enumerate(dat) if 'npoints' in var][0]
            points = dat[i0 + 1:i0 + int(params['npoints'] + 1)]

            # get points 
            p_spl = [ re.split('\t| ',t) for t in points]
            # remove empty entries
            p_spl = [[p for p in p_spli if p.split()] for p_spli in p_spl]
            # xy
            x = [float(p[0] ) for p in p_spl]
            y = [float(p[1]) for p in p_spl]
            # if names are there list of names else index
            name = [p[2].strip() if len(p)>2 else str(ip) for ip, p in enumerate(p_spl)]

            params['points'] = dict(zip(name, zip(x, y)))
        else:
            params['points'] = {}

        self.params = params

    def load_grid(self):
        """_summary_
        """        
        # only works currently for xbeach type grids (xfile, yfile)

        # load obligatory files
        self.grd['x'] = np.loadtxt(os.path.join(self.model_path, self.params['xfile']))
        if self.grd['x'].ndim==1:
            self.grd['x']= self.grd['x'][np.newaxis, ...] 

        if (self.grd['x'].shape != (self.params['ny']+1, self.params['nx']+1)):
            print('warning: x grid not of size specified in params.txt')

        self.grd['z'] = np.loadtxt(os.path.join(self.model_path, self.params['depfile']))
        if self.grd['z'].ndim==1:
            self.grd['z']= self.grd['z'][np.newaxis, ...] 
        if (self.grd['z'].shape != (self.params['ny'] + 1, self.params['nx'] + 1)):
            print('warning: z grid not of size specified in params.txt')

        ## set posdwn
        if 'posdwn' in self.params:
            if int(self.params['posdwn']) == 1:
                self.grd['z'] = -1*self.grd['z']

        # read y
        if self.params['ny'] > 0:
            self.grd['y'] = np.loadtxt(os.path.join(self.model_path, self.params['yfile']))
            if (self.grd['y'].shape != (self.params['ny'] + 1, self.params['nx'] + 1)):
                print('warning: y grid not of size specified in params.txt')
        # struct
        if ('struct' in self.params) == 1:
            if self.params['struct']==1:
                self.grd['ne'] = np.loadtxt(os.path.join(self.model_path, self.params['ne_layer']))
                assert np.atleast_2d(self.grd['ne']).shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'ne grid not of correct size'

    def get_waves(self):
        """_summary_
        """
        # waves boundary
        if self.params['wbctype'] == 'jonstable':
            dat = np.loadtxt(os.path.join(self.model_path, self.params['bcfile']))
            assert dat.shape[1] == 7, \
                'columns of jonstable should exactly be: Hm0, Tp, mainang, gammajsp, s, duration, dtbc'
            self.waves_boundary['Hm0'] = dat[:, 0]
            self.waves_boundary['Tp'] = dat[:, 1]
            self.waves_boundary['mainang'] = dat[:, 2]
            self.waves_boundary['gammajsp'] = dat[:, 3]
            self.waves_boundary['s'] = dat[:, 4]

            if self.globalstarttime is None:
                self.waves_boundary['time'] = np.cumsum(dat[:, 5])
            else:
                loctime = np.cumsum(dat[:, 5])-dat[0, 5]
                globtime = np.array([np.timedelta64(int(x), 's') for x in loctime]) + self.globalstarttime
                self.waves_boundary['time'] = globtime

        elif self.params['wbctype'] == 'jons':
            print('not yet written')
            pass
        elif self.params['wbctype'] == 'params':
            self.waves_boundary['Hm0'] = np.sqrt(2)*float(self.params['Hrms'])
            self.waves_boundary['Tp'] = float(self.params['Trep'])
            self.waves_boundary['mainang'] = float( self.params['dir0'])
            self.waves_boundary['s'] = float(self.params['m']/2)

        elif self.params['wbctype'] == 'ts_nonh':
            dat = pd.read_csv(os.path.join(self.model_path, 'boun_U.bcf'), skiprows=2, header=0, sep=' ')
            for key in dat.keys():
                if key=='t':
                    if self.globalstarttime is None:
                        self.waves_boundary['time'] = dat['t'].values
                    else:
                        globtime = np.array([np.timedelta64(int(1e6*x), 'us') for x in dat['t'].values]) + self.globalstarttime
                        self.waves_boundary['time'] = globtime
                else:
                    self.waves_boundary[key] = dat[key].values
        else:
            print('not possible')
            pass

        return self.waves_boundary

    def get_vegetation(self):
        """_summary_
        """        
        pass

    def get_tide(self):
        """_summary_
        """
        if 'zs0file' in self.params:
            dat = np.loadtxt(os.path.join(self.model_path, self.params['zs0file']))

            if int(self.params['tideloc']) == 1:
                # assert dat.shape[1] == 2, 'tideloc=1, expected 2 cols'
                pass
            if int(self.params['tideloc']) == 2:
                assert dat.shape[1] == 3, 'tideloc=2, expected 3 cols'
            if int(self.params['tideloc']) == 4:
                assert dat.shape[1] == 5, 'tideloc=1, expected 5 cols'

            self.tide['time'] = dat[:, 0]
            self.tide['zs0'] = dat[:, 1:]
            self.tide['tideloc'] = self.params['tideloc']
            if 'paulrevere' in self.params:
                self.tide['paulrevere'] = self.params['paulrevere']

            if self.globalstarttime is None:
                self.tide['time'] = dat[:, 0]
            else:
                loctime = dat[:, 0]
                globtime = np.array([np.timedelta64(int(x), 's') for x in loctime]) + self.globalstarttime
                self.tide['time'] = globtime

            return self.tide
        
        else:
            print('no timevarying boundary file imposed')
            return self.params['zs0']
        

    def get_wind(self):
        """_summary_
        """        
        if 'wind' in self.params:   
            dat = np.loadtxt(os.path.join(self.model_path, self.params['windfile']))

            self.wind['time'] = dat[:, 0]
            self.wind['u10'] = dat[:, 1]
            self.wind['u10dir'] = dat[:, 2]
            
    def load_model_setup(self):
        """_summary_
        """        
        #self.load_wind()
        self.load_grid()
        self.get_waves()
        self.get_vegetation()
        self.get_tide()

    def load_output_coordinates(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))

        # grid
        x = ds.variables['globalx'][:]
        y = ds.variables['globaly'][:]
        if len(self.AOI) > 0:
            assert len(self.AOI) == 4, 'AOI should be specified as [i_y0, i_yend, i_x0, i_xend]'

            x = x[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
            y = y[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]

        if self.plot_km_coords:
            self.var['globalx'] = x / 1e3
            self.var['globaly'] = y / 1e3
        else:
            self.var['globalx'] = x
            self.var['globaly'] = y

        self.var['gridang'] = np.arctan2(y[0, -1] - y[0, 0], x[0, -1] - x[0, 0]) # np.pi/2+np.arctan2(y[0, -1] - y[0, 0], x[0, -1] - x[0, 0]) # np.arctan2(y[0, 0] - y[-1, 0], x[0, 0] - x[-1, 0])

        # global variable time
        if self.globalstarttime is None:
            self.var['globaltime'] = ds.variables['globaltime'][:]
        else:
            self.var['globaltime'] = np.array([np.timedelta64(int(1e3*x), 'ms') for x in ds.variables['globaltime'][:].data]) \
                                      + self.globalstarttime
        # mean variable time
        if self.params['nmeanvar'] > 0:
            if self.globalstarttime is None:
                self.var['meantime'] = ds.variables['meantime'][:]
            else:
                data = ds.variables['meantime'][:]
                data = data[data.mask == False]
                if len(data)>0:
                    self.var['meantime'] = np.array(
                        [np.timedelta64(int(1e3*x), 'ms') for x in data.data.flatten()]) \
                                            + self.globalstarttime
                else:
                    self.var['meantime'] = []
        # point variable time
        if self.params['npointvar'] > 0:
            if self.globalstarttime is None:
                self.var['pointtime'] = ds.variables['pointtime'][:]
            else:
                data = ds.variables['pointtime'][:]
                data = data[data.mask == False]
                self.var['pointtime'] = np.array(
                    [np.timedelta64(int(1e6*x), 'us') for x in data.data.flatten()]) \
                                         + self.globalstarttime

            station_list = []
            dat = ds.variables['station_id'][:]
            nb, nstat = dat.shape
            for ib in range(nb) :
                sn = ''
                for istat in range(nstat):
                    sn += dat.data[ib, istat].decode('UTF-8')
                station_list.append(sn)
            station_list = [x.strip() for x in station_list]
            self.var['station_id'] = station_list

            self.var['station_x'] = ds.variables['pointx'][:]
            self.var['station_y'] = ds.variables['pointy'][:]

        cross = path_distance(x[0, :], y[0, :])
        self.var['cross'] = cross + self._cross_offset

        self.var['along'] = path_distance(x[:, 0], y[:, 0])

        # self.var['localy'], self.var['localx'] = np.meshgrid(self.var['cross'], self.var['along'])
        # self.var['localx'] = np.flipud(self.var['localx'])  # to plot
        lx, ly = rotate_grid(self.var['globalx'], self.var['globaly'], self.var['gridang'])
        
        self.var['localx'] = lx - lx[0,0]
        self.var['localy'] = ly - ly[0,0]        

        if self.params['npointvar'] > 0: 
            # point output coordinates in local coordinates
            self.var['station_x_local'] = []
            self.var['station_y_local'] = []
            for sx, sy in zip(self.var['station_x'], self.var['station_y']):
                iy, ix = np.unravel_index(((self.var['globalx'][:] - sx) ** 2 + \
                                        (self.var['globaly'][:] - sy) ** 2).argmin(),
                                        [int(self.params['ny'] + 1), int(self.params['nx'] + 1)])
                self.var['station_x_local'].append(self.var['localx'][iy, ix])
                self.var['station_y_local'].append(self.var['localy'][iy, ix])


        return

    def add_cross_offset(self, offset):
        """
        adds an offset value to the computed cross-shore local coordinate system (such that it does not start at zero at offshore boundary)
        """
        self._cross_offset = offset
        self.load_output_coordinates()

    def load_modeloutput(self, var):
        """_summary_

        Args:
            var (_type_): _description_
        """        
        if var in self.var:
            print("Variable already loaded")
            return

        if '_mean' in var:
            assert sum([var[:-5] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'.format(var)
        elif any(substr in var for substr in ['min', '_max', '_var']):
            assert sum([var[:-4] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'.format(var)
        elif 'point_' in var:
            assert sum([var[6:] in x for x in self.params['pointvar']]) > 0, '{} not in xb output'.format(var)
        else:

        # Check if the variable is one of the outputted global variables from the model
        if var not in self.params['globalvar']:
            raise IndexError("Variable is not an available globalvar output." + 
                             " The available outputs are {}".format(self.params["globalvar"]))

        # if not yet present, load coordinates
        if self.var == {}:
            print('loading model output coordinates from file')
            self.load_output_coordinates()

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))
        print('loading variable {} from file'.format(var))
        dat = ds.variables[var][:]

        #mean and point output might not be available if the eor is not reached. Therefore cut
        if 'point_' in var and len(np.atleast_1d(dat.mask)) > 1:
            if len(dat.shape)==2:
                dat = dat[~dat.mask[:, 0], :]
            elif len(dat.shape)==3:
                dat = dat[~dat.mask[:, 0, 0], 0, :]
        elif dat.ndim == 3 and any(substr in var for substr in ['_mean', '_min', '_max', '_var']) and len(np.atleast_1d(dat.mask)) > 1:
            dat = dat[~dat.mask[:, 0, 0], :, :]
        elif dat.ndim == 4 and any(substr in var for substr in ['_mean', '_min', '_max', '_var']) and len(np.atleast_1d(dat.mask)) > 1:
            dat = dat[~dat.mask[:, 0, 0, 0], :, :, :]
        elif dat.ndim == 2 and '_mean' in var and len(np.atleast_1d(dat.mask)) > 1:
            # 1D case
            pass

        if not ('point_' in var):
            if len(self.AOI) > 0:
                print('slicing map output to AOI')
                if len(dat.shape) == 2:
                    self.var[var] = dat[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                elif len(dat.shape) == 3:
                    self.var[var] = dat[:, self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                elif len(dat.shape) == 4:
                    self.var[var] = dat[:, :, self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                else:
                    print('>4D variable reading not yet implemented')
                    pass
            else:
                self.var[var] = dat

        else:
            self.var[var] = dat

        #read long name and units from netcdf
        self.long_name[var] = ds.variables[var].long_name
        self.units[var] = ds.variables[var].units
        
    def get_modeloutput(self, var):
        """_summary_

        Args:
            var (_type_): _description_

        Returns:
            _type_: _description_
        """        
        self.load_modeloutput(var)
        return self.var[var]

    def get_modeloutput_by_station(self, var, station, isedlayer=None):
        """_summary_

        Args:
            var (_type_): _description_
            station (_type_): _description_

        Returns:
            _type_: _description_
        """
        if len(self.var) == 0:
            self.load_output_coordinates()

        truelist = [station in x for x in self.var['station_id']]
        index = next((i for i, e in enumerate(truelist) if e), None)
        assert index is not None, 'station not found in output'

        self.load_modeloutput(var)

        #in case the modelrun is not yet finished, make sure the time and his len is equal
        lent = len(self.var['pointtime'])
        lentvar = len(self.var[var][:, 0])
        Nt = np.min([lent, lentvar])

        if len(self.var[var].shape)==2:
            return self.var['pointtime'][:Nt-1], self.var[var][:Nt-1, index]
        elif len(self.var[var].shape)==3:
            if isedlayer==None:
                print('no isedlayer specified, using first layer')
                isedlayer=0
            return self.var['pointtime'][:Nt-1], self.var[var][:Nt-1, isedlayer, index]

    def _mdates_concise_subplot_axes(self, ax):
        """_summary_

        Args:
            ax (_type_): _description_
        """        
        major_locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
        formatter = mdates.ConciseDateFormatter(major_locator)
        # see whether we can iterate through ax, if not than cast in list so we can
        try:
            [axi.xaxis.set_major_formatter(formatter) for axi in ax]
        except:
            ax.xaxis.set_major_formatter(formatter)

    def fig_check_tide_bc(self):
        """_summary_

        Returns:
            _type_: _description_
        """        

        assert len(self.AOI) == 0, 'can only check the tide if zs0 is loaded on the entire model domain, so without AOI'
        assert 'zs0file' in self.params, 'No tidal signal available'

        # get model output
        self.load_modeloutput('zs_mean')
        zs = self.var['zs_mean']

        # get model input
        if self.tide == {}:
            self.get_tide()

        zs0_tide = self.tide['zs0']
        if self.globalstarttime is None:
            t_tide = self.tide['time'] / 3600
            t = self.var['meantime'] / 3600
        else:
            t_tide = self.tide['time']
            t = self.var['meantime']

        lent = min([len(zs), len(t)])
        tideloc = self.params['tideloc']

        fig, ax = plt.subplots()
        plt.plot(t[:lent-1], zs[:lent-1, 0, 0], label='xb (1, 1)')
        plt.plot(t_tide, zs0_tide[:, 0], linestyle=':', label='bc (1, 1)')

        if tideloc == 2:
            if 'paulrevere' in self.params:
                if self.params['paulrevere'] == 1:  # two sea conditions, no land conditions
                    plt.plot(t[:lent-1], zs[:lent-1, -1, 0], label='xb (ny, 1)')
                    plt.plot(t_tide, zs0_tide[:, 1], linestyle=':', label='bc (ny, 1)')
                else:  # one offshore condition, one land condition
                    plt.plot(t[:lent-1], zs[:lent-1, 0, -1], label='xb (1,nx)')
                    plt.plot(t_tide, zs0_tide[:, 1], linestyle=':', label='bc (1, nx)')
            else:  # one offshore condition, one land condition
                plt.plot(t, zs[:, 0, -1], label='xb (1,nx)')
                plt.plot(t_tide, zs0_tide[:, 1], linestyle=':', label='bc (1, nx)')
        elif tideloc == 4:
            plt.plot(t[:lent-1], zs[:lent-1, -1, 0], label='xb (ny, 1)')
            plt.plot(t_tide, zs0_tide[:, 1], linestyle=':', label='bc (ny, 1)')
            plt.plot(t[:lent-1], zs[:lent-1, -1, -1], label='xb (ny, nx)')
            plt.plot(t_tide, zs0_tide[:, 2], linestyle=':', label='bc (ny, nx)')
            plt.plot(t[:lent-1], zs[:lent-1, 0, -1], label='xb (1, nx)')
            plt.plot(t_tide, zs0_tide[:, 3], linestyle=':', label='bc (1, nx)')

        plt.legend()
        plt.xlabel('t [hr]')
        plt.ylabel('zs offshore')

        if self.globalstarttime is not None:
            # Define the date format
            ax.xaxis.set_major_formatter(DateFormatter("%m-%d"))
            ax.xaxis.set_minor_formatter(DateFormatter("%H:%M"))
            # Ensure a major tick for each week using (interval=1) and minor tick each 6 hours
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_minor_locator(IndexLocator(base=6 / 24, offset=0))

        plt.tight_layout()
        if self.save_fig:
            
            plt.savefig(os.path.join(self.model_path, 'fig', 'wl_bc_check.png'), dpi=200)
        return fig, ax
    
    def _fig_map_var(self, dat, label, figsize=None, figax = None, **kwargs):
        """_summary_

        Args:
            dat (_type_): _description_
            label (_type_): _description_
            figsize (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """        

        if self.plot_localcoords:
            x = self.var['localx']
            y = self.var['localy']
        else:
            x = self.var['globalx']
            y = self.var['globaly']

        if figax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig, ax = figax
        
        im = ax.pcolor(x, y, dat, **kwargs)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        if np.max(np.abs(dat))<0.1:
            fmt = lambda x, pos: '{:.1e}'.format(x)
        else:
            fmt = lambda x, pos: '{:.1f}'.format(x)

        if 'vmin' in kwargs and 'vmax' in kwargs:
            fig.colorbar(im, cax=cax, orientation='vertical', label=label, extend='both', format=FuncFormatter(fmt))
        elif 'vmax'in kwargs:
            fig.colorbar(im, cax=cax, orientation='vertical', label=label, extend='max', format=FuncFormatter(fmt))
        else:
            fig.colorbar(im, cax=cax, orientation='vertical', label=label, format=FuncFormatter(fmt))

        if self.plot_localcoords:
            if self.plot_km_coords:
                ax.set_xlabel('cross shore [km]')
                ax.set_ylabel('along shore [km]')
            else:
                ax.set_xlabel('cross shore [m]')
                ax.set_ylabel('along shore [m]')
        else:
            if self.plot_km_coords:
                ax.set_xlabel('x [km]')
                ax.set_ylabel('y [km]')
            else:
                ax.set_xlabel('x [m]')
                ax.set_ylabel('y [m]')

        ax.set_aspect('equal')
        fig.tight_layout()

        return fig, ax

    def fig_map_var(self, var, label=None, it=np.inf, figsize=None, figax=None, japlot_mpi_boundaries = False, **kwargs):
        """_summary_

        Args:
            var (_type_): _description_
            label (_type_): _description_
            it (_type_, optional): _description_. Defaults to np.inf.
            figsize (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        self.load_modeloutput(var)

        if np.isinf(it):
            it = len(self.var['globaltime']) - 1
        assert it <= len(self.var['globaltime']) - 1, 'it should be <= {}'.format(len(self.var['globaltime']) - 1)

        # data = self.var[var][it, :, :]

        if (self.plot_localcoords is False) or (not (var in self.vector_vars)):
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                data = self.var[var][it, itype, :, :]
            else:
                data = self.var[var][it, :, :]

        else:
            print('this variable is rotated to gridori in the postprocessing scripts')
            dat1, dat2 = self._get_var_vector_pair(var)
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                dat1 = dat1[it, itype, :, :]
                dat2 = dat2[it, itype, :, :]
                
            else:
                dat1 = dat1[it, :, :]
                dat2 = dat2[it, :, :]

            data_cross, data_along = rotate_grid(dat1.data, dat2.data, self.var['gridang'])

            # select the cross- or alongshore component
            if 'u' in var: 
                data = np.where(~dat1.mask, data_cross, np.nan)  # make 0d again after rotation operation
            else:
                data = np.where(~dat1.mask, data_along, np.nan)  # make 0d again after rotation operation    

            if label is None:
                try:
                    label = str(var)  + ' [' + self.units[var] + ']'
                except:
                    print('unit of {} missing in units'.format(var))
                    label = str(var)


        fig, ax = self._fig_map_var(data, label, figsize, figax=figax, **kwargs)
        
        time_unit = "sec"

        if japlot_mpi_boundaries:
            for iy in self.mpi_iy:
                ax.plot(self.grd['x'][int(iy)-1, :], self.grd['y'][int(iy)-1, :], 'k', linewidth=0.5)
            for ix in self.mpi_ix:
                ax.plot(self.grd['x'][:,int(ix)-1], self.grd['y'][:,int(ix)-1], 'k', linewidth=0.5)     

        if self.globalstarttime is None:
            ax.set_title('{:.1f} {}'.format(self.var['globaltime'][it], time_unit))
        else:
            ax.set_title('{}'.format(self.var['globaltime'][it]))

        plt.tight_layout()
        if self.save_fig:
                
            folder = os.path.join(self.model_path, 'fig', 'map_{}'.format(var[0]))
            if not os.path.exists(folder):
                os.mkdir(folder)
            plt.savefig(os.path.join(folder, 'map_{}_it_{}.png'.format(var, it)), dpi=200)

        fig.tight_layout()

        return fig, ax

    def _get_var_vector_pair(self, var):

        if var in ['ue', 've']: 
            dat1 = self.get_modeloutput('ue')
            dat2 = self.get_modeloutput('ve')
        elif var in ['u', 'v']:
            dat1 = self.get_modeloutput('u')
            dat2 = self.get_modeloutput('v')
        elif var in ['Subg', 'Svbg']:
            dat1 = self.get_modeloutput('Subg')
            dat2 = self.get_modeloutput('Svbg')
        elif var in ['Susg', 'Svsg']:
            dat1 = self.get_modeloutput('Susg')
            dat2 = self.get_modeloutput('Svsg')

        elif var in ['ue_mean', 've_mean']: 
            dat1 = self.get_modeloutput('ue_mean')
            dat2 = self.get_modeloutput('ve_mean')
        elif var in ['u_mean', 'v_mean']:
            dat1 = self.get_modeloutput('u_mean')
            dat2 = self.get_modeloutput('v_mean')
        elif var in ['Subg_mean', 'Svbg_mean']:
            dat1 = self.get_modeloutput('Subg_mean')
            dat2 = self.get_modeloutput('Svbg_mean')
        elif var in ['Susg_mean', 'Svsg_mean']:
            dat1 = self.get_modeloutput('Susg_mean')
            dat2 = self.get_modeloutput('Svsg_mean')

        elif var in ['ue_max', 've_max']: 
            dat1 = self.get_modeloutput('ue_max')
            dat2 = self.get_modeloutput('ve_max')
        elif var in ['u_max', 'v_max']:
            dat1 = self.get_modeloutput('u_max')
            dat2 = self.get_modeloutput('v_max')
        elif var in ['Subg_max', 'Svbg_max']:
            dat1 = self.get_modeloutput('Subg_max')
            dat2 = self.get_modeloutput('Svbg_max')
        elif var in ['Susg_max', 'Svsg_max']:
            dat1 = self.get_modeloutput('Susg_max')
            dat2 = self.get_modeloutput('Svsg_max')

        if var in ['ue_min', 've_min']: 
            dat1 = self.get_modeloutput('ue_min')
            dat2 = self.get_modeloutput('ve_min')
        elif var in ['u_min', 'v_min']:
            dat1 = self.get_modeloutput('u_min')
            dat2 = self.get_modeloutput('v_min')
        elif var in ['Subg_min', 'Svbg_min']:
            dat1 = self.get_modeloutput('Subg_min')
            dat2 = self.get_modeloutput('Svbg_min')
        elif var in ['Susg_min', 'Svsg_min']:
            dat1 = self.get_modeloutput('Susg_min')
            dat2 = self.get_modeloutput('Svsg_min')

        return dat1, dat2
    
    def fig_map_diffvar(self, var, label, it0=0, itend=np.inf, clim=None, figsize=None, **kwargs):
        """_summary_

        Args:
            var (_type_): _description_
            label (_type_): _description_
            it0 (int, optional): _description_. Defaults to 0.
            itend (_type_, optional): _description_. Defaults to np.inf.
            clim (_type_, optional): _description_. Defaults to None.
            figsize (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """        
        assert itend > it0, 'itend must be larger than it0'
        assert it0 >= 0, 'it0 should be >= 0'

        if clim is None:
            kwargs = {'cmap': 'RdBu', 'norm': colors.CenteredNorm()}
        else:
            kwargs = {'cmap': 'RdBu', 'vmin': clim[0], 'vmax': clim[1]}

        self.load_modeloutput(var)

        if np.isinf(itend):
            itend = len(self.var['globaltime']) - 1
        assert itend <= len(self.var['globaltime']) - 1, 'itend should be <= {}'.format(len(self.var['globaltime']) - 1)

        var0 = self.var[var][it0, :, :]
        varend = self.var[var][itend, :, :]

        if (self.plot_localcoords is False) or (not (var in self.vector_vars)):
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                var0 = self.var[var][it0, itype, :, :]
                varend = self.var[var][itend, itype, :, :]
            else:
                var0 = self.var[var][it0, :, :]
                varend = self.var[var][itend, :, :]

        else:
            print('this variable is rotated to gridori in the postprocessing scripts')
            dat1, dat2 = self._get_var_vector_pair(var)
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                var01 = dat1[it0, itype, :, :]
                varend1 = dat1[itend, itype, :, :]    

                var02 = dat2[it0, itype, :, :]
                varend2 = dat2[itend, itype, :, :]    
                
            else:
                var01 = dat1[it0, :, :]
                varend1 = dat1[itend, :, :]    

                var02 = dat2[it0, :, :]
                varend2 = dat2[itend, :, :]    

            data_cross1, data_along1 = rotate_grid(var01.data, var02.data, self.var['gridang'])
            data_cross2, data_along2 = rotate_grid(varend1.data, varend2.data, self.var['gridang'])

            # select the cross- or alongshore component
            if 'u' in var: 
                var0 = np.where(~var01.mask, data_cross1, np.nan).flatten()  # make 0d again after rotation operation
                varend = np.where(~var01.mask, data_cross2, np.nan).flatten()  # make 0d again after rotation operation
            else:
                var0 = np.where(~var01.mask, data_along1, np.nan).flatten()  # make 0d again after rotation operation
                varend = np.where(~var01.mask, data_along2, np.nan).flatten()  # make 0d again after rotation operation
        
        fig, ax = self._fig_map_var(varend - var0, label, figsize, **kwargs)

        #make title
        if self.globalstarttime is None:
            ax.set_title('{:.1f}Hr - {:.1f}Hr'.format(self.var['globaltime'][itend],
                                                  self.var['globaltime'][it0]))
        else:
            ax.set_title('{} - {}'.format(self.var['globaltime'][itend],
                                                      self.var['globaltime'][it0]))

        plt.tight_layout()

        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'difmap_{}_it_{}-{}.png'.format(var, itend, it0)), dpi=200)
        return fig, ax

    def fig_cross_var(self,var, it, iy=None, itype=None, coord=None, plot_ref_bathy=True, remove_dry_points=False, figax = None, zmin=-25, ylim=None, fmt='.-', tight_layout=True):
        """_summary_

        Args:
            var (_type_): _description_
            it (_type_): _description_
            iy (_type_, optional): _description_. Defaults to None.
            coord (_type_, optional): _description_. Defaults to None.
            plot_ref_bathy (bool, optional): _description_. Defaults to True.
            zmin (int, optional): _description_. Defaults to -25.

        Returns:
            _type_: _description_
        """        
        assert not ((iy is None) & (coord is None)), 'specify either an alongshore index iy or a coordinate coord'

        try:
            self.load_modeloutput(var)
        except:
            print('var not found as output specified in params. Will try to continue to see if computed earlier')
        self.load_modeloutput('zb')

        x = self.var['globalx']
        y = self.var['globaly']
        if any(substr in var for substr in ['_mean', '_min', '_max', '_var']):
            t = self.var['meantime'][it]
        else:
            t = self.var['globaltime'][it]

        if iy is None:
            iy, _ = np.unravel_index(((x - coord[0]) ** 2 + (y - coord[1]) ** 2).argmin(), x.shape)

        if (self.plot_localcoords is False) or (not (var in self.vector_vars)):
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                data = self.var[var][it, itype, iy, :]
            else:
                data = self.var[var][it, iy, :]

        else:
            print('this variable is rotated to gridori in the postprocessing scripts')
            dat1, dat2 = self._get_var_vector_pair(var)
            if len(self.var[var].shape)==4:
                if itype is None:
                    print('using first sed type for printing. For others, set itype>=1')
                    itype = 0
                dat1 = dat1[it, itype, iy, :]
                dat2 = dat2[it, itype, iy, :]
                
            else:
                dat1 = dat1[it, iy, :]
                dat2 = dat2[it, iy, :]

            data_cross, data_along = rotate_grid(dat1.data, dat2.data, self.var['gridang'])

            # select the cross- or alongshore component
            if 'u' in var: 
                data = np.where(~dat1.mask, data_cross, np.nan).flatten()  # make 0d again after rotation operation
            else:
                data = np.where(~dat1.mask, data_along, np.nan).flatten()  # make 0d again after rotation operation       
            
        cross = self.var['cross']     

        if remove_dry_points:
            if any(substr in var for substr in ['_mean', '_min', '_max', '_var']):
                try:
                    self.load_modeloutput('zs_min')
                    self.load_modeloutput('zb_mean')
                    data = np.where((self.var['zs_min'][it, iy, :]-self.var['zb_mean'][it, iy, :]).flatten()>=0.01, data, np.nan)
                except:
                    print('zs_min or zb_mean not saved on file, so no dry points are removed from the plot')
            else:
                try:
                    self.load_modeloutput('zs')
                    self.load_modeloutput('zb')
                    data = np.where((self.var['zs'][it, iy, :]-self.var['zb'][it, iy, :]).flatten()>=0.01, data, np.nan)    
                except:
                    print('zs or zb not saved on file, so no dry points are removed from the plot')            
                                
            

        if figax is None:
            fig, ax1 = plt.subplots(figsize=[5, 3])
            ax1.plot(cross, data, 'k.-')
        else:
            fig, ax1 = figax[0], figax[1]
            ax1.plot(cross, data, fmt)
       

        if self.plot_km_coords:
            ax1.set_xlabel('cross shore [km]')
        else:
            ax1.set_xlabel('cross shore [m]')

        if plot_ref_bathy:

            if '_mean' in var:
                itz = np.where(self.var['globaltime']>=t)[0][0]
            else:
                itz = it

            z = self.var['zb'][itz, iy, :]

            ax2 = ax1.twinx()
            ax1.set_zorder(ax2.get_zorder() + 1)  # move ax in front
            ax1.patch.set_visible(False)

            ax2.fill_between(cross, z, zmin, color='lightgrey')
            ax2.set_ylim([zmin, z.max() + 1])
            ax2.set_ylabel('$z_b$ [m+NAP]', color='lightgrey')

            ax = [ax1, ax2]
        else:
            ax = ax1

        ax1.set_xlim([cross.min(), cross.max()])

        if not ylim is None:
            ax1.set_ylim(ylim)

        ax1.set_ylabel(var + ' [' + self.units[var] + ']', color='k')
        ax1.set_title('time: {}'.format(t))
        if tight_layout:
            plt.tight_layout()

        if self.save_fig:
            folder = os.path.join(self.model_path, 'fig', '{}_iy{}'.format(var, iy))
            if not os.path.exists(folder):
                os.mkdir(folder)
            plt.savefig(os.path.join(folder, '{}_iy{}_it{}.png'.format(var, iy, it)),
                        dpi=200)

        return fig, ax
        
    def set_var(self, varname, units, var):
        """_summary_

        Args:
            varname (_type_): _description_
            var (_type_): _description_
        """
        self.var[varname] = var
        self.units[varname] = units
        return

    def fig_profile_change(self, it=-1, iy=None, coord=None, zmin=-25, zmax=12):
        """_summary_

        Args:
            iy (_type_, optional): _description_. Defaults to None.
            coord (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """        
        if iy is None:
            assert coord is not None, 'if no iy index is specified, a coordinate needs to be specified (xi,yi)'
            iy, _ = np.unravel_index(((self.var['globalx'][:] - coord[0]) ** 2 + (self.var['globaly'][:] - coord[1]) ** 2).argmin(), self.var['globalx'].shape)


        zs = self.get_modeloutput('zs')
        zb = self.get_modeloutput('zb')
        cross = self.var['cross']  # because we are sure that after getting the above three variables this one is initialized

        # only load the ne layer if one is in place
        if 'struct' in self.params: # check if struct is in params
            if int(self.params['struct']) == 1:
                self.load_grid()
                if len(self.AOI) > 0:
                    ne = self.grd['ne'][self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                else:
                    ne = self.grd['ne']
                ne = zb[0, :, :]-ne

        fig, ax = plt.subplots()
        plt.plot(cross, np.nanmax(zs, axis=0)[iy, :], color='blue', label='zs-max')
        plt.plot(cross, np.nanmin(zs, axis=0)[iy, :], color='royalblue', label='zs-min')
        plt.plot(cross, zb[0, iy, :], color='k', label='pre')
        plt.plot(cross, zb[it, iy, :], color='r', label='post')

        if 'struct' in self.params: # check if struct is in params
            if int(self.params['struct']) == 1:
                plt.fill_between(cross, zb[0, iy, :], ne[iy, :], color='lightgrey', label='erodible')
                plt.fill_between(cross, ne[iy, :], zmin, color='grey', label='non-erodible')
        else:
            plt.fill_between(cross, zb[0, iy, :], zmin, color='lightgrey', label='erodible')

        plt.title('profile iy = {}'.format(int(iy)))
        plt.legend()
        plt.xlabel('cross-shore [m]')
        plt.ylabel('[m+NAP]')
        plt.xlim([cross[0], cross[-1]])
        plt.ylim([zmin, zmax])
        plt.grid(linestyle=':', color='grey', linewidth=0.5)

        plt.tight_layout()
        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'profile_change_iy{}.png'.format(iy)), dpi=200)

        return fig, ax

    def fig_map_quiver(self, var=['ue', 've'], label='ue [m/s]', it=np.inf, streamspacing=50,
                       figsize=None, vmax=None, vmin=None, ifrac = 0, **kwargs):
        """plots map plots of map output, only works for rectilinear grids (that can be of varying grid resolution).
        Does not work for curvilinear grids!

        Args:
            var (list, optional): _description_. Defaults to ['ue', 've'].
            label (str, optional): _description_. Defaults to 'ue [m/s]'.
            it (_type_, optional): _description_. Defaults to np.inf.
            streamspacing (int, optional): spacing between quiver arrows in meters. Defaults to 50.
            figsize (_type_, optional): _description_. Defaults to None.
            vmax (_type_, optional): _description_. Defaults to None.
            vmin (_type_, optional): _description_. Defaults to None.
            ifrac (int, optional): _description_. Defaults to 0.

        Returns:
            _type_: _description_
        """        
        # toolbox specific import
        from .general.geometry import rotate_grid
            
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        ja_plot_localcoords = self.plot_localcoords

        if ja_plot_localcoords is False:
            print('can only plot quiver in local coordinates so changing the setting of ja_plot_localcoords')
            self.set_plot_localcoords(True)

        self.load_modeloutput(var[0])
        self.load_modeloutput(var[1])

        if self.plot_km_coords:
            # reduce the streamspacing to km's internally
            streamspacing = streamspacing/1e3
        
        if np.isinf(it):
            it = len(self.var['globaltime']) - 1
        assert it <= len(self.var['globaltime']) - 1, 'it should be <= {}'.format(len(self.var['globaltime']) - 1)


        x = np.flipud(self.var['localx'].data)
        y = np.flipud(self.var['localy'].data)
        if len(self.var[var[0]].shape)==3:
            u = np.flipud(self.var[var[0]][it, :, :].data)
            v = np.flipud(self.var[var[1]][it, :, :].data)
            u[u < -100] = 0
            v[v < -100] = 0
            data = np.sqrt((self.var[var[0]][it, :, :]) ** 2 + (self.var[var[1]][it, :, :]) ** 2).squeeze()

        elif len(self.var[var[0]].shape)==4:
            u = np.flipud(self.var[var[0]][it, ifrac, :, :].data)
            v = np.flipud(self.var[var[1]][it, ifrac, :, :].data)
            data = np.sqrt((self.var[var[0]][it, ifrac, :, :]) ** 2 + (self.var[var[1]][it, ifrac, :, :]) ** 2).squeeze()

        u, v = rotate_grid(u, v, self.var['gridang'])

        fu = interp2d(x[:, 0], y[0, :], u.T)
        fv = interp2d(x[:, 0], y[0, :], v.T)

        xt = np.arange(x[0, 0], x[-1, 0], streamspacing)
        yt = np.arange(y[0, 0], y[0, -1], streamspacing)
        X,Y = np.meshgrid(xt,yt)
        ut = fu(xt, yt)
        vt = fv(xt, yt)

        if 'color' not in kwargs:
            kwargs['color'] = 'white'
        
        fig, ax = self._fig_map_var(data, label,figsize=figsize, **{'cmap': 'jet', 'vmax': vmax, 'vmin': vmin})
        ax.quiver(X, Y, ut, vt, units='xy', pivot='middle', **kwargs)
        # ax.streamplot(X, Y, ut, vt, color='k')

        fig.tight_layout()
        if self.globalstarttime is None:
            ax.set_title('t = {:.1f}Hr'.format(self.var['globaltime'][it] ))
        else:
            ax.set_title('t = {}'.format(self.var['globaltime'][it]))
        if self.save_fig:
            folder = os.path.join(self.model_path, 'fig', 'map_quiver_{}'.format(var[0]))
            if not os.path.exists(folder):
                os.mkdir(folder)
            plt.savefig(os.path.join(folder, 'map_quiver_{}_it_{}.png'.format(var[0], it)), dpi=200)

        if ja_plot_localcoords is False:
            self.set_plot_localcoords(False)
        return fig, ax
    
    def fig_map_contour_var(self, var, levels, label=None, it=np.inf, figax=None, figsize=None, **kwargs):
        """
        Generates a contour plot of a specified variable at a given time step.

        Parameters:
        - var: str
            The variable name to plot.
        - levels: array-like
            The levels at which to draw the contour lines.
        - label: str, optional
            The label for the color bar and plot title. Defaults to the variable name.
        - it: int, optional
            The time index for the plot. Defaults to the last time step if np.inf.
        - figax: tuple, optional
            A tuple containing the figure and axis objects. If None, a new figure and axis are created.
        - figsize: tuple, optional
            The size of the figure. Used only if figax is None.
        - kwargs: dict, optional
            Additional keyword arguments passed to plt.contour.

        Returns:
        - fig: matplotlib.figure.Figure
            The figure object containing the plot.
        - ax: matplotlib.axes.Axes
            The axis object containing the plot.
        """

        # Define a nested class for formatting contour levels
        # This class is used to format contour level labels, ensuring that trailing zeros are removed 
        # from the string representation of the level if present.
        class nf(float):
            def __repr__(self):
                s = f'{self:.1f}'
                return f'{self:.0f}' if s[-1] == '0' else s

        # Load the model output data for the specified variable
        self.load_modeloutput(var)

        # Set the time index to the last time step if not specified
        if np.isinf(it):
            it = len(self.var['globaltime']) - 1
        assert it <= len(self.var['globaltime']) - 1, 'it should be <= {}'.format(len(self.var['globaltime']) - 1)

        # Extract the data for the specified variable and time step
        data = self.var[var][it, :, :]
        if label is None:
            label = str(var)

        # Select coordinate system
        if self.plot_localcoords:
            x = self.var['localx']
            y = self.var['localy']
        else:
            x = self.var['globalx']
            y = self.var['globaly']

        # Create a new figure and axis if not provided
        if figax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig, ax = figax

        # Plot the contour
        CS = ax.contour(x, y, data, levels=levels, **kwargs)
        
        # Format the contour levels
        CS.levels = [nf(val) for val in CS.levels]
        ax.clabel(CS, CS.levels, inline=True, fontsize=10)

        # Set axis labels based on coordinate system
        if self.plot_localcoords:
            if self.plot_km_coords:
                ax.set_xlabel('along shore [km]')
                ax.set_ylabel('cross shore [km]')
            else:
                ax.set_xlabel('along shore [m]')
                ax.set_ylabel('cross shore [m]')
        else:
            if self.plot_km_coords:
                ax.set_xlabel('x [km]')
                ax.set_ylabel('y [km]')
            else:
                ax.set_xlabel('x [m]')
                ax.set_ylabel('y [m]')

        ax.set_aspect('equal')
        fig.tight_layout()

        # Set the plot title
        if self.globalstarttime is None:
            ax.set_title(f'{var} - t = {self.var["globaltime"][it]:.1f} Hr')
        else:
            ax.set_title(f'{var} - t = {self.var["globaltime"][it]}')

        # Save the figure if required
        if self.save_fig:
            folder = os.path.join(self.model_path, 'fig', f'map_{var[0]}')
            if not os.path.exists(folder):
                os.mkdir(folder)
            plt.savefig(os.path.join(folder, f'map_{var}_it_{it}.png'), dpi=200)

        fig.tight_layout()

        return fig, ax