import numpy as np
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


class XBeachModelAnalysis():
    '''
    XBeach model analysis class
    '''

    def __init__(self, fname, model_path):
        self.fname = fname
        self.model_path = model_path
        self.get_params()
        # self.params = None
        self.grd = {}
        self.waves_boundary = {}
        self.dat = {}
        self.tide = {}
        self.wind = {}
        self.var = {}
        self.save_fig = False
        self.plot_localcoords = False
        self.plot_km_coords = False
        self.AOI = []
        self.globalstarttime = None

    def __repr__(self):
        return self.fname

    def set_save_fig(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.save_fig = yesno
        plotdir = os.path.join(self.model_path, 'fig')
        if not os.path.isdir(plotdir):
            os.mkdir(plotdir)

    def set_plot_localcoords(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_localcoords = yesno

    def set_plot_km_coords(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_km_coords = yesno

    def set_aoi(self, AOI):
        self.var = {}  # drop variables from memory because they need to be reloaded with appropriate AOI
        self.AOI = AOI

    def set_globalstarttime(self, tstart):
        assert type(tstart) is str, 'tstart must be given as a string of the format 2021-10-11T13:00:00'
        self.globalstarttime = np.datetime64(tstart)

    def get_params(self):
        '''
        todo: better to read these from the xbeach.log instead of the params.txt
        '''

        assert os.path.exists(os.path.join(self.model_path, 'params.txt')), 'missing params.txt'

        f = open(os.path.join(self.model_path, 'params.txt'), 'r')
        dat = f.read().split('\n')

        # read params from params file
        params = {}
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
            x = [float(t[0]) for t in points]
            y = [float(t[1]) for t in points]
            name = [t[2].strip() for t in points]
            params['points'] = dict(zip(name, zip(x, y)))
        else:
            params['points'] = {}

        self.params = params

    def load_grid(self):
        # only works currently for xbeach type grids (xfile, yfile)

        # read params.txt if this is not done yet
        if self.params is None:
            self.get_params()

        # load obligatory files
        self.grd['x'] = np.loadtxt(os.path.join(self.model_path, self.params['xfile']))
        assert self.grd['x'].shape == (self.params['ny']+1, self.params['nx']+1), 'x grid not of correct size'

        self.grd['z'] = np.loadtxt(os.path.join(self.model_path, self.params['depfile']))
        assert self.grd['z'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'z grid not of correct size'
        if 'posdwn' in self.params:
            if int(self.params['posdwn']) == 1:
                self.grd['z'] = -1*self.grd['z']

        # read optional files
        if self.params['ny'] > 0:
            self.grd['y'] = np.loadtxt(os.path.join(self.model_path, self.params['yfile']))
            assert self.grd['y'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'y grid not of correct size'

        if int(self.params['struct']) == 1:
            self.grd['ne'] = np.loadtxt(os.path.join(self.model_path, self.params['ne_layer']))
            assert self.grd['ne'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'ne grid not of correct size'

    def get_waves(self):
        # read params.txt if this is not done yet
        if self.params is None:
            self.get_params()

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
            self.waves_boundary['time'] = np.cumsum(dat[:, 5])

        elif self.params['wbctype'] == 'jons':
            print('not yet written')
            pass
        elif self.params['wbctype'] == 'params':
            self.waves_boundary['Hm0'] = np.sqrt(2)*float(self.params['Hrms'])
            self.waves_boundary['Tp'] = float(self.params['Trep'])
            self.waves_boundary['mainang'] = float( self.params['dir0'])
            self.waves_boundary['s'] = float(self.params['m']/2)
        else:
            print('not yet written')
            pass

    def get_vegetation(self):
        pass

    def get_tide(self):
        if self.params is None:
            self.get_params()

        dat = np.loadtxt(os.path.join(self.model_path, self.params['zs0file']))

        if int(self.params['tideloc']) == 1:
            assert dat.shape[1] == 2, 'tideloc=1, expected 2 cols'
        if int(self.params['tideloc']) == 2:
            assert dat.shape[1] == 3, 'tideloc=2, expected 3 cols'
        if int(self.params['tideloc']) == 4:
            assert dat.shape[1] == 5, 'tideloc=1, expected 5 cols'

        self.tide['time'] = dat[:, 0]
        self.tide['zs0'] = dat[:, 1:]
        self.tide['tideloc'] = self.params['tideloc']
        if 'paulrevere' in self.params:
            self.tide['paulrevere'] = self.params['paulrevere']

    def get_wind(self):
        if self.params is None:
            self.get_params()

        if self.params['wind'] != 1:
            print('no wind forcing was imposed')
            return
        
            
            
        dat = np.loadtxt(os.path.join(self.model_path, self.params['windfile']))

        self.wind['time'] = dat[:, 0]
        self.wind['u10'] = dat[:, 1]
        self.wind['u10dir'] = dat[:, 2]

            
    def load_model_setup(self):
        self.get_params()
        self.load_grid()
        self.get_waves()
        self.get_vegetation()
        self.get_tide()

    def load_output_coordinates(self):

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))

        # global variable time
        if self.globalstarttime is None:
            self.var['globaltime'] = ds.variables['globaltime'][:]
        else:
            self.var['globaltime'] = np.array([np.timedelta64(int(x), 's') for x in ds.variables['globaltime'][:].data]) \
                                      + self.globalstarttime
        # mean variable time
        if self.params['nmeanvar'] > 0:
            if self.globalstarttime is None:
                self.var['meantime'] = ds.variables['meantime'][:]
            else:
                data = ds.variables['meantime'][:]
                data = data[data.mask == False]
                self.var['meantime'] = np.array(
                    [np.timedelta64(int(x), 's') for x in data.data.flatten()]) \
                                         + self.globalstarttime
        # point variable time
        if self.params['npointvar'] > 0:
            if self.globalstarttime is None:
                self.var['pointtime'] = ds.variables['pointtime'][:]
            else:
                data = ds.variables['pointtime'][:]
                data = data[data.mask == False]
                self.var['pointtime'] = np.array(
                    [np.timedelta64(int(x), 's') for x in data.data.flatten()]) \
                                         + self.globalstarttime

            station_list = []
            dat = ds.variables['station_id'][:]
            nb, nstat = dat.shape
            for ib in range(nb) :
                sn = ''
                for istat in range(nstat):
                    sn += dat.data[ib,istat].decode('UTF-8')
                station_list.append(sn)
            station_list = [x.strip() for x in station_list]
            self.var['station_id'] = station_list

        x = ds.variables['globalx'][:]
        y = ds.variables['globaly'][:]
        if len(self.AOI) > 0:
            assert len(self.AOI) == 4, 'AOI should be specified as [x0, xend, y0, yend]'

            x = x[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
            y = y[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]

        if self.plot_km_coords:
            self.var['globalx'] = x/1e3
            self.var['globaly'] = y/1e3
        else:
            self.var['globalx'] = x
            self.var['globaly'] = y

        self.var['gridang'] = np.arctan2(y[0, 0]-y[-1, 0], x[0, 0]-x[-1, 0])

        def path_distance(polx, poly):
            '''
            computes the distance along a polyline
            ----------
            polx : TYPE array
                X COORDINATES.
            poly : TYPE array
                Y COORDINATES.

            Returns: TYPE array
                PATHDISTANCE
            -------
            python alternative to matlab function.

            '''
            dx = np.diff(polx)
            dy = np.diff(poly)

            dr = np.sqrt(dx ** 2 + dy ** 2)

            pathdistance = np.insert(np.cumsum(dr), 0, 0, axis=0)

            return pathdistance

        self.var['cross'] = path_distance(x[0, :], y[0, :])
        self.var['along'] = path_distance(x[:, 0], y[:, 0])

        self.var['localy'], self.var['localx'] = np.meshgrid(self.var['cross'], self.var['along'])
        self.var['localx'] = np.flipud(self.var['localx'])  # to plot



    def load_modeloutput(self, var):
        if var in self.var:
            return

        if '_mean' in var:
            assert sum([var[:-5] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'
        elif '_min' in var:
            assert sum([var[:-4] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'
        elif '_max' in var:
            assert sum([var[:-4] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'
        elif 'point_' in var:
            assert sum([var[6:] in x for x in self.params['pointvar']]) > 0, '{} not in xb output'
        else:
            assert sum([var in x for x in self.params['globalvar']]) > 0, '{} not in xb output'

        # if not yet present, load coordinates
        if self.var == {}:
            print('loading model output coordinates from file')
            self.load_output_coordinates()

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))
        print('loading variable {} from file'.format(var))
        dat = ds.variables[var][:]

        #mean and point output might not be availble if the eor is not reached. Therefore cut
        if 'point_' in var and len(np.atleast_1d(dat.mask)) > 1:
            dat = dat[~dat.mask[:, 0], :]
        elif len(dat.shape) == 3 and '_mean' in var and len(np.atleast_1d(dat.mask)) > 1:
            dat = dat[~dat.mask[:, 0, 0], :, :]
        elif len(dat.shape) == 2 and '_mean' in var and len(np.atleast_1d(dat.mask)) > 1:
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

    def get_modeloutput(self, var):
        self.load_modeloutput(var)
        return self.var[var]

    def get_modeloutput_by_station(self, var, station):

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

        return self.var['pointtime'][:Nt-1], self.var[var][:Nt-1, index]

    def _mdates_concise_subplot_axes(self, ax):
        major_locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
        formatter = mdates.ConciseDateFormatter(major_locator)
        # see whether we can iterate through ax, if not than cast in list so we can
        try:
            [axi.xaxis.set_major_formatter(formatter) for axi in ax]
        except:
            ax.xaxis.set_major_formatter(formatter)


    def fig_check_tide_bc(self):

        assert len(self.AOI) == 0, 'can only check the tide if zs0 is loaded on the entire model domain, so without AOI'


        # get model output
        self.load_modeloutput('zs')
        zs = self.var['zs']



        # get model input
        if self.tide == {}:
            self.get_tide()

        zs0_tide = self.tide['zs0']
        if self.globalstarttime is None:
            t_tide = self.tide['time'] / 3600
            t = self.var['globaltime'] / 3600
        else:
            t_tide = np.array(
                [np.timedelta64(int(x), 's') for x in self.tide['time']]) \
                     + self.globalstarttime
            t = self.var['globaltime']

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

        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'wl_bc_check.png'), dpi=200)
        return fig, ax

    def _fig_map_var(self, dat, label, figsize=None, **kwargs):

        if self.plot_localcoords:
            x = self.var['localx']
            y = self.var['localy']
        else:
            x = self.var['globalx']
            y = self.var['globaly']

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.pcolor(x, y, dat, **kwargs)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        if np.max(np.abs(dat))<0.001:
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

        return fig, ax

    def fig_map_var(self, var, label, it=np.inf):

        self.load_modeloutput(var)

        if np.isinf(it):
            it = len(self.var['globaltime']) - 1
        assert it <= len(self.var['globaltime']) - 1, 'it should be <= {}'.format(len(self.var['globaltime']) - 1)

        data = self.var[var][it, :, :]
        fig, ax = self._fig_map_var(data, label)

        if self.globalstarttime is None:
            ax.set_title('t = {:.1f}Hr'.format(self.var['globaltime'][it] ))
        else:
            ax.set_title('t = {}'.format(self.var['globaltime'][it]))
        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'map_{}_it_{}.png'.format(var,it)), dpi=200)
        return fig, ax

    def fig_map_diffvar(self, var, label, it0=0, itend=np.inf, clim=None):
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

        var0 = self.var['zb'][it0, :, :]
        varend = self.var['zb'][itend, :, :]

        fig, ax = self._fig_map_var(varend - var0, label, **kwargs)

        #make title
        if self.globalstarttime is None:
            ax.set_title('{:.1f}Hr - {:.1f}Hr'.format(self.var['globaltime'][itend],
                                                  self.var['globaltime'][it0]))
        else:
            ax.set_title('{} - {}'.format(self.var['globaltime'][itend],
                                                      self.var['globaltime'][it0]))

        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'difmap_{}_it_{}-{}.png'.format(var, itend, it0)), dpi=200)
        return fig, ax


    def fig_profile_change(self, iy=None, coord=None):
        if iy is None:
            assert coord is not None, 'if no iy index is specified, a coordinate needs to be specified (xi,yi)'

        zs = self.get_modeloutput('zs')
        zb = self.get_modeloutput('zb')
        cross = self.var['cross']  # because we are sure that after getting the above three variables this one is initialized

        # only load the ne layer if one is in place
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
        plt.plot(cross, zb[-1, iy, :], color='r', label='post')

        if int(self.params['struct']) == 1:
            plt.fill_between(cross, zb[0, iy, :], ne[iy, :], color='lightgrey', label='erodible')
            plt.fill_between(cross, ne[iy, :], -25, color='grey', label='non-erodible')
        else:
            plt.fill_between(cross, zb[0, iy, :], -25, color='lightgrey', label='erodible')

        plt.title('profile iy = {}'.format(int(iy)))
        plt.legend()
        plt.xlabel('cross-shore [m]')
        plt.ylabel('[m+NAP]')
        plt.xlim([cross[0], cross[-1]])
        plt.ylim([-25, 12])
        plt.grid(linestyle=':', color='grey', linewidth=0.5)

        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'profile_change_iy{}.png'.format(iy)), dpi=200)

        return fig, ax

    def fig_map_quiver(self, var=['ue', 've'], label='ue [m/s]', it=np.inf, streamspacing=50,
                       figsize=None, vmax=None, vmin=None, ifrac = 0, **kwargs):
        '''
        plots map plots of map output, only works for rectilinear grids (that can be of varying grid resolution).
        Does not work for curvilinear grids!
        '''
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        ja_plot_localcoords = self.plot_localcoords

        if ja_plot_localcoords is False:
            self.set_plot_localcoords(True)

        self.load_modeloutput(var[0])
        self.load_modeloutput(var[1])

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

        u, v = xb.rotate_grid(u, v, self.var['gridang'])

        fu = interp2d(x[:, 0], y[0, :], u.T)
        fv = interp2d(x[:, 0], y[0, :], v.T)

        xt = np.arange(x[0, 0], x[-1, 0], streamspacing)
        yt = np.arange(y[0, 0], y[0, -1], streamspacing)
        X,Y = np.meshgrid(xt,yt)
        ut = fu(xt, yt)
        vt = fv(xt, yt)

        fig, ax = self._fig_map_var(data, label,figsize=figsize, **{'cmap': 'jet', 'vmax': vmax, 'vmin': vmin})
        ax.quiver(X, Y, ut, vt, color='w', units='xy', pivot='middle', **kwargs)
        # ax.streamplot(X, Y, ut, vt, color='k')

        fig.tight_layout()
        if self.globalstarttime is None:
            ax.set_title('t = {:.1f}Hr'.format(self.var['globaltime'][it] ))
        else:
            ax.set_title('t = {}'.format(self.var['globaltime'][it]))
        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'map_quiver_{}_it_{}.png'.format(var[0], it)), dpi=200)

        if ja_plot_localcoords is False:
            self.set_plot_localcoords(False)
        return fig, ax
