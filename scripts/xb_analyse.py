import numpy as np
import os
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors


class XBeachModelAnalysis():
    '''
    XBeach model analysis class
    '''

    def __init__(self, fname, model_path):
        self.fname = fname
        self.model_path = model_path
        self.params = None
        self.grd = {}
        self.waves_boundary = {}
        self.dat = {}
        self.tide = {}
        self.var = {}
        self.save_fig = False
        self.plot_localcoords = False
        self.plot_km_coords = False
        self.AOI = []

    def __repr__(self):
        return self.fname

    def set_save_fig(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.save_fig = yesno

    def set_plot_localcoords(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_localcoords = yesno

    def set_plot_km_coords(self, yesno):
        assert type(yesno) is bool, 'input type must be bool'
        self.plot_km_coords = yesno

    def set_aoi(self, AOI):
        self.var = {}  # drop variables from memory because they need to be reloaded with appropriate AOI
        self.AOI = AOI

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
        i0 = [i for i, var in enumerate(dat) if 'nglobalvar' in var][0]
        params['globalvar'] = dat[i0+1:i0+int(params['nglobalvar']+1)]
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
                self.waves_boundary['tíme'] = np.cumsum(dat[:, 5])

            elif self.params['wbctype'] == 'jons':
                print('not yet written')
                pass
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

    def load_model_setup(self):
        self.get_params()
        self.load_grid()
        self.get_waves()
        self.get_vegetation()
        self.get_tide()

    def load_output_coordinates(self):

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))

        self.var['globaltime'] = ds.variables['globaltime'][:]
        if self.params['nmeanvar'] > 0:
            self.var['meantime'] = ds.variables['meantime'][:]
        if self.params['npointvar'] > 0:
            self.var['pointtime'] = ds.variables['pointtime'][:]

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
            pass

        elif 'mean' in var:
            assert sum([var[5:] in x for x in self.params['meanvar']]) > 0, '{} not in xb output'
            pass
        elif 'point' in var:
            assert sum([var[6:] in x for x in self.params['pointvar']]) > 0, '{} not in xb output'
            pass
        else:
            assert sum([var in x for x in self.params['globalvar']]) > 0, '{} not in xb output'
            pass

        # if not yet present, load coordinates
        if self.var == {}:
            print('loading model output coordinates from file')
            self.load_output_coordinates()

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))
        print('loading variable {} from file'.format(var))
        dat = ds.variables[var][:]

        if not ('point' in var):
            if len(self.AOI) > 0:
                print('slicing map output to AOI')
                if len(dat.shape) == 2:
                    self.var[var] = dat[self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                elif len(dat.shape) == 3:
                    self.var[var] = dat[:, self.AOI[0]:self.AOI[1], self.AOI[2]:self.AOI[3]]
                else:
                    print('4D variable reading not yet implemented')
                    # todo: >3D variable reading implementing
                    pass
            else:
                self.var[var] = dat
        else:
            self.var[var] = dat

    def get_modeloutput(self, var):
        self.load_modeloutput(var)
        return self.var[var]

    def get_modeloutput_by_station(self, var, station):

        truelist = [station in x for x in self.var['station_id']]
        index = next((i for i, e in enumerate(truelist) if e), None)
        assert index is not None, 'station not found in output'

        self.load_modeloutput(var)

        return np.ma.vstack([self.var['pointtime'], self.var[var][:, index]]).T

    def fig_check_tide_bc(self):

        assert len(self.AOI) == 0, 'can only check the tide if zs0 is loaded on the entire model domain, so without AOI'

        # get model input
        if self.tide == {}:
            self.get_tide()
        t_tide = self.tide['time']
        zs0_tide = self.tide['zs0']

        # get model output
        self.load_modeloutput('zs')
        zs = self.var['zs']
        t = self.var['globaltime'] / 3600
        tideloc = self.params['tideloc']

        fig, ax = plt.subplots()
        plt.plot(t, zs[:, 0, 0], label='xb (1, 1)')
        plt.plot(t_tide / 3600, zs0_tide[:, 0], linestyle=':', label='bc (1, 1)')

        if tideloc == 2:
            if 'paulrevere' in self.params:
                if self.params['paulrevere'] == 1:  # two sea conditions, no land conditions
                    plt.plot(t, zs[:, -1, 0], label='xb (ny, 1)')
                    plt.plot(t_tide / 3600, zs0_tide[:, 1], linestyle=':', label='bc (ny, 1)')
                else:  # one offshore condition, one land condition
                    plt.plot(t, zs[:, 0, -1], label='xb (1,nx)')
                    plt.plot(t_tide / 3600, zs0_tide[:, 1], linestyle=':', label='bc (1, nx)')
            else:  # one offshore condition, one land condition
                plt.plot(t, zs[:, 0, -1], label='xb (1,nx)')
                plt.plot(t_tide / 3600, zs0_tide[:, 1], linestyle=':', label='bc (1, nx)')
        elif tideloc == 4:
            plt.plot(t, zs[:, -1, 0], label='xb (ny, 1)')
            plt.plot(t_tide / 3600, zs0_tide[:, 1], linestyle=':', label='bc (ny, 1)')
            plt.plot(t, zs[:, -1, -1], label='xb (ny, nx)')
            plt.plot(t_tide / 3600, zs0_tide[:, 2], linestyle=':', label='bc (ny, nx)')
            plt.plot(t, zs[:, 0, -1], label='xb (1, nx)')
            plt.plot(t_tide / 3600, zs0_tide[:, 3], linestyle=':', label='bc (1, nx)')

        plt.legend()
        plt.xlabel('t [hr]')
        plt.ylabel('zs offshore')
        if self.save_fig:
            plt.savefig(os.path.join(self.model_path, 'fig', 'wl_bc_check.png'), dpi=200)
        return fig, ax

    def _fig_map_var(self, dat, label, **kwargs):

        if self.plot_localcoords:
            x = self.var['localx']
            y = self.var['localy']
        else:
            x = self.var['globalx']
            y = self.var['globaly']

        fig, ax = plt.subplots()
        im = ax.pcolor(x, y, dat, **kwargs)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical', label=label)

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
        ax.set_title('t = {:.1f}Hr'.format(self.var['globaltime'][it] / 3600))

        return fig, ax

    def fig_map_diffvar(self, var, label, it0=0, itend=np.inf):
        assert itend > it0, 'itend must be larger than it0'
        assert it0 >= 0, 'it0 should be >= 0'

        self.load_modeloutput(var)

        if np.isinf(itend):
            itend = len(self.var['globaltime']) - 1
        assert itend <= len(self.var['globaltime']) - 1, 'itend should be <= {}'.format(len(self.var['globaltime']) - 1)

        var0 = self.var['zb'][it0, :, :]
        varend = self.var['zb'][itend, :, :]

        fig, ax = self._fig_map_var(varend - var0, label, **{'cmap': 'RdBu', 'norm': colors.CenteredNorm()})
        ax.set_title('{:.1f}Hr - {:.1f}Hr'.format(self.var['globaltime'][itend] / 3600,
                                                  self.var['globaltime'][it0] / 3600))

        return fig, ax

    def plot_profile(self):
        # plot cross-sections
        x1, y1 = 117421.6, 560053.6
        x2, y2 = 115469.4, 558176.2
        iy1, ix1 = np.unravel_index(((x - x1) ** 2 + (y - y1) ** 2).argmin(), x.shape)
        iy2, ix2 = np.unravel_index(((x - x2) ** 2 + (y - y2) ** 2).argmin(), x.shape)
        ne = zb[0, :, :] - np.loadtxt(rundir + 'ne.txt')  # see erodible layer

    def fig_profile_change(self, iy=None, coord=None):
        if iy is None:
            assert coord is not None, 'if no iy index is specified, a coordinate needs to be specified (xi,yi)'

        zs = self.get_modeloutput('zs')
        zb = self.get_modeloutput('zb')
        cross = self.var['zs']  # because we are sure that after getting the above three variables this one is initialized

        # only load the ne layer if one is in place
        if int(self.params['struct']) == 1:
            self.load_grid()
            ne = self.grid['zb']

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
        return fig