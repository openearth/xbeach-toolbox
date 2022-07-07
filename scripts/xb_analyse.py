import numpy as np
import os
import netCDF4 as nc
import matplotlib.pyplot as plt

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

    def __repr__(self):
        return self.fname

    def save_fig(self,yesno):
        self.save_fig = yesno

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
        if len(ixlist)>0:
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
            params['points'] = dict(zip(name,zip(x,y)))
        else:
            params['points'] = {}

        self.params = params

    def get_grid(self):
        ## only works currently for xbeach type grids (xfile, yfile)

        #read params.txt if this is not done yet
        if self.params is None:
            self.get_params()

        # load obligatory files
        self.grd['x'] = np.loadtxt(os.path.join(self.model_path, self.params['xfile']))
        assert self.grd['x'].shape == (self.params['ny']+1, self.params['nx']+1), 'x grid not of correct size'

        self.grd['z'] = np.loadtxt(os.path.join(self.model_path, self.params['depfile']))
        assert self.grd['z'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'z grid not of correct size'
        if 'posdwn' in self.params:
            if int(self.params['posdwn']) == 1:
                self.grid['z'] = -1*self.grid['z']

        #read optional files
        if self.params['ny']>0:
            self.grd['y'] = np.loadtxt(os.path.join(self.model_path, self.params['yfile']))
            assert self.grd['y'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'y grid not of correct size'

        if int(self.params['struct']) == 1:
            self.grd['ne'] = np.loadtxt(os.path.join(self.model_path, self.params['ne_layer']))
            assert self.grd['ne'].shape == (self.params['ny'] + 1, self.params['nx'] + 1), 'ne grid not of correct size'

    def get_waves(self):
        # read params.txt if this is not done yet
        if self.params is None:
            self.get_params()

            ## waves boundary
            if self.params['wbctype'] == 'jonstable':
                dat = np.loadtxt(os.path.join(self.model_path, self.params['bcfile']))

                assert dat.shape[1] == 7, 'columns of jonstable should exactly be: Hm0, Tp, mainang, gammajsp, s, duration, dtbc'

                self.waves_boundary['Hm0'] = dat[:,0]
                self.waves_boundary['Tp'] = dat[:, 1]
                self.waves_boundary['mainang'] = dat[:, 2]
                self.waves_boundary['gammajsp'] = dat[:, 3]
                self.waves_boundary['s'] = dat[:, 4]
                self.waves_boundary['tíme'] = np.cumsum(dat[:,5])

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
        self.get_grid()
        self.get_waves()
        self.get_vegetation()
        self.get_tide()

    def read_output_coordinates(self):

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))

        self.var['globaltime'] = ds.variables['globaltime'][:]
        if self.params['nmeanvar']>0:
            self.var['meantime'] = ds.variables['meantime'][:]
        if self.params['npointvar']>0:
            self.var['pointtime'] = ds.variables['pointtime'][:]
        self.var['x'] = ds.variables['globalx'][:]
        self.var['y'] = ds.variables['globaly'][:]

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

            pathDistance = np.insert(np.cumsum(dr), 0, 0, axis=0)

            return pathDistance

        self.var['cross'] = path_distance(self.var['x'][0, :], self.var['y'][0, :])
        self.var['along'] = path_distance(self.var['x'][:, 0], self.var['y'][:, 0])


    def read_modeloutput(self, var):
        if var in self.var:
            pass

        elif 'mean' in var:
            assert sum([var in x for x in self.params['meanvar']]) > 0, '{} not in xb output'
            pass
        elif 'point' in var:
            assert sum([var in x for x in self.params['pointvar']]) > 0, '{} not in xb output'
            pass
        else:
            assert sum([var in x for x in self.params['globalvar']]) > 0, '{} not in xb output'
            pass

        # if not yet present, load coordinates
        if self.var == {}:
            print('loading model output coordinates from file')
            self.read_output_coordinates()

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))
        print('loading variable {} from file'.format(var))
        self.var[var] = ds.variables[var][:]

    def fig_check_tide_bc(self):

        # get model input
        if self.tide == {}:
            self.get_tide()
        t_tide = self.tide['time']
        zs0_tide = self.tide['zs0']

        # get model output
        self.read_modeloutput('zs')
        zs = self.var['zs']
        t = self.var['globaltime'] / 3600
        tideloc = self.params['tideloc']

        fig, ax = plt.subplots()
        plt.plot(t, zs[:, 0, 0], label='xb (1, 1)')
        plt.plot(t_tide / 3600, zs0_tide[:,0], linestyle=':', label='bc (1, 1)')

        if tideloc == 2:
            if 'paulrevere' in self.params:
                if self.params['paulrevere'] == 1:  # two sea conditions, no land conditions
                    plt.plot(t, zs[:, -1, 0], label='xb (ny, 1)')
                    plt.plot(t_tide / 3600, zs0_tide[:, 1], linestyle=':', label='bc (ny, 1)')
                else:  #one offshore condition, one land condition
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


