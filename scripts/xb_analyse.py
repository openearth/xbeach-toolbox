import numpy as np
import os
import netCDF4 as nc

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

    def __repr__(self):
        return self.fname

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
        i0 = [i for i, var in enumerate(dat) if 'nmeanvar' in var][0]
        params['meanvar'] = dat[i0+1:i0+int(params['nmeanvar']+1)]
        # point variables
        ixlist = [i for i, var in enumerate(dat) if 'npointvar' in var]
        if len(ixlist)>0:
            i0 = ixlist[0]
            params['pointvar'] = dat[i0+1:i0+int(params['npointvar']+1)]
        else:
            params['pointvar'] = []


        # output points
        if 'npointvar' in params:
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
        self.var['x'] = ds.variables['globalx'][:]
        self.var['y'] = ds.variables['globaly'][:]


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
        if self.var is {}:
            self.read_output_coordinates()

        ds = nc.Dataset(os.path.join(self.model_path, 'xboutput.nc'))
        self.var[var] = ds.variables[var][:]




