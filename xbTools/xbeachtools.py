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
import json
import matplotlib.pyplot as plt

# import rotate grid from geometry
from .general.geometry import rotate_grid
from .general.deg2uv import deg2uv

class XBeachModelSetup():
    '''
    XBeach model setup class
    ''' 
    
    def __init__(self,fname):
        self.fname      = fname
        ## by default set wbctype and wavemodel to None
        self.wbctype    = None
        self.wavemodel  = None
        
        ## set default values
        self.model_path = None
        self.friction_layer = None
        self.wavefriction_layer = None
        self.struct = None
        
    def __repr__(self):
        """_summary_

        Returns:
            _type_: _description_
        """        
        return self.fname
    
    def set_params(self,input_par_dict):
        """_summary_

        Args:
            input_par_dict (_type_): _description_
        """        
        ## set wavemodel. Default is Surfbeat
        if 'Wavemodel' not in input_par_dict:
            print('No wavemodel defined. Wavemodel is set to Surfbeat')
            self.wavemodel = 'surfbeat'
        else:
            self.wavemodel = input_par_dict['Wavemodel']
        
        ## set wbctype
        if 'wbctype' in input_par_dict:
            self.wbctype = input_par_dict['wbctype'] 


        ## load parameters and categories
        f           = open(os.path.join(os.path.dirname(__file__), 'par.json'),'r')
        par_dict    = json.loads(f.read())
        ## create input dict
        self.input_par = {}
        self.input_par['par'] = {}
        ## loop over input parameters 
        value_added = False
        for input_par in input_par_dict:
            ## loop over categories
            for par_category in par_dict:
                ## if input parameter is in category, add parameter
                if input_par in par_dict[par_category]:
                    ## create category if not exist
                    if not par_category in self.input_par:
                        self.input_par[par_category] = {}
                    ## add parameter and value                    
                    self.input_par[par_category][input_par] = input_par_dict[input_par]
                    value_added = True
            if not value_added:
                self.input_par['par'][input_par] = input_par_dict[input_par]
            
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
        ## set values
        self.posdwn = posdwn
        self.xori   = xori
        self.yori   = yori
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
        nebed : input of sandy layer thickness
        struct : optional. The default is 1.

        Returns
        -------
        None.

        '''
        self.friction = friction
        self.friction_layer = friction_layer

    def set_wavefriction(self, wavefriction, wavefriction_layer = 1):      
        '''
        function to set friction layer for the xbeach model

        Parameters
        ----------
        nebed : input of sandy layer thickness
        struct : optional. The default is 1.

        Returns
        -------
        None.

        '''
        self.wavefriction = wavefriction
        self.wavefriction_layer = wavefriction_layer

    def set_waves(self,wbctype, input_struct):
        """_summary_

        Args:
            wbctype (_type_): _description_
            input_struct (_type_): _description_
        """        
        self.wbctype = wbctype
        ##
        if wbctype=='jonstable':
            required_par = ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc']
        elif wbctype=='parametric':
            required_par = ['Hm0','Tp','mainang','gammajsp','s','fnyq']
        else:
            assert False, 'Wrong wbctype'
        
        self.waves_boundary  = {}
        for item in required_par:
            assert item in input_struct, '{} missing'.format(item)
            self.waves_boundary[item] =  input_struct[item]
            
    def set_vegetation(self):
        """_summary_
        """        
        pass

    def set_tide(self):
        """_summary_
        """        
        pass
        
    def load_model_setup(self,path):
        """_summary_

        Args:
            path (_type_): _description_
        """        
        ## todo
        pass    

    def write_model(self, path, figure=True):
        """_summary_

        Args:
            path (_type_): _description_
            figure (bool, optional): _description_. Defaults to True.
        """        
        self.model_path = path
        path_params = os.path.join(path,'params.txt')
        
        assert os.path.exists(path), '{} does not exist'.format(path)
        
        
        current_date    = datetime.today().strftime('%Y-%m-%d %HH:%mm')
        user            =  os.path.basename(os.path.expanduser('~'))
        
        tabnumber = 20
        
        ## waves boundary
        if self.wbctype=='parametric':
            if 'Wave boundary condition parameters' in self.input_par:
                self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonswap.txt'
            else:
               self.input_par['Wave boundary condition parameters'] = {}
               self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonswap.txt'
            required_par = ['Hm0','Tp','mainang','gammajsp','s','fnyq']
            with open(os.path.join(path,'jonswap.txt'),'w') as f:
                for par in required_par:
                    f.write('{}\t= {}\n'.format(par,self.waves_boundary[par]).expandtabs(tabnumber))
                
        elif self.wbctype=='jonstable':
            if 'Wave boundary condition parameters' in self.input_par:
                self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonstable.txt'
            else:
               self.input_par['Wave boundary condition parameters'] = {}
               self.input_par['Wave boundary condition parameters']['bcfile'] = 'jonstable.txt'                
            required_par = ['Hm0','Tp','mainang','gammajsp','s','duration','dtbc']
            with open(os.path.join(path,'jonstable.txt'),'w') as f:
                for ii in range(len(self.waves_boundary['Hm0'])):
                    for par in required_par:
                        f.write('{} '.format(self.waves_boundary[par][ii]))
                    f.write('\n')
        
        
        ## create params
        with open(path_params,'w') as f:
            ## meta data
            f.write('%% XBeach model: {} \n'.format(self.fname))
            f.write('%% Params created on {} \n'.format(current_date))
            f.write('%% Params created by {} \n'.format(user))
            f.write('\n')

            ## general
            f.write('%% General \n')
            f.write('\n')
            if self.wavemodel!=None:
                f.write('wavemodel\t= {}\n'.format(self.wavemodel).expandtabs(tabnumber))
            if self.wbctype!=None:
                f.write('wbctype\t= {}\n'.format(self.wbctype).expandtabs(tabnumber))
            f.write('\n')
            
            ## grid
            f.write('%% Grid \n')
            f.write('\n')
            f.write('vardx\t= {}\n'.format(self.vardx).expandtabs(tabnumber))
            f.write('posdwn\t={}\n'.format(self.posdwn).expandtabs(tabnumber))
            f.write('nx\t= {}\n'.format(self.nx).expandtabs(tabnumber))
            f.write('ny\t= {}\n'.format(self.ny).expandtabs(tabnumber))
            f.write('xori\t= {}\n'.format(self.xori).expandtabs(tabnumber))
            f.write('yori\t= {}\n'.format(self.yori).expandtabs(tabnumber))
            f.write('alfa\t= {}\n'.format(self.alfa).expandtabs(tabnumber)) 
            f.write('xfile\t= x.grd\n'.expandtabs(tabnumber))
            if not self.ygr is None:
                f.write('yfile\t= y.grd\n'.expandtabs(tabnumber))
            f.write('depfile\t= bed.dep\n'.expandtabs(tabnumber))
            f.write('thetamin\t= {}\n'.format(self.thetamin).expandtabs(tabnumber))
            f.write('thetamax\t= {}\n'.format(self.thetamax).expandtabs(tabnumber))
            f.write('thetanaut\t= {}\n'.format(self.thetanaut))
            f.write('dtheta\t= {}\n'.format(self.dtheta).expandtabs(tabnumber))
            f.write('dtheta_s\t= {}\n'.format(self.dtheta).expandtabs(tabnumber))
            f.write('\n')
            
            ## write input vars
            for par_category in self.input_par:
                ## skip category starting with _
                if par_category[0]=='_':
                    continue
                
                ## write meta
                f.write('%% {} \n'.format(par_category))
                f.write('\n')
                for par in self.input_par[par_category]:
                    f.write('{}\t= {}\n'.format(par,self.input_par[par_category][par]).expandtabs(tabnumber))
                f.write('\n')
            ## write output variables
            if '_Output' in self.input_par:
                f.write('%% Output variables \n')
                f.write('\n')
                for par in self.input_par['_Output']:
                    dummy = self.input_par['_Output'][par]
                    f.write('{}\t= {}\n'.format(par,len(dummy)).expandtabs(tabnumber))
                    assert type(dummy)==list, 'expected a list for {}'.format(par)
                    for item in dummy:
                        f.write('{}\n'.format(item))
                    f.write('\n')
    
        ## write grid x
        with open(os.path.join(path,'x.grd'),'w') as f:
            for ii in range(self.ny+1):
                for jj in range(self.nx+1):
                    f.write('{} '.format(self.xgr[ii,jj]))
                f.write('\n')

        if not self.ygr is None:
            ## write grid y
            with open(os.path.join(path,'y.grd'),'w') as f:
                for ii in range(self.ny+1):
                    for jj in range(self.nx+1):
                        f.write('{} '.format(self.ygr[ii,jj]))
                    f.write('\n')

       ## write dep
        with open(os.path.join(path,'bed.dep'),'w') as f:
            for ii in range(self.ny+1):
                for jj in range(self.nx+1):
                    f.write('{} '.format(self.zgr[ii,jj]))
                f.write('\n')             
                
        ## write ne-layer
        if self.struct != None:
            with open(os.path.join(path,'ne_bed.dep'),'w') as f:
                for ii in range(self.ny+1):
                    for jj in range(self.nx+1):
                        f.write('{} '.format(self.nebed[ii,jj]))
                    f.write('\n')   
                
        if self.friction_layer != None:
            with open(os.path.join(path,'friction.dep'),'w') as f:
                for ii in range(self.ny+1):
                    for jj in range(self.nx+1):
                        f.write('{} '.format(self.friction[ii,jj]))
                    f.write('\n')  
                    
        if self.wavefriction_layer != None:
            with open(os.path.join(path,'wavefriction.dep'),'w') as f:
                for ii in range(self.ny+1):
                    for jj in range(self.nx+1):
                        f.write('{} '.format(self.wavefriction[ii,jj]))
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
        elif self.wbctype=='paramtric':
            print('wbctype=paramtric cannot be plotted')
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
            plt.plot(np.squeeze(self.xgr),np.squeeze(self.zgr)*self.posdwn)
            plt.ylabel('z')
            plt.subplot(2,1,2)
            plt.plot(np.squeeze(self.xgr)[1:],np.diff(np.squeeze(self.xgr)))
            plt.xlabel('x')
            plt.ylabel('dx')
        else:
            plt.subplot(2,2,1)
            plt.pcolor(self.xgr,self.ygr,self.zgr*self.posdwn)
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
            plt.pcolor(X_world+self.xori,Y_world+self.yori,self.zgr*self.posdwn)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
            plt.colorbar()
            plt.title('World coordinates - output')

        plt.suptitle(self.fname)

        if self.struct == 1:
            fig2 = plt.figure()
        
            plt.pcolor(self.xgr,self.ygr,self.nebed)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
            plt.title('ne_bed.dep (positive)')
            plt.axis('scaled')
            plt.grid('on')
                
        if self.friction_layer == 1:
            fig3 = plt.figure()
        
            plt.pcolor(self.xgr,self.ygr,self.friction)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
            plt.title('friction.dep (positive)')#+self.input_par['Flow parameters']['bedfriction'])
            plt.axis('scaled')
            plt.grid('on')
                
        if self.wavefriction_layer == 1:
            fig4 = plt.figure()
        
            plt.pcolor(self.xgr,self.ygr,self.wavefriction)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
            plt.title('wavefriction.dep (positive) - fw')
            plt.axis('scaled')
            plt.grid('on')

        if save_path!=None:
            fig1.savefig(os.path.join(save_path,'domain.png'),dpi=250)
            if self.struct == 1:
                fig2.savefig(os.path.join(save_path,'ne_bed.png'),dpi=250)
            if self.friction_layer == 1:
                fig3.savefig(os.path.join(save_path,'friction.png'),dpi=250)
            if self.wavefriction_layer == 1:
                fig4.savefig(os.path.join(save_path,'wavefriction.png'),dpi=250)