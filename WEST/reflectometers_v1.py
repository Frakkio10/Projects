#%%
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd

from pywed import tsmat, tsbase
import statistics_west as sw
import equimap
import imas_west  # convenient package for WEST IMAS data
from pywest import shot_to_campaign
from DBS.io.utils import suppress_output
from scipy.io import loadmat 

from matlabtools import Struct

#Code Description:
#flag = 0 means getting the data from imas_west
#flag = 1 means getting the data from file 


class Profile(Struct):
    """Abstract class for profiles."""
    
    def __init__(self, name):
        self.name = name
    
    # @abstractmethod
    def plot(self, ax):
        pass
    
    def select_twindow(self, twindow, time, *quantities, axis=1):
        """Return time array and quantities within twindow. Assumes time is along the <axis>."""

        cond = (time > twindow[0]) & (time < twindow[1])
        indices = np.arange(len(time))
        ind = indices[cond]
    
        return [time[ind]] + [np.take(q, ind, axis=axis) for q in quantities]


class DREFRAP(Profile):
    """DREFRAP reflectometer density profile."""
    
    
    '''initialization function to extract the raw data either from IMAS or from a file:
            -flag = 0 (set by default): extract the data from imas 
            -flag = 1: extract data from a file given by the reflectometry team for example
        then if twindow is specified it select the data in that time window    
    '''
    
    def __init__(self, shot, flag = 0, ne_max = None, drho = 0.0, R_max = 0.0,  twindow=None, verbose=False):
        super().__init__('DREFRAP')
        self.machine = 'WEST'
        self.rho_bord = 1.3
        
        #Obtaining the density data from IMAS
        if flag == 0:
            if verbose:
                data = imas_west.get(shot, 'reflectometer_profile')
            else:
                with suppress_output():
                    data = imas_west.get(shot, 'reflectometer_profile')
            if len(data.time)==0:
                warnings.warn(f'No DREFRAP data found on IMAS for shot {shot}.')
            # time, density profiles and associated radial locations:
            t_ignitron = sw.read_signal_shot(shot, 't_ignitron',shot_to_campaign(shot)).values[0]
            time = data.time  - t_ignitron
            ne = data.channel[0].n_e.data
            R = data.channel[0].position.r 
            Z = data.channel[0].position.z
            Phi = data.channel[0].position.phi
            
            
        #to select the indeces of the raw data based on the time interval of interest
        if twindow is not None:
            time, ne, R, Z, Phi = self.select_twindow(twindow, time, ne, R, Z, Phi, axis=1)
        
        self.shot      = shot
        self.data      = data
        self.t_ignitron= t_ignitron
        self.time      = time
        self.twindow   = twindow
        self.ne        = ne
        self.R         = R
        self.Z         = Z
        self.Phi       = Phi
        self.flag      = flag
        self.drho      = drho
        self.ne_max    = ne_max
        self.R_max     = R_max 
        
    def get_ne(self, twindow=None,  averaged=False, verbose=True):
        
        if self.ne_max is None:
            self.fig, self.ax = plt.subplots()
            self.ax.plot(self.R, self.ne*1e-19, 'silver', lw = 0.5)
            self.ax.set_xlabel(r'R [m]')
            self.ax.set_ylabel(r'n_e [$10^{19}$ $m^{-3}$]')
            plt.show()
            

            #asking the user the maximum density profile to select which data can be used for interpolation 
            self.set_max = input('Do you want to select a maximum density and maximum radius? (y/n)')
            if self.set_max == 'y':
                ne_max = float(input('Maximum Density to select (in 1e19) : '))
                ne_max = ne_max * 1e19
                self.ne_max = ne_max
                R_max = float(input('Profiles have to reach at least this R value (m): '))
                self.R_max = R_max
            elif self.set_max == 'n':
                self.ne_max = 0
                self.R_max = 0
                
                
        else:
            self.ne_max = self.ne_max * 1e19
            self.set_max = 'y'

        
        if twindow is not None:
           self = DREFRAP(shot=self.shot, flag = self.flag, drho = self.drho, ne_max = self.ne_max, R_max = self.R_max, twindow=twindow)

        if not averaged:
            
            #select the data that satisfied the choice of ne_max and R_max 
            if self.set_max == 'y':
                cond = (np.max(self.ne, axis = 0) > self.ne_max) & ( np.max(self.R, axis = 0) > self.R_max)
                ind = np.where(cond)[0]
                ne = self.ne[:, ind]
                R = self.R[:, ind]
                self.ne = ne
                self.R = R
                
            # interpolation from cylindrical to rho coordinates:
            t, R, Phi,Z = self.get_flattened_data(self.time, self.R, self.Phi,self.Z)

            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm') 
            else:
                rho_psi = equimap.get(self.shot, t.mean() + self.t_ignitron, R, Phi, Z, 'rho_pol_norm') 
            
            rho_psi = rho_psi.reshape(*self.R.shape, order='F')
            
            t, R, ne = self.time, self.R, self.ne
            
            data = dict({'shot':self.shot, 't':t, 'rho_psi':rho_psi, 'ne':ne, 'R':R, 'Z':Z,\
                'drho':self.drho, 'ne_max':self.ne_max, 'R_max':self.R_max})
            
            self.data = data

        else: #averaged
            
            if self.set_max == 'y':
                cond = (np.max(self.ne, axis = 0) > self.ne_max) & ( np.max(self.R, axis = 0) > self.R_max)
                ind = np.where(cond)[0]
                self.ne_sel = self.ne[:, ind]
                self.R_sel = self.R[:, ind]
            else:
                self.R_sel = self.R
                self.ne_sel = self.ne
            
            t, Phi, Z = self.get_flattened_data(self.time, self.Phi,self.Z)
            
            R_int, ne_interp = self.interpolate(nmesh=100)
            
            Phi = Phi.mean() * np.ones_like(R_int)
            Z = Z.mean() * np.ones_like(R_int)
            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R_int,Phi,Z,'rho_pol_norm')
            else:
                rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R_int,Phi,Z,'rho_pol_norm')
            
            ne = ne_interp.mean(axis=1)
            self.ne_int = ne
            ne_std = ne_interp.std(axis = 1)            
            self.rho_tmp = np.flip(rho_psi)
            ne_tmp = ne


            if self.rho_tmp[-1] < self.rho_bord:
                rho_add = np.linspace(self.rho_tmp[-1] + (self.rho_bord - self.rho_tmp[-1])/10, self.rho_bord, 100)
                rho_end = self.rho_tmp[-1]
                ne_add = np.linspace(ne_tmp[0]/10 , 0, 100)
                rho_psi = np.append(self.rho_tmp,rho_add)
                ne_tmp = np.append(ne_add, ne_tmp)
                rho_psi = np.flip(rho_psi)
                ne = ne_tmp
            
            header_txt = f'Reflectometry profile averaged over {len(t)} profiles' if averaged else ''
            if twindow is None:
                twindow = [t.min(), t.max()]

            data = dict({'header':header_txt, 'shot':self.shot, 'twindow':twindow, 't':t, 'rho_psi':rho_psi + self.drho, 'ne':ne, 'ne_std':ne_std, 'R':R_int, 'Z':Z, \
                'drho':self.drho, 'ne_max':self.ne_max, 'R_max':self.R_max, 'ne_int':self.ne_int})
            
            self.data = data
            
        return data
    
    def plot_ne(self, twindow=None, xcoord='rho_psi', averaged=False, ax=None, axis_labels=True, 
                shot_annotation=True, errorbars=False, errorband=False, xshift=0, unit_factor=11, **kwargs):

        verbose = kwargs.pop('verbose', False)
        data = self.data
        
        rho_psi = data['rho_psi']
        R = data['R']
        ne = data['ne']
        if averaged:
            ne_std = data['ne_std']

        if ax is None:
            fig, ax = plt.subplots()
        
        if xcoord == 'rho_psi':
            x = rho_psi + xshift
            xlabel = r'$\rho_\psi$'
        elif xcoord == 'R':
            x = R + xshift
            xlabel = r'$R$ [m]'
        else:
            raise ValueError('xcoord must be either "rho_psi" or "R"')
        
        fmt = kwargs.get('fmt', '.')
        kwargs.pop('fmt', None)

        if averaged:
            if errorbars:
                ax.errorbar(np.flip(self.rho_tmp), data['ne_int'] * unit_factor, yerr=ne_std * unit_factor, fmt=fmt, **kwargs)
            if errorband:
                ax.fill_between(np.flip(self.rho_tmp), (data['ne_int'] - ne_std) * unit_factor, (data['ne_int'] + ne_std) * unit_factor, alpha=0.3)
            
            ax.plot(np.flip(self.rho_tmp), data['ne_int'] * unit_factor, fmt, **kwargs)
        else:
            ax.plot(x, ne * unit_factor, fmt, **kwargs)

        if axis_labels:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'$n$ [m$^{-3}$]')

        if shot_annotation:
            # from DBS.src.visualize import annotate_shot
            # annotate_shot(self.shot, ax=ax, machine=self.machine)
            pass

    plot = plot_ne # shorthand

    def get_flattened_data(self, *quantities):
        return [q.reshape(-1,1, order='F').flatten() for q in quantities]


    def get_indices_within_nsigma(X, nsigma=1, center=np.median):
        return ( X < (center(X) + nsigma * np.std(X)) ) & ( (X > (center(X) - nsigma * np.std(X))) )
    

    def select_on_common_R_mesh(R, ne, nmesh, ascending_R=True):

        R_min = np.max(np.min(R, axis=0)) 
        R_max = np.min(np.max(R, axis=0)) 

        if ascending_R:
            R_mesh = np.linspace(R_max, R_min, nmesh)

            if R[0,0] > R[-1,0]:
                R = np.flip(R, axis=0)
                ne = np.flip(ne, axis=0)
        else:
            R_mesh = np.linspace(R_min, R_max, nmesh)
                
        return R_mesh, R, ne

    def interpolate(self, nmesh=100):
        """ interpolates the profile on a regular grid """
        
        R_mesh, R, ne = DREFRAP.select_on_common_R_mesh(self.R_sel, self.ne_sel, nmesh=nmesh, ascending_R=True)
        
        from scipy.interpolate import pchip_interpolate #, interp1d
        from scipy.interpolate import CubicSpline

        def pchip_interpolation(X, Y, X_new):
            """ wrapper for pchip interpolation along first axis of the 2d arrays X and Y onto the 1d array X_new """
            Y_new = np.zeros((X_new.shape[0], Y.shape[1]))
            for i in range(X.shape[1]):
                Y_new[:,i] = pchip_interpolate(X[:,i], Y[:,i], X_new, axis=0)
                
            return Y_new
        
        ne_new = pchip_interpolation(R, ne, R_mesh)

        return R_mesh, ne_new
# %%
if __name__ == '__main__':
    nprof = DREFRAP(57558, flag = 0,  drho = 0.0, twindow=[7.6, 7.8], verbose = False)    
    data = nprof.get_ne(averaged = True, verbose = False)
    
    #nprof2.plot_ne(errorband = False, errorbars = False, averaged=True, verbose = True)
    fig, ax = plt.subplots(1, 2, figsize = (10, 4))
    ax[0].plot(nprof.R, nprof.ne*1e-19, 'silver', lw = 0.5)
    ax[0].plot(data['R'], data['ne_int']*1e-19, 'b')
    ax[0].set_xlabel('R [m]')
    ax[0].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
    ax[1].plot(data['rho_psi'], data['ne']*1e-19, 'b')
    ax[1].set_xlabel(r'$\rho_\psi$')
    ax[1].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')

# %%
