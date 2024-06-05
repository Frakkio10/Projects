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

# from abc import ABC, abstractmethod
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

        #cond = (time > twindow[0]) & (time < twindow[1]) #original SR
        cond = (time > twindow[0]) & (time < twindow[1])
        indices = np.arange(len(time))
        ind = indices[cond]
    

        return [time[ind]] + [np.take(q, ind, axis=axis) for q in quantities]


class DREFRAP(Profile):
    """DREFRAP reflectometer density profile."""
    
    
    def __init__(self, shot, flag, drho = 0.0, ne_max = None, R_max = None,  twindow=None, verbose=False):
        super().__init__('DREFRAP')
        self.machine = 'WEST'
        self.rho_bord = 1.3
        
        #Obtaining the density data 
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
            
        elif flag == 1:
            gnr_path_DREFRAP = '/Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_{}_prof.mat'
            filemat_ne = gnr_path_DREFRAP.format(shot)
            try:
                data = loadmat(filemat_ne)
            except FileNotFoundError:
                print(f'File {filemat_ne} not found')
                return
            t_ignitron = sw.read_signal_shot(shot, 't_ignitron',shot_to_campaign(shot)).values[0]
            time = np.transpose(data['tX'][0])
            #time = time  - t_ignitron
            ne = np.transpose(data['NEX'])
            R = np.transpose(data['RX'])
            a, b = np.shape(ne)
            Phi = 2.4435*np.ones([a,b]) 
            Z = np.zeros([a,b])

        mod_raf = tsmat(shot, 'DREFRAP;BALAYAGE_V;mod_raf_V')
        
        #if mod_raf == 1:
            #function 
            
        if twindow is not None:
        #if twindow is None:
            time, ne, R, Z, Phi = self.select_twindow(twindow, time, ne, R, Z, Phi, axis=1)
            # t, ne, R, Z, Phi = self.get_twindow(twindow)
        
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
            # interpolation from cylindrical to rho coordinates:
            t, R, Phi,Z = self.get_flattened_data(self.time, self.R, self.Phi,self.Z)
            
            # apply horizontal for adjustment:
            #R += self.drho

            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm') 
            else:
                rho_psi = equimap.get(self.shot, t.mean() + self.t_ignitron, R, Phi, Z, 'rho_pol_norm') 
            
            rho_psi = rho_psi.reshape(*self.R.shape, order='F')
            
            #starting the edit 
            if self.set_max == 'y':
                #cond1 = np.max(self.ne, axis = 0) > self.ne_max
                #cond2 = np.max(self.R, axis = 0) > self.R_max
                cond = (np.max(self.ne, axis = 0) > self.ne_max) & ( np.max(self.R, axis = 0) > self.R_max)
                #ind = np.where(cond1 & cond2)[0]
                ind = np.where(cond)[0]
                self.ne = self.ne[:, ind]
                self.R = self.R[:, ind]
            #finish the edit
            #t, R, ne = self.time, self.R, self.ne
            
            data = dict({'shot':self.shot, 't':t, 'rho_psi':rho_psi, 'ne':ne, 'R':R, 'Z':Z,\
                'ne_max':self.ne_max, 'R_max':self.R_max})
        
        else: #averaged
            
            #self.R += self.drho
            #starting the edit 
            if self.set_max == 'y':
                cond = (np.max(self.ne, axis = 0) > self.ne_max) & ( np.max(self.R, axis = 0) > self.R_max)
                ind = np.where(cond)[0]
                self.ne_sel = self.ne[:, ind]
                self.R_sel = self.R[:, ind]
            #finish the edit
            else:
                self.R_sel = self.R
                self.ne_sel = self.ne
            
            t, Phi, Z = self.get_flattened_data(self.time, self.Phi,self.Z)
            
            R_int, ne_interp = self.interpolate(nmesh=100)
            
            # apply horizontal shift for adjustment:
            #R += self.drho
            
            Phi = Phi.mean() * np.ones_like(R_int)
            Z = Z.mean() * np.ones_like(R_int)
            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R_int,Phi,Z,'rho_pol_norm')
            else:
                #rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm')
                rho_psi = equimap.get(self.shot,t.mean(),R_int,Phi,Z,'rho_pol_norm')         
            
            ne = ne_interp.mean(axis=1)
            self.ne_int = ne
            ne_std = ne_interp.std(axis = 1)            
            rho_tmp = np.flip(rho_psi)
            ne_tmp = ne

            if rho_tmp[-1] < self.rho_bord:
                rho_add = np.linspace(rho_tmp[-1] + (self.rho_bord - rho_tmp[-1])/10, self.rho_bord, 100)
                rho_end = rho_tmp[-1]
                ne_add = np.linspace(ne_tmp[0]/10 , 0, 100)
                
                rho_tmp = np.append(rho_tmp,rho_add)
                ne_tmp = np.append( ne_add, ne_tmp)
                rho_psi = np.flip(rho_tmp)
                ne = ne_tmp
            
            
            header_txt = f'Reflectometry profile averaged over {len(t)} profiles' if averaged else ''
            if twindow is None:
                twindow = [t.min(), t.max()]

            data = dict({'header':header_txt, 'shot':self.shot, 'twindow':twindow, 't':t, 'rho_psi':rho_psi + self.drho, 'ne':ne, 'ne_std':ne_std, 'R':R_int, 'Z':Z, \
                'drho':self.drho, 'ne_max':self.ne_max, 'R_max':self.R_max, 'ne_int':self.ne_int})

        return data
    
    def plot_ne(self, twindow=None, xcoord='rho_psi', averaged=False, ax=None, axis_labels=True, 
                shot_annotation=True, errorbars=False, errorband=False, xshift=0, unit_factor=1, **kwargs):

        verbose = kwargs.pop('verbose', True)
        data = self.get_ne(twindow=twindow, averaged=averaged, verbose=verbose)

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
        
        fmt = kwargs.get('fmt', 'o')
        kwargs.pop('fmt', None)

        if averaged:
            if errorbars:
                ax.errorbar(x, ne * unit_factor, yerr=ne_std * unit_factor, fmt=fmt, **kwargs)
            if errorband:
                ax.fill_between(x, (ne - ne_std) * unit_factor, (ne + ne_std) * unit_factor, alpha=0.3)
            
            ax.plot(x, ne * unit_factor, fmt, **kwargs)
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

        #cond_R_min = DREFRAP.get_indices_within_nsigma( np.min(R, axis=0), nsigma=2)
        #cond_R_max = DREFRAP.get_indices_within_nsigma( np.max(R, axis=0), nsigma=2) # this one is less crucial since ne ~ 0 at the edge
        #cond = cond_R_min & cond_R_max

        #R_sel = R[:,cond]
        #ne_sel = ne[:,cond]

        # print(cond_R_min.sum(), cond_R_max.sum(), cond.sum(), R.shape)
        #print("discarded fraction:  {:.2f}".format(1 - cond.sum() / R.shape[1]))

        R_min = np.max(np.min(R, axis=0)) #original from SR
        R_max = np.min(np.max(R, axis=0)) #original from SR

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
    
#%%
if __name__ == '__main__':
    nprof = DREFRAP(57558, flag = 0, drho = 0, twindow=[7.6, 7.8], verbose = False)    
    data = nprof.get_ne(verbose=False, averaged=True)

    
    fig, ax = plt.subplots(1, 2, figsize = (10, 4))
    ax[0].plot(nprof.R , nprof.ne*1e-19, 'silver', lw = 0.5)
    ax[0].plot(data['R'], data['ne_int']*1e-19, 'r')
    ax[0].set_xlabel('R [m]')
    ax[0].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
    ax[1].plot(data['rho_psi'], data['ne']*1e-19, 'r')
    ax[1].set_xlabel(r'$\rho_\psi$')
    ax[1].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
        

    

# %%
if __name__ == '__main__':
    #fig, ax = plt.subplots( figsize = (8,8))
    # %matplotlib widget
    nprof = DREFRAP(60269, flag = 1, twindow=[5.2, 5.4], drho = 0 , verbose = False)
    data0 = nprof.get_ne(verbose=False, averaged=True)
    plt.plot(data0['R'], data0['ne'])
    plt.xlim(2.5, 3.2)
    plt.show()
    nprof.plot(twindow=[5.2, 5.4], averaged=True, errorbars=True, errorband=False, unit_factor=1e-19)    
    

#%%
if __name__ == '__main__':
    #fig, ax = plt.subplots( figsize = (8,8))
    nprof = DREFRAP(60269, 1, 0, twindow=[5.2, 5.4], drho = 0 , verbose = False)
    data0 = nprof.get_ne(verbose=False, averaged=True)
    nprof = DREFRAP(60269, flag = 1, twindow=[5.2, 5.4], drho = -0.037 , verbose = False)
    data1 = nprof.get_ne(verbose=False, averaged=True)
    #nprof.plot(twindow=[8.6, 8.8], averaged=True, errorbars=True, errorband=False, unit_factor=1e-19)
    fig, ax = plt.subplots(2, 1, figsize = (6,8))
    
    ax[0].plot(data0['R'], data0['ne'], lw = 5, label = 'no shift')    
    ax[0].plot(data1['R'], data1['ne'], label = ' shift')    
    ax[0].set_xlabel('R')
    ax[1].plot(data0['rho_psi'], data0['ne'], label = 'no shift')    
    ax[1].plot(data1['rho_psi'], data0['ne'], label = 'shift') 
    ax[1].set_xlabel(r'$\rho_\psi$')

    
    ax[0].legend()   
    ax[1].legend()  
    ax[0].grid()
    ax[1].grid()
    

#%%
#plt.plot(nprof['R'].mean(axis = 1), nprof['ne'].mean(axis = 1))


# %%
if __name__ == '__main__':

    # %matplotlib widget
    prof = Profile('test')
    prof


    # %%
    # in case no DREFRAP available:
    prof = DREFRAP(60270, twindow=[6,7], verbose=True)
    # prof = DREFRAP(54896, twindow=[4,9])
    prof.R.shape, prof.time.shape, prof.ne.shape
    # len(prof.time)==0
    # prof.plot_ne()

    # %%
    fig, axs = plt.subplots(1,2,figsize=(9,3))

    shots = [54896, 54903]
    # twindow = [8,9]
    twindow = [3, 4]
    # twindow = [8, 8.5]

    legend_elements = []

    for i, (shot, config, col) in enumerate(zip(shots, ['USN', 'LSN'], ['blue', 'red'])):

        prof = DREFRAP(shot, twindow=twindow)

        # _ = prof.ne

        for j, (ax, xcoord) in enumerate(zip(axs, ['R', 'rho_psi'])):
            

            # raw profiles:
            kwargs = {
                    #   'twindow': twindow,
                    'xcoord': xcoord,
                    'ax': ax,
                    # 'fmt': '-',
                    #   'markersize': 0.1,
                    'axis_labels': True,
                    'shot_annotation': False,
                    'color':col,
                    }
            prof.plot_ne(averaged=False, alpha=0.05, fmt='.', markersize=0.5, **kwargs)
            prof.plot_ne(averaged=True, fmt='-', lw=1, errorbars=True, capsize=2, capthick=1, **kwargs)
            kwargs.pop('color')
            prof.plot_ne(averaged=True, fmt='-', lw=3, color='black', **kwargs)
            


            # kwargs = {
            #         #   'twindow': twindow,
            #           'xcoord': xcoord,
            #           'averaged': True,
            #           'ax': ax,
            #           'fmt': '-',
            #         #   'markersize': 0.1,
            #           'axis_labels': True,
            #           'shot_annotation': False,
            #           'color':col,
            #           'alpha':0.3
            #           }

            if j==0:
                legend_elements.append(plt.Line2D([0], [0], color=col, lw=2, label=f"{shot} ({config}), t={twindow[0]:.1f}-{twindow[1]:.1f}s"))


    # Create a custom legend
    # legend = plt.legend(handles=legend_elements)

    axs[0].legend(handles=legend_elements, loc='lower left')
    plt.tight_layout()

#%%
shot = 60269
gnr_path_DREFRAP = '/Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_{}_prof.mat'
filemat_ne = gnr_path_DREFRAP.format(shot)
data = loadmat(filemat_ne)


t_ignitron = sw.read_signal_shot(shot, 't_ignitron',shot_to_campaign(shot)).values[0]
time = np.transpose(data['tX'][0])
#time = time  - t_ignitron
ne = np.transpose(data['NEX'])
R = np.transpose(data['RX'])
a, b = np.shape(ne)
Phi = 2.4435*np.ones([a,b]) 
Z = np.zeros([a,b])
# %%
