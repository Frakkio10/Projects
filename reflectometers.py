#%%
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd

import statistics_west as sw
import equimap
import imas_west  # convenient package for WEST IMAS data
from pywest import shot_to_campaign
from DBS.io.utils import suppress_output


# from abc import ABC, abstractmethod
from matlabtools import Struct
#%%

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
        indices = np.arange(len(time))[cond]

        return [time[indices]] + [np.take(q, indices, axis=axis) for q in quantities]


class DREFRAP(Profile):
    """DREFRAP reflectometer density profile."""
    
    
    def __init__(self, shot, twindow=None, data=None, verbose=False):
        super().__init__('DREFRAP')
        self.machine = 'WEST'

        if data is None:

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

        if twindow is not None:
            time, ne, R, Z, Phi = self.select_twindow(twindow, time, ne, R, Z, Phi, axis=1)
            # t, ne, R, Z, Phi = self.get_twindow(twindow)

        self.shot      = shot
        self.data      = data
        self.t_ignitron= t_ignitron
        self.time      = time
        self.ne        = ne
        self.R         = R
        self.Z         = Z
        self.Phi       = Phi


    def get_ne(self, twindow=None,  averaged=False, Rshift=0, verbose=True):


        if twindow is not None:
           self = DREFRAP(shot=self.shot, twindow=twindow, data=self.data)

        if not averaged:

            # interpolation from cylindrical to rho coordinates:
            t, R, Phi,Z = self.get_flattened_data(self.time, self.R, self.Phi,self.Z)


            # apply horizontal for adjustment:
            R += Rshift

            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm')
            else:
                rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm')
                                
            rho_psi = rho_psi.reshape(*self.R.shape, order='F')

            
            t, R, ne = self.time, self.R, self.ne

            data = dict({'shot':self.shot, 't':t, 'rho_psi':rho_psi, 'ne':ne, 'R':R, 'Z':Z})
        
        else:
            t, Phi,Z = self.get_flattened_data(self.time, self.Phi,self.Z)

            R, ne_interp = self.interpolate(nmesh=50)

            # apply horizontal shift for adjustment:
            R += Rshift
            
            Phi = Phi.mean() * np.ones_like(R)
            Z = Z.mean() * np.ones_like(R)

            if not verbose:
                with suppress_output():
                    rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm')
            else:
                rho_psi = equimap.get(self.shot,t.mean() + self.t_ignitron,R,Phi,Z,'rho_pol_norm')
            # ax.errorbar(R_mesh, ne_mean, yerr=ne_std, fmt='o', color='black', markersize=3, capsize=3, capthick=1, zorder=10)
            
            ne = ne_interp.mean(axis=1)
            ne_std = ne_interp.std(axis=1)
            
            header_txt = f'Reflectometry profile averaged over {len(t)} profiles' if averaged else ''
            if twindow is None:
                twindow = [t.min(), t.max()]
                
            data = dict({'header':header_txt, 'shot':self.shot, 'twindow':twindow, 't':t, 'rho_psi':rho_psi, 'ne':ne, 'ne_std':ne_std, 'R':R, 'Z':Z, 'Rshift':Rshift, })

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
        # ax.plot(x, ne, fmt, **kwargs)

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


    def get_indices_within_nsigma(X, nsigma=2, center=np.median):
        return ( X < (center(X) + nsigma * np.std(X)) ) & ( (X > (center(X) - nsigma * np.std(X))) )
    

    def select_on_common_R_mesh(R, ne, nmesh=50, ascending_R=True):

        cond_R_min = DREFRAP.get_indices_within_nsigma( np.min(R, axis=0), nsigma=2)
        cond_R_max = DREFRAP.get_indices_within_nsigma( np.max(R, axis=0), nsigma=2) # this one is less crucial since ne ~ 0 at the edge
        cond = cond_R_min & cond_R_max

        R_sel = R[:,cond]
        ne_sel = ne[:,cond]

        # print(cond_R_min.sum(), cond_R_max.sum(), cond.sum(), R.shape)
        print("discarded fraction:  {:.2f}".format(1 - cond.sum() / R.shape[1]))

        R_min = np.max(np.min(R, axis=0)[cond])
        R_max = np.min(np.max(R, axis=0)[cond])

        if ascending_R:
            R_mesh = np.linspace(R_max, R_min, nmesh)
            if R_sel[0,0] > R_sel[-1,0]:
                R_sel = np.flip(R_sel, axis=0)
                ne_sel = np.flip(ne_sel, axis=0)
        else:
            R_mesh = np.linspace(R_min, R_max, nmesh)


        return R_mesh, R_sel, ne_sel

    def interpolate(self, nmesh=50):
        """ interpolates the profile on a regular grid """

        R_mesh, R, ne = DREFRAP.select_on_common_R_mesh(self.R, self.ne, nmesh=nmesh, ascending_R=True)

        from scipy.interpolate import pchip_interpolate #, interp1d


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

    # %matplotlib widget
    nprof = DREFRAP(57954, twindow=[5,6])
    data = nprof.get_ne(verbose=False, averaged=True)
    
    nprof.plot(averaged=True, errorbars=True, errorband=True, unit_factor=1e-19)

# %%
if __name__ == '__main__':

    # %matplotlib widget
    prof = Profile('test')
    prof


    # %%
    # in case no DREFRAP available:
    prof = DREFRAP(58155, twindow=[6,7], verbose=True)
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

