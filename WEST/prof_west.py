#%%
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.io import loadmat 

#WEST 
import imas_west
import pywest 
import pywed as pw
from pywed import tsbase
from pywed import tsmat
from matlabtools import Struct

#DBS
from DBS.io.interface import get_angle_pol
from DBS.io.interface import DataInterface

#%%
class profs_WEST():
    
    def __init__(self, shot = None):
        if shot is None:
            shot = int(input('Shot Number: '))
        self.shot = shot
        self.machine = 'west'
        self.t_ignitron = tsmat(shot, 'IGNITRON|1')[0]
    
    def DIFDOP(self, tstart, tend, isweep = None):
        if isweep is None:
            isweep = 5
        time, anglepol = get_angle_pol(self.machine, self.shot, tstart, tend, return_val='array')
        dataI = DataInterface(self.shot, isweep, channelval = 1, machine='west')
        params = dataI.params
        t_DIFDOP = params.TDIFDOP
        t_acq = t_DIFDOP + params.t0seq  
        header_txt = 'DIFDOP parameters'
        difdop = dict({'header':header_txt, 'shot':self.shot, 'time':time, 'anglepol':anglepol, 'window':[tstart, tend], 'parameters':params, 't_DIFDOP':t_DIFDOP, 't_acq':t_acq})
        self.difdop = difdop 
        return difdop
         
    def get_tt(self):
        shot, t_ignitron = self.shot, self.t_ignitron
        summary = imas_west.get(shot, 'summary')
        ece = imas_west.get(shot, 'ece')
        time = summary.time - t_ignitron
        ip = summary.global_quantities.ip.value #total plasma current [A]
        p_ohm = summary.global_quantities.power_ohm.value #Ohmic power [W]
        p_rad = summary.global_quantities.power_radiated.value #total radiated power [W]
        n_e = summary.line_average.n_e.value #[m^-3]
        p_ic = summary.heating_current_drive.power_ic.value #total IC Power [W]
        p_lh = summary.heating_current_drive.power_lh.value #total LH Power [W]
        time_ece = ece.time - t_ignitron
        T_e = ece.channel[0].t_e.data #electron temperature [eV]
        header_txt = 'Summary: time traces'
        time_traces = dict({'header':header_txt, 'shot':shot, 'time':time, 'ip':ip, 'p_ohm':p_ohm, 'p_rad':p_rad, \
            'n_e':n_e,'p_ic':p_ic, 'p_lh':p_lh, 'time_ece':time_ece, 'T_e':T_e})
        self.time_traces = time_traces 
        return time_traces
    
    def plot_tt(self, isweep = None, colors = None, ax = None, save = False):
        tt, shot = self.time_traces, self.shot
        if ax is None:
            fig, ax = plt.subplots(figsize = (15, 4))
        tstart, tend = min(tt['time']), max(tt['time'])
        difdop = self.DIFDOP(tstart, tend)

        #set plots
        ax.set_title(f'WEST shot #{shot}', fontsize = 18)
        ax.grid(c = 'silver', ls ='--', lw = 0.5)
        ax.set_xlabel('time [s]')
        
        #start plotting
        ax.plot(difdop['time'], difdop['anglepol'], 'silver', lw = 0.5, label = 'inclinometer')
        ax.plot(tt['time'], -tt['ip']/2*1e-5, 'magenta', label = r'$I_p$ [200 kA]')
        ax.plot(tt['time'], tt['n_e']*1e-19, 'b', label = r'$\overline{n_e}$ [$10^{19}$ $m^{-3}$]')
        ax.plot(tt['time'], tt['p_ohm']*1e-6, 'g', label = r'$P_{ohm}$ [MW]') 
        if np.any(tt['p_lh'] != 0):
            ax.plot(tt['time'], tt['p_lh']*1e-6, 'r', label = r'$P_{lh}$ [MW]')
        if np.any(tt['p_ic'] != 0):
            ax.plot(tt['time'], tt['p_ic']*1e-6, 'r', label = r'$P_{ic}$ [MW]') 
        if np.any(tt['T_e'] != 0):
            ax.plot(tt['time_ece'], tt['T_e']*1e-3, 'dodgerblue', label = r'$T_{e0}$ [eV]') 
        
        if isweep is not None:
            t_acq = difdop['t_acq']
            for sweep, col in zip(isweep, colors):
                ax.plot(t_acq[sweep-1]*np.ones(100), np.linspace(-0.1, max(difdop['anglepol']), 100), '--', c = col,label = 'it = %d' %sweep)

        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_xlabel('time [s]')
        ax.set_xlim(tstart - 0.5, tend + 1)
        
        if save == True:
            save_fig = gnr_path_timetraces.format(shot)
            plt.savefig(save_fig)
        plt.show()
           
        
        
# %%
if __name__ == '__main__':
    prof = profs_WEST(58333)
    tt = prof.get_tt()
    prof.plot_tt()
    prof.plot_tt(isweep = [5, 12, 24], colors = ['r', 'g', 'b'])
# %%

# %%

# %%
