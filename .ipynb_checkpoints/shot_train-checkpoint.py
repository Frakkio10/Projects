#%%
import numpy as np
import matplotlib.pyplot as plt
import imas_west
import scipy.io as sc
import pywest
import pywed as pw
import pandas as pd
from DBS.io.interface import get_angle_pol
from pywest import polview
from DBS.analysis import DBS_Profile
import DBS
# %%
def time_trace(shot):
    t_ignitron = pw.tsbase(shot, 'RIGNITRON')[0][0,0]
    summ = imas_west.get(shot, 'summary')
    ece = imas_west.get(shot, 'ece')
    
    summary_df = pd.DataFrame()
    summary_df['time'] = summ.time - t_ignitron 
    summary_df['ip'] = summ.global_quantities.ip.value 
    summary_df['n_e'] = summ.line_average.n_e.value 
    summary_df['p_ic'] = summ.heating_current_drive.power_ic.value 
    summary_df['p_lh'] = summ.heating_current_drive.power_lh.value 

    ece_df = pd.DataFrame()
    ece_df['time_ece'] = ece.time - t_ignitron
    ece_df['T_e'] = ece.t_e_central.data
    
    return summary_df, ece_df

#%% shot 57558
shot = 57558
isweep = 11

#profiles
summary, ece = time_trace(shot)

#data analysed
machine, xmode, ch = 'west', False, 1
args = (machine, shot, isweep, xmode, ch)
prof_V = DBS_Profile(*args)

machine, xmode, ch = 'west', True, 2
args = (machine, shot, isweep, xmode, ch)
prof_W = DBS_Profile(*args)

#angle
tstart, tend = min(summary.time), max(summary.time)
t, anglepol = get_angle_pol(machine, shot, tstart, tend, return_val='array')
machine, shot, tstart, tend = 'west', shot, 0.0, 10.0
t, anglepol = get_angle_pol(machine, shot, tstart, tend, return_val='array')
# %% plots: v_perp, time_trace, angle
fig, ax = plt.subplots(3,1, figsize = (10, 10))

ax[0].set_title('SHOT #%d' %shot,  size = 18)
ax[0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, '.', markersize = 8, label = 'O-mode')
ax[0].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, 'rx', markersize = 8, label = 'X-mode')
ax[0].set_xlabel(r'$\rho$')
ax[0].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[0].legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
ax[0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1].plot(summary['time'], -summary['ip']/2*1e-5, 'magenta', label = r'$I_p$ [200kA]') #200kA
ax[1].plot(summary['time'], summary['n_e']*1e-19, 'lime', label = r'$\overline{n_e}$ [$10^{19}$ $m^{-3}$]') #1e19
ax[1].plot(summary['time'], summary['p_ic']*1e-6, 'red', label = r'$P_{ICRH}$ [MW]') #MW
ax[1].plot(summary['time'], summary['p_lh']*1e-6, 'orange', label = r'$P_{LH}$ [MW]') #MW
ax[1].plot(ece['time_ece'], ece['T_e']*1e-3, 'dodgerblue', label = r'$T_{e0}$ [eV]') #eV
ax[1].grid(c = 'silver', ls ='--', lw = 0.5)
ax[1].legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
ax[1].set_xlabel('time [s]')
ax[1].set_xlim(min(summary.time) - 0.5,max(summary.time) + 0.5)

ax[2].plot(t, anglepol, label = r'$\alpha$ [rad]')
ax[2].legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
ax[2].grid(c = 'silver', ls ='--', lw = 0.5)
ax[2].set_xlabel('time [s]')

plt.show()
# %% pltos: f_dop, v_perp 
fig, ax = plt.subplots(2,2, figsize = (20,8))
fig.suptitle('SHOT #%d' %shot, size = 18)
ax[0,0].set_title('Iteration %d V-band' %isweep)
ax[0,0].plot(prof_V.f0[prof_V.validated == 1], prof_V.fDop[prof_V.validated == 1]*1e-3, '.', markersize = 8)
ax[0,0].set_xlabel('Probing Frequency [GHz]')
ax[0,0].set_ylabel('Doppler Frequency [kHz]')
ax[0,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, '.', markersize = 8)
ax[1,0].set_xlabel(r'$\rho$')
ax[1,0].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,0].grid(c = 'silver', ls ='--', lw = 0.5)

machine, shot, isweep, xmode, ch = 'west', 57558, 11, True, 2
args = (machine, shot, isweep, xmode, ch)
prof = DBS_Profile(*args)

ax[0,1].set_title('Iteration %d W-band' %isweep)
ax[0,1].plot(prof_W.f0[prof_W.validated == 1], prof_W.fDop[prof_W.validated == 1]*1e-3, '.', markersize = 8)
ax[0,1].set_xlabel('Probing Frequency [GHz]')
ax[0,1].set_ylabel('Doppler Frequency [kHz]')
ax[0,1].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,1].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, '.', markersize = 8)
ax[1,1].set_xlabel(r'$\rho$')
ax[1,1].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,1].grid(c = 'silver', ls ='--', lw = 0.5)

plt.show()
# %% shot 58333
shot = 58333
isweep = 8

#profiles
summary, ece = time_trace(shot)

#data analysed
machine, xmode, ch = 'west', False, 1
args = (machine, shot, isweep, xmode, ch)
prof_V = DBS_Profile(*args)

machine, xmode, ch = 'west', True, 2
args = (machine, shot, isweep, xmode, ch)
prof_W = DBS_Profile(*args)

#angle
tstart, tend = min(summary.time), max(summary.time)
t, anglepol = get_angle_pol(machine, shot, tstart, tend, return_val='array')
# %%
