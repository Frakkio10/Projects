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
from DBS.beamtracing import DBSbeam
import DBS
from scipy.io import loadmat

#%% Functions
gnr_path_fdop = "/Home/difdop/DBSdata/processed/fDop_estimation/west/{}_FO/_ch{}_isweep{}.mat"
gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}.mat"
gnr_path_fdop_LV = '/Home/difdop/DIFDOP_LV/processed/FitSpectre/{}{}_{}'
gnr_path_savefig = '/Home/FO278650/Bureau/Analysis/fig_analised/FO{}_{}.pdf'


#time_trace
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

#Laure's data
def compare_spectra(shot, isweep, ch, start_fr, end_fr):
    global gnr_path_fdop, gnr_path_btr, gnr_path_fdop_LV

    if ch == 1:
        mode = 0
        band = 'V'
    else:
        mode = 1
        band = 'W'

    #load LV data
    filename_fdop = gnr_path_fdop_LV.format(shot, band, isweep)
    data =  loadmat(filename_fdop)
    
    #load beamtracing data
    filename_btr = gnr_path_btr.format(shot, ch, mode, isweep)
    beam_trc = loadmat(filename_btr)['outp'][0]

    LV_df = pd.DataFrame()
    LV_df['val'] = data['table'][isweep-1,start_fr:end_fr] == 1
    LV_df['f_prob'] = data['F_pal'][0, start_fr:end_fr]
    LV_df['fD_G'] = data['fcG_fit'][isweep-1, start_fr:end_fr]
    LV_df['fD_L'] = data['fcL_fit'][isweep-1, start_fr:end_fr]
    LV_df['fD_T'] = data['fcT_fit'][isweep-1, start_fr:end_fr]
    LV_df['k_perp'] = beam_trc['k_perp'][0][0, start_fr:end_fr]*2e2
    LV_df['rho'] = beam_trc['rho'][0][0, start_fr:end_fr]
    LV_df['v_perp'] = 2*np.pi*LV_df.fD_L/LV_df.k_perp
    
    return LV_df

#DBS load data
def DBS_df(shot, ch, isweep):
    global gnr_path_fdop, gnr_path_btr
    
    if ch == 1:
        mode = 0
    else:
        mode = 1
        
    filename_fdop = gnr_path_fdop.format(shot, ch, isweep)
    filename_btr = gnr_path_btr.format(shot, ch, mode, isweep)
    
    fdop_est = loadmat(filename_fdop)['outp'][0]
    beam_trc = loadmat(filename_btr)['outp'][0]
    
    validated = fdop_est['validated'][0][0]
    
    DBS_spect_df = pd.DataFrame()
    DBS_btr_df = pd.DataFrame()
    
    val = fdop_est['validated'][0][0]
    DBS_spect_df['fDop'] =  fdop_est['fDop'][0][0, val == 1]*1e-3
    DBS_spect_df['dfDop'] =  fdop_est['dfDop'][0][0, val == 1]*1e-3
    DBS_spect_df['pbfreq'] =  fdop_est['freqGHz'][0][0, val == 1]
    DBS_btr_df['k_perp'] = beam_trc['k_perp'][0][0, val == 1]*2e2
    DBS_btr_df['rho'] = beam_trc['rho'][0][0, val == 1]
    #DBS_btr_df['v_perp'] = 2*np.pi*DBS_spect_df.fDop/DBS_btr_df.k_perp
    
    return DBS_spect_df, DBS_btr_df



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

LV_spectraV = compare_spectra(shot, isweep, 1, 0, 20)
LV_spectraW = compare_spectra(shot, isweep, 2, 0, 20)

print('hello')



# %% plots: v_perp, time_trace, angle
fig, ax = plt.subplots(3,1, figsize = (15, 10))

ax[0].set_title('SHOT #%d' %shot,  size = 18)
ax[0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO O-mode')
ax[0].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, 'bx', markersize = 8, label = 'FO X-mode')
ax[0].plot(LV_spectraV.rho[LV_spectraV.val], LV_spectraV.v_perp[LV_spectraV.val]*1e-3, 'r.', markersize = 5, label = 'LV O-mode')
ax[0].plot(LV_spectraW.rho[LV_spectraW.val], LV_spectraW.v_perp[LV_spectraW.val]*1e-3, 'rx', markersize = 5, label = 'LV X-mode')
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

save_fig = gnr_path_savefig.format(shot, 'complete')
fig.savefig(save_fig)
plt.show()

# %% pltos: f_dop, v_perp 
fig, ax = plt.subplots(2,2, figsize = (20,8))
fig.suptitle('SHOT #%d' %shot, size = 18)
ax[0,0].set_title('Iteration %d V-band' %isweep)
ax[0,0].plot(prof_V.f0[prof_V.validated == 1], prof_V.fDop[prof_V.validated == 1]*1e-3, '.', markersize = 8, label = 'Francesco')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_G[LV_spectraV.val]*1e-3, 'x', markersize = 3, label = 'Gaussian')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_L[LV_spectraV.val]*1e-3, 's', markersize = 3, label = 'Lorentian')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_T[LV_spectraV.val]*1e-3, '^', markersize = 3, label = 'Taylor')
ax[0,0].legend()
ax[0,0].set_xlabel('Probing Frequency [GHz]')
ax[0,0].set_ylabel('Doppler Frequency [kHz]')
ax[0,0].grid(c = 'silver', ls ='--', lw = 0.5)


ax[1,0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO')
ax[1,0].plot(LV_spectraV.rho[LV_spectraV.val], LV_spectraV.v_perp[LV_spectraV.val]*1e-3, 'rx', markersize = 5,  label = 'LV')
ax[1,0].set_xlabel(r'$\rho$')
ax[1,0].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[0,1].set_title('Iteration %d W-band' %isweep)
ax[0,1].plot(prof_W.f0[prof_W.validated == 1], prof_W.fDop[prof_W.validated == 1]*1e-3, '.', markersize = 8, label = 'Francesco')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_G[LV_spectraW.val]*1e-3, 'x', markersize = 3, label = 'Gaussian')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_L[LV_spectraW.val]*1e-3, 's', markersize = 3, label = 'Lorentian')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_T[LV_spectraW.val]*1e-3, '^', markersize = 3, label = 'Taylor')
ax[0,1].legend()
ax[0,1].set_xlabel('Probing Frequency [GHz]')
ax[0,1].set_ylabel('Doppler Frequency [kHz]')
ax[0,1].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,1].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO')
ax[1,1].plot(LV_spectraW.rho[LV_spectraW.val], LV_spectraW.v_perp[LV_spectraW.val]*1e-3, 'rx', markersize = 5,  label = 'LV')
ax[1,1].legend()
ax[1,1].set_xlabel(r'$\rho$')
ax[1,1].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,1].grid(c = 'silver', ls ='--', lw = 0.5)

print(shot)
save_fig = gnr_path_savefig.format(shot, 'fdop_vperp')
fig.savefig(save_fig)
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

LV_spectraV = compare_spectra(shot, isweep, 1, 0, 20)
LV_spectraW = compare_spectra(shot, isweep, 2, 0, 20)
# %% plots: v_perp, time_trace, angle
fig, ax = plt.subplots(3,1, figsize = (15, 10))

fig.suptitle('SHOT #%d' %shot,  size = 18)
ax[0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO O-mode')
ax[0].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, 'bx', markersize = 8, label = 'FO X-mode')
ax[0].plot(LV_spectraV.rho[LV_spectraV.val], LV_spectraV.v_perp[LV_spectraV.val]*1e-3, 'r.', markersize = 5, label = 'LV O-mode')
ax[0].plot(LV_spectraW.rho[LV_spectraW.val], LV_spectraW.v_perp[LV_spectraW.val]*1e-3, 'rx', markersize = 5, label = 'LV X-mode')
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
ax[2].set_xlabel('time[s]')
ax[2].legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
ax[2].grid(c = 'silver', ls ='--', lw = 0.5)

print(shot)
save_fig = gnr_path_savefig.format(shot, 'complete')
fig.savefig(save_fig)
plt.show()
# %% plots: f_dop, v_perp 
fig, ax = plt.subplots(2,2, figsize = (20,8))
fig.suptitle('SHOT #%d' %shot, size = 18)
ax[0,0].set_title('Iteration %d V-band' %isweep)
ax[0,0].plot(prof_V.f0[prof_V.validated == 1], prof_V.fDop[prof_V.validated == 1]*1e-3, '.', markersize = 8, label = 'Francesco')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_G[LV_spectraV.val]*1e-3, 'x', markersize = 3, label = 'Gaussian')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_L[LV_spectraV.val]*1e-3, 's', markersize = 3, label = 'Lorentian')
ax[0,0].plot(LV_spectraV.f_prob[LV_spectraV.val], LV_spectraV.fD_T[LV_spectraV.val]*1e-3, '^', markersize = 3, label = 'Taylor')
ax[0,0].legend()
ax[0,0].set_xlabel('Probing Frequency [GHz]')
ax[0,0].set_ylabel('Doppler Frequency [kHz]')
ax[0,0].grid(c = 'silver', ls ='--', lw = 0.5)


ax[1,0].plot(prof_V.rho_psi[prof_V.validated == 1], prof_V.v_perp[prof_V.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO')
ax[1,0].plot(LV_spectraV.rho[LV_spectraV.val], LV_spectraV.v_perp[LV_spectraV.val]*1e-3, 'rx', markersize = 5,  label = 'LV')
ax[1,0].set_xlabel(r'$\rho$')
ax[1,0].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[0,1].set_title('Iteration %d W-band' %isweep)
ax[0,1].plot(prof_W.f0[prof_W.validated == 1], prof_W.fDop[prof_W.validated == 1]*1e-3, '.', markersize = 8, label = 'Francesco')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_G[LV_spectraW.val]*1e-3, 'x', markersize = 3, label = 'Gaussian')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_L[LV_spectraW.val]*1e-3, 's', markersize = 3, label = 'Lorentian')
ax[0,1].plot(LV_spectraW.f_prob[LV_spectraW.val], LV_spectraW.fD_T[LV_spectraW.val]*1e-3, '^', markersize = 3, label = 'Taylor')
ax[0,1].legend()
ax[0,1].set_xlabel('Probing Frequency [GHz]')
ax[0,1].set_ylabel('Doppler Frequency [kHz]')
ax[0,1].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,1].plot(prof_W.rho_psi[prof_W.validated == 1], prof_W.v_perp[prof_W.validated == 1]*1e-3, 'b.', markersize = 8, label = 'FO')
ax[1,1].plot(LV_spectraW.rho[LV_spectraW.val], LV_spectraW.v_perp[LV_spectraW.val]*1e-3, 'rx', markersize = 5,  label = 'LV')
ax[1,1].legend()
ax[1,1].set_xlabel(r'$\rho$')
ax[1,1].set_ylabel(r'$v_{\perp}$ [km/s]')
ax[1,1].grid(c = 'silver', ls ='--', lw = 0.5)

print(shot)
save_fig = gnr_path_savefig.format(shot, 'fdop_vperp')
fig.savefig(save_fig)
plt.show()
# %%

#(machine, shot, time=None, xmode=False, isweep=None, channelval=None, twindow=None, ifreqs='all', verbose=True, plot=True, ask_before_run=False, load_if_existing=False, ax=None, **kwargs)
output, beam3d_interface = DBSbeam('west',58333, isweep = 8, channelval = 1, ifreqs = 'all')
[beam3d_interface.launcher[i].anglepol for i in range(len(beam3d_interface.launcher))]
# %%
t, anglepol = get_angle_pol('west', 58333, 0, 100, return_val='array')

plt.plot(t, anglepol)
# %%
beam3d_interface.launcher
# %%

