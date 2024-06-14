#%%
#general python 
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.io import loadmat 

#WEST 
import imas_west
import pywest 
import pywed as pw
from pywest import polview

#DBS 
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam
from DBS.beamtracing.src.west.io import retreive_west_density_data
from DBS.io.interface import get_angle_pol
from DBS.io.interface import DataInterface
from DBS.io.interface import extract_DBS_files

#fro
from WEST.custom_beamtracing import custom_beam3d
from WEST.reflectometers_v1 import DREFRAP

gnr_path_fdop = "/Home/difdop/DBSdata/processed/fDop_estimation/west/{}_FO/_ch{}_isweep{}.mat"
gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}_inc.mat"

# %% functions

machine = 'west'

def get_tt(shot):
    t_ignitron = pw.tsbase(shot, 'RIGNITRON')[0][0,0]
    summ = imas_west.get(shot, 'summary')
    ece = imas_west.get(shot, 'ece')
    
    summary_df = pd.DataFrame()
    summary_df['time'] = summ.time - t_ignitron 
    summary_df['ip'] = summ.global_quantities.ip.value 
    summary_df['p_ohm'] = summ.global_quantities.power_ohm.value 
    summary_df['n_e'] = summ.line_average.n_e.value 
    summary_df['p_ic'] = summ.heating_current_drive.power_ic.value 
    summary_df['p_lh'] = summ.heating_current_drive.power_lh.value 

    ece_df = pd.DataFrame()
    ece_df['time_ece'] = ece.time - t_ignitron
    ece_df['T_e'] = ece.t_e_central.data

    return summary_df, ece_df

def plot_tt(shot, sweeps, save = False):
    global machine
    
    summary, ece = get_tt(shot)
    tstart, tend = min(summary.time), max(summary.time)

    t, anglepol = get_angle_pol(machine, shot, tstart, tend,t_averaged = False, use_inclinometer = False)
    #t, anglepol_inc = get_angle_pol(machine, shot, tstart, tend,t_averaged = False, use_inclinometer = True)

    dataI = DataInterface(shot, isweep = 5, channelval = 1, machine='west')
    params = dataI.params
    t_DIFDOP = params.TDIFDOP
    t_acq = t_DIFDOP + params.t0seq  
    
    fig, ax = plt.subplots( figsize = (15, 4))
    ax.set_title('SHOT #%d' %shot,  size = 18)
    #ax.plot(t , anglepol_inc, 'silver', lw = 0.8, label = r'inclinometre [°]')
    ax.plot(t , anglepol, 'skyblue', lw = 0.8, label = r'verin [°]')

    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    ax.grid(c = 'silver', ls ='--', lw = 0.5)
    ax.set_xlabel('time [s]')
    ax.plot(summary['time'], -summary['ip']/2*1e-5, 'magenta', label = r'$I_p$ [200kA]') #200kA
    ax.plot(summary['time'], summary['n_e']*1e-19, 'b', label = r'$\overline{n_e}$ [$10^{19}$ $m^{-3}$]') #1e19
    ax.plot(summary['time'], summary['p_ic']*1e-6, 'red', label = r'$P_{ICRH}$ [MW]') #MW
    ax.plot(summary['time'], summary['p_lh']*1e-6, 'r', label = r'$P_{LH}$ [MW]') #MW
    #ax.plot(summary['time'], summary.p_ohm*1e-6, 'coral', label = r'$P_{Ohm}$ [MW]') #MW
    ax.plot(ece['time_ece'], ece['T_e']*1e-3, 'dodgerblue', label = r'$T_{e0}$ [eV]') #eV
    
    colors = ['r', 'g', 'b']
    for sweep, col in zip(sweeps, colors):
        ax.plot(t_acq[sweep-1]*np.ones(100), np.linspace(-0.1, max(anglepol), 100), '--', c = col,label = 'it = %d' %sweep)

    ax.grid(c = 'silver', ls ='--', lw = 0.5)
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    ax.set_xlabel('time [s]')
    ax.set_xlim(min(summary.time) - 0.5,max(summary.time) + 0.5)
    #ax.set_xlim(4.3, 5)
    
    if save == True:
        save_fig1 = gnr_path_timetraces1.format(shot)
        save_fig2 = gnr_path_timetraces2.format(shot)
        plt.savefig(save_fig1)
        plt.savefig(save_fig2)
    plt.show()

def DBS_prof(shot, isweep, channelval, modex):
    global gnr_path_fdop, gnr_path_btr
        
    filename_fdop = gnr_path_fdop.format(shot, channelval, isweep)
    filename_btr = gnr_path_btr.format(shot, channelval, modex, isweep)

    try:
        fdop = loadmat(filename_fdop)['outp'][0]
    except FileNotFoundError:
        print(f'File {filename_fdop} not found')
        return 
    try:
        beamtr = loadmat(filename_btr)['outp'][0]
    except FileNotFoundError:
        print(f'File {filename_btr} not found')
        return        
        
    profile = pd.DataFrame()
    
    profile['ifreq'] = np.arange(1, fdop['time'][0][0].size + 1)
    profile['f0'] = fdop['freqGHz'][0][0]
    profile['t'] = fdop['time'][0][0]
    profile['rho_psi'] = beamtr['rho'][0][0]
    profile['k_perp'] = beamtr['k_perp'][0][0]*2e2
    profile['fDop'] = fdop['fDop'][0][0]*1e-3
    profile['dfDop'] = fdop['dfDop'][0][0]*1e-3
    profile['fDop_min'] = fdop['fDop_min'][0][0]*1e-3
    profile['fDop_max'] = fdop['fDop_max'][0][0]*1e-3
    profile['v_perp'] = 2 * np.pi * profile.fDop / profile.k_perp
    profile['dv_perp'] = 2 * np.pi * profile.dfDop / profile.k_perp
    profile['v_perp_low'] = 2 * np.pi * profile.fDop_min / profile.k_perp
    profile['v_perp_up'] = 2 * np.pi * profile.fDop_max / profile.k_perp
    profile['validated'] = fdop['validated'][0][0]
    profile['shot'] = shot*np.ones(fdop['time'][0][0].size)
    profile['isweep'] = isweep*np.ones(fdop['time'][0][0].size)

    return profile

def plot_prof(profile, color, font, lbl, ax = None, errorbar = False):

    validated = profile.validated
    if ax == None:
        fig, ax = plt.subplots()
    
    k_perp = np.mean(profile.k_perp[validated == 1])*1e-2
    dk_perp = np.std(profile.k_perp[validated == 1])*1e-2

    if errorbar == True:
        combined = zip(profile.rho_psi[validated == 1], profile.v_perp[validated == 1], profile.dv_perp[validated == 1])
        sort_combined = sorted(combined, key= lambda pair: pair[0])
        rho_psi, v_perp, dv_perp = zip(*sort_combined)
        ax.errorbar(rho_psi, -1*np.array(v_perp), yerr = -1*np.array(dv_perp), c = color, fmt = font, label = r'$k_{\perp}$ = %.1f $\pm$ %.1f $cm^{-1}$' %(k_perp, dk_perp))
    
    else:
        combined = zip(profile.rho_psi[validated == 1], profile.v_perp[validated == 1])
        sort_combined = sorted(combined, key= lambda pair: pair[0])
        rho_psi, v_perp= zip(*sort_combined)
        ax.plot(rho_psi, -1*np.array(v_perp), font, c = color, label = r'$k_{\perp}$ = %.1f $\pm$ %.1f $cm^{-1}$' %(k_perp, dk_perp))

def conc(profile_O, profile_X):
    #Omode
    validated = profile_O.validated
    combined = zip(profile_O.rho_psi[validated == 1], profile_O.v_perp[validated == 1], profile_O.dv_perp[validated == 1])
    sort_combined = sorted(combined, key= lambda pair: pair[0])
    rho_psiO, v_perpO, dv_perpO = zip(*sort_combined)
    
    #Xmode
    validated = profile_X.validated
    combined = zip(profile_X.rho_psi[validated == 1], profile_X.v_perp[validated == 1], profile_X.dv_perp[validated == 1])
    sort_combined = sorted(combined, key= lambda pair: pair[0])
    rho_psiX, v_perpX, dv_perpX = zip(*sort_combined)
    rho_psi = np.concatenate([rho_psiO, rho_psiX])
    v_perp = np.concatenate([v_perpO, v_perpX])
    dv_perp = np.concatenate([dv_perpO,dv_perpX])

    return rho_psi, v_perp, dv_perp
# %% example for the plots 

#time traces
shot = 58108
_ = extract_DBS_files(shot, machine=machine,  extract_all=True)
isweeps = [9, 12, 19]
plot_tt(shot, isweeps, save = False)

_#%%
#equilibrium 
shots = [57558]
colors = ['r']
time = [5]
figure, ax = plt.subplots()
for shot, col, t  in zip(shots, colors, time):
    polview(shot, t, ax=ax, colors=col)
ax.legend(shots)
#%%
#density profile 

fig, ax = plt.subplots()
nprof = DREFRAP(shot, flag = 0, Rshift = -0.0278, twindow=[4.8, 5.0], verbose = True)
data = nprof.get_ne(verbose=False, averaged=False)
nprof.plot(ax = ax, averaged=True, errorbars=True, errorband=True, unit_factor=1e-19)
ax.set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]', fontsize = 12)
ax.set_xlabel(r'$\rho_{\psi}$', fontsize = 12)

# %% Custom beamtracing
t_start, dt = 10.8, 0.2
machine, shot, twindow, modex, channelval = 'west', 60269, [t_start, t_start + dt], 1, 2
sweep, Rshift = 40, 0
flag = 1 #flag = 0 for data from imas_west
outp, interface = custom_beam3d(machine, shot, sweep, Rshift, twindow, modex, channelval, flag, verbose = True)


# %%
shot, channelval, modex = 58333, 2, 1
isweeps = [8, 24]
colors3 = ['r', 'b']

fig, ax = plt.subplots(figsize = (10, 4))
ax.set_title(f'SHOT #{shot}', fontsize = 20)
for isweep, col in zip(isweeps, colors3):
    profile = DBS_prof(shot, isweep, channelval, modex)
    #profile.v_per = -profile.v_perp
    profile_O = DBS_prof(shot, isweep, 1, 0)
    profile_X = DBS_prof(shot, isweep, 2, 1)
    rho_psi, v_perp, dv_perp = conc(profile_O, profile_X)
    #plot_prof(profile, col, 'o', ax = ax, lbl = None, errorbar = True) 
    #ax.plot(rho_psi, v_perp, 'o', c = col, label = f'{isweep}')
    #ax.errorbar(rho_psi, v_perp, yerr = dv_perp, c = col, fmt = 's', alpha = 0.8, label = f'it {isweep}')
    ax.plot(rho_psi, v_perp, 'o-', c = col, alpha = 0.8, label = f'it {isweep}' )

ax.legend(loc = 'center left', bbox_to_anchor = (0.8, 0.2))
ax.set_xlabel(r'$\rho_{\psi}$', fontsize = 16)
ax.set_ylabel(r'$v_{\perp}$ [km/s]', fontsize = 16)
#ax.plot(np.linspace(-0.96,1.02, 100), np.zeros(100), 'k')
#ax.plot(np.ones(100), np.linspace(-10, 5, 100), 'k')
#ax.set_xlim(1.0150, 1.0175)
#ax.set_ylim(-0.5, 0.5)
ax.set_xlim(1m 1)
ax.grid(c = 'silver', ls ='--', lw = 0.5)
#ax.plot(np.linspace(1.02, 1.03, 100), np.zeros(100), 'k')
#ax.plot(1.0278*np.ones(100), np.linspace(-1, 1, 100))

#%%
shot, channelval, modex = 60269, 2, 1
isweeps = [12, 29, 40]
colors3 = ['r', 'g', 'b']

fig, ax = plt.subplots(figsize = (10, 4))
ax.set_title(f'SHOT #{shot}', fontsize = 20)
for isweep, col in zip(isweeps, colors3):
    profile = DBS_prof(shot, isweep, channelval, modex)
    ax.plot(profile.f0, profile.fDop, '.', c = col, label = int(profile.isweep[0]))
    
ax.legend(loc = 2)
#%%
shot, channelval, modex = 58333, 1, 0
isweeps = [4, 17, 24]
colors3 = ['r', 'g', 'b']
t_start = [11.8, 19.6, 23.8]
dt = 0.6

fig, ax = plt.subplots(figsize = (6, 4))
ax.set_title(f'SHOT #{shot}', fontsize = 20)

for isweep, col, ti in zip(isweeps, colors3, t_start):
    profile = DBS_prof(shot, isweep, channelval, modex )

    t, anglepol = get_angle_pol(machine, shot, ti, ti + dt, return_val='array')
    theta = np.mean(anglepol)
    validated = profile.validated
    combined = zip(profile.rho_psi[validated == 1], profile.k_perp[validated == 1])
    sort_combined = sorted(combined, key= lambda pair: pair[0])
    rho_psi, k_perp= zip(*sort_combined)
    if isweep == 17:
        p = np.polyfit(rho_psi[:-1], np.array(k_perp)[:-1]*1e-2, 1)
        ax.plot(np.linspace(0.65, 0.92), p[0]*np.linspace(0.65, 0.92) + p[1], '--', lw = 1.5,  c = col, label = 'y = %.1f x + %.1f' %(p[0], p[1]))
    else:
        p = np.polyfit(rho_psi[:-2], np.array(k_perp)[:-2]*1e-2, 1)
        ax.plot(np.linspace(0.65, 0.92), p[0]*np.linspace(0.65, 0.92)+ p[1], '--', lw = 1.5, c = col, label = 'y = %.1fx + %.1f'%(p[0], p[1]))

    ax.plot(rho_psi, np.array(k_perp)*1e-2, 'o-', markersize = 5, c = col,   label = r'$\theta$ = %.1f' %theta)

ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
#ax.set_xlim(0.95, 1.03)
#ax.set_xlim(0.6, 0.95)
ax.set_xlabel(r'$\rho_{\psi}$', fontsize = 12)
ax.set_ylabel(r'$k_{\perp}$', fontsize = 12)
# %%

shot, isweep = 58333, 24
profile_X = DBS_prof(shot, isweep, 2, 1)
profile_O = DBS_prof(shot, isweep, 1, 0)


# %%
dt = 0.6
ti = 17.8
t, anglepol = get_angle_pol(machine, 58333, ti, ti + dt, return_val='array')
theta = np.mean(anglepol)
print(theta)
#%%
DBS_prof(57558, 5, 2, 1)
#%%