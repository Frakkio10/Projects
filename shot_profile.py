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
from DBS.beamtracing.src.west.io import retreive_west_density_data

import DBS
from scipy.io import loadmat
from DBS.io.interface import DataInterface
from DBS.io.interface import extract_DBS_files

gnr_path_timetraces1 = '/Home/FO278650/Bureau/Analysis/timetraces/FO{}.pdf'
gnr_path_timetraces2 = '/Home/FO278650/Bureau/Analysis/timetraces/FO{}.jpg'
gnr_path_probfreq = '/Home/FO278650/Bureau/Analysis/probingF/FO{}.jpg'

gnr_path_fdop = "/Home/difdop/DBSdata/processed/fDop_estimation/west/{}_FO/_ch{}_isweep{}.mat"
gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}.mat"



# %%
machine = 'west'
def plot_tt(shot, summary, ece, tstart, tend, sweeps, save = False):
    global machine, gnr_path_timetraces
    
    t, anglepol = get_angle_pol(machine, shot, tstart, tend, return_val='array')
    dataI = DataInterface(shot, isweep = 5, channelval = 1, machine='west')
    params = dataI.params
    t_DIFDOP = params.TDIFDOP
    t_acq = t_DIFDOP + params.t0seq  
    
    fig, ax = plt.subplots( figsize = (15, 4))
    ax.set_title('SHOT #%d' %shot,  size = 18)
    ax.plot(t , anglepol, 'silver', lw = 0.8, label = r'$\alpha$ [rad]')
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    ax.grid(c = 'silver', ls ='--', lw = 0.5)
    ax.set_xlabel('time [s]')
    ax.plot(summary['time'], -summary['ip']/2*1e-5, 'magenta', label = r'$I_p$ [200kA]') #200kA
    ax.plot(summary['time'], summary['n_e']*1e-19, 'lime', label = r'$\overline{n_e}$ [$10^{19}$ $m^{-3}$]') #1e19
    ax.plot(summary['time'], summary['p_ic']*1e-6, 'red', label = r'$P_{ICRH}$ [MW]') #MW
    ax.plot(summary['time'], summary['p_lh']*1e-6, 'orange', label = r'$P_{LH}$ [MW]') #MW
    ax.plot(summary['time'], summary.p_ohm*1e-6, 'teal', label = r'$P_{Ohm}$ [MW]') #MW
    #ax.plot(ece['time_ece'], ece['T_e']*1e-3, 'dodgerblue', label = r'$T_{e0}$ [eV]') #eV
    
    for sweep in sweeps:
        ax.plot(t_acq[sweep-1]*np.ones(100), np.linspace(min(anglepol), max(anglepol), 100), '--', label = 'it = %d' %sweep)

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
    
def plot_fprob(shot, isweep, save = False):
    global machine, gnr_path_probfreq
    
    fig, ax = plt.subplots(1,2, figsize = (15,4))
    fig.suptitle('#SHOT: %d' %shot)
    #Vband
    ax[0].set_title('V-band O-mode')
    dataI = DataInterface(shot, isweep, channelval=1, machine = 'west')
    params = dataI.params
    ax[0].plot(np.linspace(1, len(params.F), len(params.F)), params.F, '.' )
    #Wband
    ax[1].set_title('W-band X-mode')
    dataI = DataInterface(shot, isweep, channelval=2, machine='west')
    params = dataI.params
    ax[1].plot(np.linspace(1, len(params.F), len(params.F)), params.F, '.' )
    ax[0].set_ylabel('Probing Frequency [GHz]')
    ax[1].set_ylabel('Probing Frequency [GHz]')

    if save == True:
        save_fig1 = gnr_path_probfreq.format(shot)
        plt.savefig(save_fig1)
        
    plt.show()
    
def time_trace(shot):
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

def sel(k, dk, sweeps):
    global profiles
    
    for sweep in sweeps:
        profW = profiles[sweep]
        val = np.zeros(len(profW))
        k_perp = np.array(profW.k_perp)
        
        for i in range(0, len(profW)):
            if (k_perp[i] >= k - dk) & (k_perp[i] <= k + dk):
                val[i] = 1
        
        if sweep == sweeps[0]:
            selected = profW[val == 1]
        else:
            selected = pd.concat([selected, profW[val == 1]])
            
    return selected 

#%%
def DBS_df(shot, ch, isweep, sort = False):
    global gnr_path_fdop, gnr_path_btr
    
    if ch == 1:
        mode = 0
    else:
        mode = 1
        
    filename_fdop = gnr_path_fdop.format(shot, ch, isweep)
    filename_btr = gnr_path_btr.format(shot, ch, mode, isweep)
    try:
        fdop_est = loadmat(filename_fdop)['outp'][0]
        val = fdop_est['validated'][0][0]

    beam_trc = loadmat(filename_btr)['outp'][0]
    
    except FileNotFoundError:
        print(f'File {filename_fdop} not found')
        return 
        
    DBS_spect_df = pd.DataFrame()
    
    val = fdop_est['validated'][0][0]
    
    if sort == True:
        #sort fprob and fdop
        combined = zip(fdop_est['freqGHz'][0][0, val == 1], fdop_est['fDop'][0][0, val == 1]*1e-3, fdop_est['dfDop'][0][0, val == 1]*1e-3)
        sort_combined = sorted(combined, key= lambda pair: pair[0])
        f_prob, f_Dop, df_Dop = zip(*sort_combined)
        DBS_spect_df['fDop'] =  f_Dop
        DBS_spect_df['dfDop'] =  df_Dop
        DBS_spect_df['pbfreq'] =  f_prob
        
        return DBS_spect_df

    
    DBS_spect_df = pd.DataFrame()
    DBS_btr_df = pd.DataFrame()
    
    val = fdop_est['validated'][0][0]
    DBS_spect_df['fDop'] =  fdop_est['fDop'][0][0, val == 1]*1e-3
    DBS_spect_df['dfDop'] =  fdop_est['dfDop'][0][0, val == 1]*1e-3
    DBS_spect_df['pbfreq'] =  fdop_est['freqGHz'][0][0, val == 1]
    DBS_btr_df['k_perp'] = beam_trc['k_perp'][0][0, val == 1]*2e2
    DBS_btr_df['rho_psi'] = beam_trc['rho'][0][0, val == 1]
    DBS_btr_df['v_perp'] = 2*np.pi*DBS_spect_df.fDop/DBS_btr_df.k_perp
    
    return DBS_spect_df, DBS_btr_df
    #return DBS_spect_df
#%%
#ne_reflec=out_reflec.channel{1}.n_e.data;
#phi_reflec=out_reflec.channel{1}.position.phi.data;r_reflec=out_reflec.channel{1}.position.r.data;z_reflec=out_reflec.channel{1}.position.z.data;
out_reflec = imas_west.get(60270, 'reflectometer_profile')
ne_reflec = out_reflec.channel[0].n_e.data
pos = out_reflec.channel[0].position
r_refl, z_refl, phi_refl = pos.r, pos.z, pos.phi

#%% example
shot = 58333
sweeps = [4, 18, 24]
_ = extract_DBS_files(shot, machine=machine,  extract_all=True)

summary, ece = time_trace(shot)
tstart, tend = min(summary.time), max(summary.time)
plot_tt(shot, summary, ece, tstart, tend, sweeps, save = False )
#plot_fprob(shot, 20, save = False)


#%%
df_60270_21_1 = DBS_df(shot, 1, 21)
df_60270_21_2 = DBS_df(shot, 2, 21)
df_60270_28_1 = DBS_df(shot, 1, 28)
df_60270_28_2 = DBS_df(shot, 2, 28)


fig, ax = plt.subplots()

ax.plot(df_60270_21_1.pbfreq, df_60270_21_1.fDop, 'o', c = 'b')
ax.plot(df_60270_21_2.pbfreq, df_60270_21_2.fDop, 'x', c = 'b')
ax.plot(df_60270_28_1.pbfreq, df_60270_28_1.fDop, 'o', c = 'r')
ax.plot(df_60270_28_2.pbfreq, df_60270_28_2.fDop, 'x', c = 'r')

#%% 57715 sweep [5, 8, 9, 11] freqs [1:20] ch 2

#shot info 

machine, xmode, ch = 'west', False, 1
sweeps = [8, 17, 25]
#saving the analyzed profiles
profiles = {}
for sweep in sweeps:
    profile_name = f'profile_{sweep}'
    args = (machine, shot, sweep , xmode, ch)
    profile_data = DBS_Profile(*args)[0:20]
    profiles[sweep] = profile_data

#%%#defining the k_perp to analyze
k_arr = np.arange(1, 3, 0.5)*1e3
dk = 0.3*1e3
selected = {}
for k in k_arr:
    selected_name = f'selected_{k*1e-2}'
    selected_data = sel(k, dk, sweeps)
    selected[k*1e-2] = selected_data

#%%
for k in [1000, 2000]:
    val = selected[k*1e-2].validated
    #plt.plot(np.sort(selected[k*1e-2].rho_psi[val == 1]), np.sort(selected[k*1e-2].v_perp[val == 1]*1e-3), '-o', label = k*1e-2)
    plt.plot(selected[k*1e-2].rho_psi[val == 1], selected[k*1e-2].v_perp[val == 1]*1e-3, 'o', label = k*1e-2)

plt.legend()
#%%
def plot_prof(prof, color, font, lbl, ax = None,  err = False, save = False):
    combined = zip(prof.rho_psi[prof.validated == 1], prof.v_perp[prof.validated == 1]*1e-3, prof.dv_perp[prof.validated == 1]*1e-3)
    sort_combined = sorted(combined, key= lambda pair: pair[0])
    rho_psi, v_perp, dv_perp = zip(*sort_combined)
    
    if err == True:
        ax.errorbar(rho_psi, v_perp, yerr = dv_perp, c = color, fmt = font, label = lbl)
    else:
        ax.plot(rho_psi, v_perp, c = color, marker = font, label = lbl)
        
    if save == True:
        plt.savefig()

#%%
machine = 'west'
shot = 58333
isweep = 24
proff = DBS_df(shot, 2, isweep, sort = False)


#%% READ PROFILE
machine = 'west'
shot = 58333
isweep = 4
args = (machine, shot, isweep, False, 1)
proff_Omode4 = DBS_Profile(*args)
args = (machine, shot, isweep, True, 2)
proff_Xmode4 = DBS_Profile(*args)

isweep = 24
args = (machine, shot, isweep, True, 2)
proff_Xmode24 = DBS_Profile(*args)

isweep = 17
args = (machine, shot, isweep, True, 2)
proff_Xmode17 = DBS_Profile(*args)
#%%
fig, ax = plt.subplots(figsize = (10,4))
#plot_prof(proff_Omode4, 'dodgerblue', 'o', 'Omode 4', ax, err = False)
plot_prof(proff_Xmode4, 'dodgerblue', 'x', 'Xmode 4', ax, err = False)
plot_prof(proff_Xmode24, 'r', 's', 'Xmode 24', ax, err = False)
plot_prof(proff_Xmode17, 'g', 'o', 'Xmode 17', ax, err = False)


ax.legend()
ax.grid(c = 'silver', ls ='--', lw = 0.5)
#sort



#%% PLOT PROFILE 

fig, ax = plt.subplots(figsize = (10, 4))
#ax.errorbar(proff_Omode.rho_psi[proff_Omode.validated == 1], proff_Omode.v_perp[proff_Omode.validated == 1], yerr = proff_Omode.dv_perp[proff_Omode.validated == 1], fmt = 'o-', label = 'O-mode')
#ax.errorbar(proff_Xmode.rho_psi[proff_Xmode.validated == 1], proff_Xmode.v_perp[proff_Xmode.validated == 1], yerr = proff_Xmode.dv_perp[proff_Xmode.validated == 1], fmt = 'x-', label = 'X-mode')
#ax.plot(proff_Omode.rho_psi[proff_Omode.validated == 1], proff_Omode.v_perp[proff_Omode.validated == 1], 'bo', label = 'O-mode')
#ax.plot(proff_Xmode.rho_psi[proff_Xmode.validated == 1], proff_Xmode.v_perp[proff_Xmode.validated == 1], 'bx', label = 'X-mode')
#ax.errorbar(rho_psi_O, v_perp_O, yerr = dv_perp_O, fmt = 'o')
#ax.set_ylim(-12e3, 2e3)
#%% plot information 

t_start = [14.2, 19.6, 24.4]

fig, ax = plt.subplots(figsize = (6, 4))
for i in range(0, len(sweeps)):
    t, anglepol = get_angle_pol(machine, shot, t_start[i], t_start[i] + 0.5, return_val='array')
    ax.plot(profiles[sweeps[i]].rho_psi, profiles[sweeps[i]].k_perp*1e-2, '.', label = r'it = %d $\alpha$ = %.1f rad' %(sweeps[i], np.mean(anglepol)))


ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$k_\perp$ [$cm^{-1}$]')
ax.legend()

plt.show()

#%%

from pywest import polview
#%%

shots = [58333, 57558]
time = [10, 7.8]
figure, ax = plt.subplots()
for shot, col, t  in zip(shots, ['r', 'b','g'], time):
    polview(shot, t, ax=ax, colors=col)
ax.legend(shots)


# %% to go back on this 
#dataI = DataInterface(shot, isweep='False', channelval=1, machine='west')

shot = 58333
dataI = DataInterface(shot, isweep = 5, channelval=2, machine='west')
params = dataI.params
t_difdop = params.TDIFDOP
t_acq = t_difdop + params.t0seq         
t, anglepol = get_angle_pol(machine, shot, tstart, tend, return_val='array')
fig, ax = plt.subplots(figsize = (15,4))
ax.plot(t  , anglepol, 'silver', label = r'$\alpha$ [rad]')

sweeps = [19]
#    ax.plot(t_acq[sweeps[i]-1]*np.ones(100), np.linspace(min(anglepol), max(anglepol), 100), '--', label = 'it = %d' %sweeps[i])

#for i in range(0, len(t_acq), 10):
#   ax.plot(t_acq[i]*np.ones(100), np.linspace(min(anglepol), max(anglepol), 100), 'gray', '-.')
for i in range(0, len(sweeps)):
    ax.plot(t_acq[sweeps[i]-1]*np.ones(100), np.linspace(min(anglepol), max(anglepol), 100), '--', label = 'it = %d' %sweeps[i])
ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))

ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.set_xlabel('time [s]')
ax.set_xlim(0, 13)

#%%
time = t + t_difdop
idx = (time > 8.4) & (time < 8.6)
angle = np.mean(anglepol[idx])

#%%
shot = 58333
t_start = 23.8
dt = 0.6
t, anglepol = get_angle_pol(machine, shot, t_start, t_start + dt, return_val='array')

print(np.mean(anglepol))
#plt.plot(t, anglepol)
# %%
shot = 58333
fig, ax = plt.subplots(1,2, figsize = (8,4))

fig.suptitle('#SHOT: %d' %shot)
dataI = DataInterface(shot, isweep = 20, channelval=1, machine='west')
params = dataI.params
ax[0].plot(np.linspace(1, len(params.F), len(params.F)), params.F, '.' )
dataI = DataInterface(shot, isweep = 20, channelval=2, machine='west')
params = dataI.params
ax[1].plot(np.linspace(1, len(params.F), len(params.F)), params.F, '.' )

ax[0].set_ylabel('Probing Frequency [GHz]')


# %%

# %%
from DBS.beamtracing import DensityProf1d
machine, shot = 'west', 57715
time = np.array([4.7, 6.5, 7.1, 8.3])


#%%
for t in time:
    dens_prof = DensityProf1d.from_shot('west', shot, t , channelval = 2)
    plt.plot(dens_prof.rho_psi, dens_prof.ne, label = time)
    
plt.legend()
# %%
dataI = DataInterface(57715, isweep = 8, channelval=2, machine='west')
params = dataI.params
# %%
time = params.t0seq
# %%
dens_prof = DensityProf1d.from_shot('west', 58333, isweep = 8, channelval = 2)

# %%
shot=58333
isweep=8
channelval=2
machine='west'

dataI = DataInterface(shot, isweep, channelval, machine)
params = dataI.params
dtseq = np.diff(dataI.params.t0seq)[0]
t_difdop = params.TDIFDOP

t0 = dataI.params.t0seq[isweep-1]
time = np.round(t0 + dtseq / 2, 2)
twindow = [np.round(t0, 2), np.round(t0 + dtseq, 2)]
# %%
shot = 57715
twindow = [6.2, 6.8]
density = retreive_west_density_data(shot, twindow, override=False, verbose=True)
# %%

from pywest import reflectometers 
from DBS.beamtracing.src.west import io 

#dens = retreive_west_density_data(60269, [5.2, 5.4])
DREFRAP = reflectometers.DREFRAP

nprof = DREFRAP(58333, [14.2, 14.8])

dens = retreive_west_density_data(58333, [14.2, 14.8])

# %%
dens
# %%
plt.plot(dens.rho_psi, dens.ne)
# %%
#%%
shot = 58333
interf = imas_west.get(shot, 'interferometer')
refle = imas_west.get(shot, 'reflectometer_profile')
bolo = imas_west.get(shot, 'bolometer')
ece = imas_west.get(shot, 'ece')
equib = imas_west.get(shot, 'equilibrium')
icrf = imas_west.get(shot, 'ic_antennas')
lhrf = imas_west.get(shot, 'lh_antennas')
# %%

# %%
P_LH1 = lhrf.antenna[0].power_launched.data
P_LH2 = lhrf.antenna[1].power_launched.data

time = lhrf.antenna[0].power_launched.time



plt.plot(time - 30, P_LH1*1e-6)
plt.plot(time - 30, P_LH2*1e-6)
plt.plot(time - 30, (P_LH1 + P_LH2)*1e-6)

plt.ylim(0,10)
plt.xlim(20, 20.3)
# %%
dens = refle.channel[0].n_e.data
rad = refle.channel[0].position.r

for time in [0, 25, 50, 75, 100]:
    plt.plot(rad[:, time], dens[:,time], label = time)
    
plt.legend()
plt.show()
# %% comparison of shot 58333 and 60269
shot = 58333
#machine, xmode, channel = 'west', [False, True], [1,2]
machine, xmode, channels = 'west', True,  [1,2]
sweeps = [8, 17, 25]
#saving the analyzed profiles
profiles = {}

for sweep in sweeps:
    for ch in channels:
        profile_name = f'profile_{ch}_{sweep}'
        profile_data = DBS_df(shot, ch, sweep, sort = True)
        profiles[sweep, ch] = profile_data
# %%
color = ['r', 'g', 'b']
fig, ax = plt.subplots(figsize = (10, 4))
i = 0
for sweep in sweeps:
    if type(profiles[sweep, 1]) is not pd.core.frame.DataFrame:
        continue
    ax.errorbar(profiles[sweep, 1].pbfreq, profiles[sweep, 1].fDop, yerr = profiles[sweep, 1].dfDop, fmt = '-o', c = color[i])
    ax.errorbar(profiles[sweep, 2].pbfreq, profiles[sweep, 2].fDop, yerr = profiles[sweep, 2].dfDop, fmt = '-o', c = color[i])
    i += 1


# %%
shot = 60269
#machine, xmode, channel = 'west', [False, True], [1,2]
machine, xmode, channels = 'west', True,  [1,2]
sweeps = [29, 40]
#saving the analyzed profiles
profiles = {}

for sweep in sweeps:
    for ch in channels:
        profile_name = f'profile_{ch}_{sweep}'
        profile_data = DBS_df(shot, ch, sweep, sort = True)
        profiles[sweep, ch] = profile_data
# %%
color = ['r', 'g', 'b']
fig, ax = plt.subplots(figsize = (10, 4))
i = 0
for sweep in sweeps:
    ax.errorbar(profiles[sweep, 1].pbfreq, profiles[sweep, 1].fDop, yerr = profiles[sweep, 1].dfDop, fmt = '-o', c = color[i])
    ax.errorbar(profiles[sweep, 2].pbfreq, profiles[sweep, 2].fDop, yerr = profiles[sweep, 2].dfDop, fmt = '-x', c = color[i])
    i += 1
    
#%%

# %% Beamtracing 
from DBS.beamtracing.DBSbeam import _DBSbeam 


machine, shot, isweep, xmode, channelval = 'west', 57558, 19, 1, 1

output, beam3d_interface = _DBSbeam(machine, shot, isweep = 11, xmode = 0, channelval = 1, verbose = True, plot = True, ifreqs = 'all')

# %%
from pywest import polview

shots = [55562, 59824, 59785]
figure, ax = plt.subplots()
for shot, col in zip(shots, ['r', 'b','g']):
    polview(shot, time=6.5, ax=ax, colors=col)
ax.legend(shots)