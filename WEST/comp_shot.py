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
from pywest import polview
from matlabtools import Struct

#DBS
from DBS.io.interface import get_angle_pol
from DBS.io.interface import DataInterface
from DBS.beamtracing.interface import Beam3dInterface
from WEST.summary_WEST import TIMETRACE
from WEST.reflectometers_v1 import DREFRAP



#matlab
import matlab.engine
#eng = matlab.engine.start_matlab()
#eng.quit()
# %%
FO_fdop_gnr = '/Home/difdop/DBSdata/processed/fDop_estimation/west/{}_FO/_ch{}_isweep{}.mat'
FO_btra_gnr = '/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}_{}_dr{}.mat'
#FO_btra_gnr = '/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}_{}.mat'
gnr_path_fdop_LV = '/Home/difdop/DIFDOP_LV/processed/FitSpectre/{}{}_{}'
gnr_path_btra_LV = '/Home/difdop/DIFDOP_LV/processed/FitSpectre/{}{}_{}'
HC_beamtr_shift_gnr = '/Home/difdop/WEST_traceray_2022/{}/{}WEST_EQ005_{}_reflecto_{}_freq_all_dr-00{}.mat'
HC_beamtr_gnr = '/Home/difdop/WEST_traceray_2022/{}/{}WEST_EQ005_{}_reflecto_{}_freq_all_HC.mat'
#LV_beamtr_gnr = '/Home/difdop/WEST_traceray_2022/{}/{}WEST_EQ005_{}_reflecto_{}_freq_all_dr{}.mat'
LV_beamtr_gnr = '/Home/difdop/WEST_traceray_2022/{}/{}WEST_{}_reflecto_{}_freq_all_dr{}.mat'

#LV_beamtr_gnr = '/Home/difdop/WEST_traceray_2022/{}/{}WEST_EQ005_{}_reflecto_{}_freq_all_dr{}.mat.mat'
HC_fDop_gnr = '/Home/difdop/DIFDOP/processed/FitSpectre/{}FdF_DIFDOP_{}_HC.mat'
LV_fdop_gnr = '/Home/difdop/DIFDOP_LV/processed/FitSpectre/{}FdF_DIFDOP_X.mat'

save_fig_gnr = '/Home/FO278650/Bureau/analyzed_shot/comparison/shot{}_it{}_xmode.pdf'


# %%
def read_file(shot, xmode, isweep, user, angle = 'ver', ch = None, shift = None, flag = 0, filename_btra = None):
    global FO_fdop_gnr, FO_btra_gnr

    if user == 'FO':

        user_txt = 'Francesco Orlacchio'
        filename_fdop = FO_fdop_gnr.format(shot, ch, isweep)
        if filename_btra is None:
            filename_btra = FO_btra_gnr.format(shot, ch, xmode, isweep, angle, str(float(shift)).split('.')[1])
        #filename_btra = FO_btra_gnr.format(shot, ch, xmode, isweep, angle)

        fDop = loadmat(filename_fdop)['outp'][0]
        btra = loadmat(filename_btra)['outp'][0]
        plasma = btra['plasma'][0]
        beami = btra['beami'][0][0]
        
        header_txt = f'shot {shot} isweep {isweep} Xmode {xmode}'
        outp = dict({'user':user_txt, 'header':header_txt,\
            'fDop':fDop['fDop'][0][0], 'dfDop':fDop['dfDop'][0][0], 'val':fDop['validated'][0][0], \
            'v_perp':2*np.pi*fDop['fDop'][0][0]/(btra['k_perp'][0][0]*2e2), 'dv_perp':2*np.pi*fDop['dfDop'][0][0]/(btra['k_perp'][0][0]*2e2), \
            'k_perp': btra['k_perp'][0][0]*2, 'rho': btra['rho'][0][0], 'freqGHz': btra['freqGHz'][0][0], \
            'rho_psi': plasma['neprofilerho'][0][0][0], 'ne':plasma['neprofile'][0][0][0], \
            'angle':beami['thetadirgolareldeg'][0][0][0]\
        })
        
    elif user == 'HC':
        user_txt = 'HUGO CORVOYSIER'
        if ch == 1:
            mode = 'O'
        else:
            mode = 'X'
            
        if shift is not None:
            filename_btra = HC_beamtr_shift_gnr.format(shot, shot, xmode, isweep, shift)
            btra = loadmat(filename_btra)['out'][0][0][flag][0][0][0][0]
        else:
            filename_btra = HC_beamtr_gnr.format(shot, shot, xmode, isweep)
            btra = loadmat(filename_btra)['out'][0][0][flag][0][0][0][0]

        filename_fdop = HC_fDop_gnr.format(shot, mode)
        fDop = loadmat(filename_fdop)
        
        fdop = np.mean(np.array([fDop['fcL_fit'], fDop['fcG_fit'], fDop['fcT_fit']]), axis = 0)
        header_txt = f'shot {shot} isweep {isweep} Xmode {xmode}'
        outp = dict({'user':user_txt, 'header':header_txt,\
            'fDop':fdop[isweep-1], 'val':fDop['table'][isweep-1, 0:59], 'v_perp':2*np.pi*fdop[isweep-1]/(btra['k_perp'][0][0][0:59]*2e2),\
            'fL':fDop['fcL_fit'][isweep-1], 'fG':fDop['fcG_fit'][isweep-1], 'fT':fDop['fcT_fit'][isweep-1], \
            'dfL':fDop['dfL_fit'][isweep-1], 'dfG':fDop['dfG_fit'][isweep-1], 'dfT':fDop['dfT_fit'][isweep-1], \
            'k_perp': btra['k_perp'][0][0][0:59], 'rho':btra['rhodif'][0][0][0:59], 'freqGHz':btra['f'][0][0][0:59], \
            'rho_psi': btra['rho'][0][0], 'ne':btra['ne'][0][0], 'angle':fDop['difdop_anglepol'][0][isweep-1]\
        })
    elif user == 'LV':

        user_txt = 'LAURE VERMARE'
        #filename_btra = LV_beamtr_gnr.format(shot, shot, xmode, isweep)
        if filename_btra is None:
            if shift == 0.0:
                filename_btra = LV_beamtr_gnr.format(shot, shot, xmode, isweep, str(shift).split('.')[0] )
            else:
                filename_btra = LV_beamtr_gnr.format(shot, shot, xmode, isweep, '-'+str(shift).split('.')[1] )
        
        if isweep == 19:
            btra = loadmat(filename_btra)['out'][0][0][0][0]['d19'][0][0]
        elif isweep == 11:
            btra = loadmat(filename_btra)['out'][0][0][0][0][0][0][0]
        
        filename_fdop = LV_fdop_gnr.format(shot)
        fDop = loadmat(filename_fdop)
        
        fdop = np.mean(np.array([fDop['fcL_fit'], fDop['fcG_fit'], fDop['fcT_fit']]), axis = 0)
        header_txt = f'shot {shot} isweep {isweep} Xmode {xmode}'
        outp = dict({'user':user_txt, 'header':header_txt,\
            'fDop':fdop[isweep-1], 'val':fDop['table'][isweep-1], 'v_perp':2*np.pi*fdop[isweep-1]/(btra['k_perp'][0][0][0:59]*2e2),\
            'fL':fDop['fcL_fit'][isweep-1], 'fG':fDop['fcG_fit'][isweep-1], 'fT':fDop['fcT_fit'][isweep-1], \
            'dfL':fDop['dfL_fit'][isweep-1], 'dfG':fDop['dfG_fit'][isweep-1], 'dfT':fDop['dfT_fit'][isweep-1], \
            'k_perp': btra['k_perp'][0][0][0:59], 'rho':btra['rhodif'][0][0][0:59], 'freqGHz':btra['f'][0][0][0:59], \
            'rho_psi': btra['rho'][0][0], 'ne':btra['ne'][0][0], \
            'angle':fDop['angle'][0][isweep-1], 'drho': btra['drho'][0][0][0]\
        })
        
    return outp



#%%
summ = TIMETRACE(58333)
tt = summ.get_tt()
summ.plot_tt(isweep = [8, 24], colors = ['r', 'g', 'b'])


shots = [58333]
colors = ['r']
time = [10]
figure, ax = plt.subplots()
for shot, col, t  in zip(shots, colors, time):
    polview(shot, t, ax=ax, colors=col)
ax.legend(shots)
#%% 58108

shot, xmode, isweep, rshift, flag = 58108, 1, 9, 1, 2

out_9x = read_file(shot, 1, 9, 'FO', 'ver',ch = 2, shift = -0.0265)
out_9o = read_file(shot, 0, 9, 'FO', 'ver',ch = 1, shift = -0.0265)
out_18x = read_file(shot, 1, 18, 'FO', 'ver',ch = 2, shift = -0.024)
out_18o = read_file(shot, 0, 18, 'FO', 'ver',ch = 1, shift = -0.024)
out_12x = read_file(shot, 1, 12, 'FO', 'ver',ch = 2, shift = -0.0238)
out_12o = read_file(shot, 0, 12, 'FO', 'ver',ch = 1, shift = -0.0238)

fig, ax = plt.subplots(2, 2, figsize = (15, 11))

fig.suptitle('SHOT %d ' %(shot), size = 18)
ax[0,0].set_title('Density Profile', fontsize = 14)
ax[0,0].plot(out_9x['rho_psi'], out_9x['ne']*1e-19, 'r', label = 'it: 9')
ax[0,0].plot(out_18x['rho_psi'], out_18x['ne']*1e-19, 'b', label = 'it: 18')
ax[0,0].plot(out_12x['rho_psi'], out_12x['ne']*1e-19, 'g', label = 'it: 12')
ax[0,0].legend()
ax[0,0].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[0,0].set_ylabel(r'$n_e$ [$10^{-19}$ $m^{-3}$]', fontsize = 12)
ax[0,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[0,1].set_title(r'$k_{\perp}$ vs $\rho_\psi$', fontsize = 14)
ax[0,1].plot(out_9x['rho'][out_9x['val'] == 1], out_9x['k_perp'][out_9x['val'] == 1], 'rx', label = r'9-Xmode')
ax[0,1].plot(out_9o['rho'][out_9o['val'] == 1], out_9o['k_perp'][out_9o['val'] == 1], 'r.', label = r'9-Omode')
ax[0,1].plot(out_18x['rho'][out_18x['val'] == 1], out_18x['k_perp'][out_18x['val'] == 1], 'bx', label = r'18-Xmode')
ax[0,1].plot(out_18o['rho'][out_18o['val'] == 1], out_18o['k_perp'][out_18o['val'] == 1], 'b.', label = r'18-Omode')
ax[0,1].plot(out_12x['rho'][out_12x['val'] == 1], out_12x['k_perp'][out_12x['val'] == 1], 'gx', label = r'12-Xmode')
ax[0,1].plot(out_12o['rho'][out_12o['val'] == 1], out_12o['k_perp'][out_12o['val'] == 1], 'g.', label = r'12-Omode')
ax[0,1].legend()
ax[0,1].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[0,1].set_ylabel(r'$k_{\perp}$ $[cm^{-1}]$', fontsize = 12)
ax[0,1].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,0].set_title(r'velocity profile', fontsize = 14)
ax[1,0].plot(out_9x['rho'][out_9x['val'] == 1], out_9x['v_perp'][out_9x['val'] == 1]*1e-3, 'r.', label = 'it: 9')
ax[1,0].plot(out_9o['rho'][out_9o['val'] == 1], out_9o['v_perp'][out_9o['val'] == 1]*1e-3, 'r.')
ax[1,0].plot(out_18x['rho'][out_18x['val'] == 1], out_18x['v_perp'][out_18x['val'] == 1]*1e-3, 'b.', label = 'it: 18')
ax[1,0].plot(out_18o['rho'][out_18o['val'] == 1], out_18o['v_perp'][out_18o['val'] == 1]*1e-3, 'b.')
ax[1,0].plot(out_12x['rho'][out_12x['val'] == 1], out_12x['v_perp'][out_12x['val'] == 1]*1e-3, 'g.', label = 'it: 12')
ax[1,0].plot(out_12o['rho'][out_12o['val'] == 1], out_12o['v_perp'][out_12o['val'] == 1]*1e-3, 'g.')

ax[1,0].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[1,0].set_ylabel(r'$v_{\perp}$ $[km/s]$', fontsize = 12)
ax[1,0].legend(loc = 4)
ax[1,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,1].set_title(r'$f_{Dop}$ estimation', fontsize = 14)
ax[1,1].plot(out_12x['freqGHz'][out_12x['val'] == 1], out_12x['fDop'][out_12x['val'] == 1]*1e-3, 'rx', label = 'Xmode')
ax[1,1].plot(out_9o['freqGHz'][out_9o['val'] == 1], out_9o['fDop'][out_9o['val'] == 1]*1e-3, 'r.', label = 'Omode')
ax[1,1].set_xlabel(r'$f_{prob}$ [GHz]', fontsize = 12)
ax[1,1].set_ylabel(r'$f_{Dop}$ [kHz]', fontsize = 12)
ax[1,1].legend(loc = 1)
ax[1,1].grid(c = 'silver', ls ='--', lw = 0.5)
#
save_fig = save_fig_gnr.format(shot, isweep, xmode)
fig.savefig(save_fig)
plt.show()
#ax[1,1].set_xlim(97, 100)
# %% SHIFT

shot, xmode, isweep = 58108, 1, 12
out_FO_noshift = read_file(shot, xmode, isweep, 'FO', 'ver',ch = 2, shift = 0)

plt.plot(out_FO_noshift['rho'][out_FO_noshift['val'] == 1], out_FO_noshift['v_perp'][out_FO_noshift['val'] == 1]*1e-2, 'k-', label = 'noshift')
plt.legend()
plt.grid()
plt.xlim(1.02, 1.03)
plt.show()
#%% velocity profile 57558

shot, isweep, flag = 57558, 19, 0
Rshifts = [ -0.0308, -0.02, -0.01]
isweeps = [ 4, 11, 19]
colors = ['r', 'b', 'g']

#summ = TIMETRACE(58333)
#tt = summ.get_tt()
#summ.plot_tt(isweep = isweeps, colors = colors)

fig, ax = plt.subplots(figsize = (8, 4))

for sweep, sh, col in zip(isweeps, Rshifts, colors):
    out_FO_O = read_file(shot, 0, sweep, 'FO', 'ver', ch = 1, shift = sh)
    out_FO_X = read_file(shot, 1, sweep, 'FO', 'ver', ch = 2, shift = sh)
    drho = out_FO_X['drho']
    ax.plot(out_FO_O['rho'], out_FO_O['k_perp'],'.', c = col)
    ax.plot(out_FO_X['rho'], out_FO_X['k_perp'],'.',  c = col, label = f'it: {sweep}')
    #ax.errorbar(out_FO_O['rho'][out_FO_O['val'] == 1], out_FO_O['v_perp'][out_FO_O['val'] == 1]*1e-3, yerr = out_FO_O['dv_perp'][out_FO_O['val'] == 1]*1e-3, c = col, fmt = 'o')
    #ax.errorbar(out_FO_X['rho'][out_FO_X['val'] == 1], out_FO_X['v_perp'][out_FO_X['val'] == 1]*1e-3, yerr = out_FO_X['dv_perp'][out_FO_X['val'] == 1]*1e-3, c = col, fmt = 'o', label = f'it: {sweep}')

file = '/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch2_modex1_isweep12_inc_dr02.mat'
out = read_file(shot, 1, sweep, 'FO', 'ver', ch = 2, shift = sh, filename_btra = file)
ax.plot(out['rho'], out['k_perp'], 'r.')


ax.set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax.set_ylabel(r'$k_\perp$ [$cm^{-1}$]')
ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.legend(loc = 2)
plt.show()


#%%
file = '/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch1_modex0_isweep11_ver_dr01.mat'
out_1_O = read_file(shot, 0, sweep, 'FO', 'ver', ch = 1, shift = sh, filename_btra = file)
file = '/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch2_modex1_isweep19_ver_dr0.mat'
out_2_X = read_file(shot, 1, sweep, 'FO', 'ver', ch = 2, shift = sh, filename_btra = file)
file = '/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch1_modex0_isweep19_ver_dr0.mat'
out_2_O = read_file(shot, 0, sweep, 'FO', 'ver', ch = 1, shift = sh, filename_btra = file)

#%%
t_start = [4.6, 6, 7.6]
ne_max = [5.5, 6, 6]
dt = 0.2

fig, ax = plt.subplots(figsize = (6, 4))

for ti, sweep,sh, col, ne in zip(t_start, isweeps, Rshifts, colors, ne_max):
    
    nprof = DREFRAP(shot, flag = 0, Rshift = sh, ne_max = ne, R_max = 0,  twindow=[ti, ti + dt], verbose = False)    
    data = nprof.get_ne(verbose=False, averaged=True)
    ax.plot(data['rho_psi'], data['ne']*1e-19, c = col, label = f'it: {sweep}')
    
ax.set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax.set_ylabel(r'$n_e$ [$10^{19} m^{-3}$]')
ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.legend(loc = 0)

shots = [57558]
colors = ['r']

figure, ax = plt.subplots()
for  shot, col, t  in zip(shots, colors, t_start):
    polview(shot, t, ax=ax, colors=col)
    
ax.set_xlabel('R [m]', fontsize = 12)
ax.set_ylabel('Z [m]', fontsize = 12)
# %% DIFFERENT ITERATION
shot = 58333
isweeps = [ 8, 24]
colors = ['r', 'b']
rshifts = [-0.035, -0.035]
ang = ['inc']
#%%
fig, ax = plt.subplots(figsize = (8, 4))

ax.set_title(f'Velocity profile SHOT - {57558}', fontsize = 16)
for sweep, col in zip(isweeps, colors):
    shot, xmode, flag = 58333, 1, 0
    out_FO_X = read_file(shot, xmode, sweep, 'FO', 'inc', shift = -0.035, ch = 2)
    #out_FO_O = read_file(shot, 0, sweep, 'FO', 'inc', shift = 0 ,ch = 1)

    ax.plot(out_FO_X['rho'][out_FO_X['val'] == 1], out_FO_X['v_perp'][out_FO_X['val'] == 1]*1e-3,'.', c = col, label = f'it = {sweep}')
    #ax.plot(out_FO_O['rho'][out_FO_O['val'] == 1], out_FO_O['v_perp'][out_FO_O['val'] == 1]*1e-3,'.', c = col)


ax.legend(loc = 'center left', bbox_to_anchor = (0.8, 0.2))
ax.set_xlabel(r'$\rho_{\psi}$', fontsize = 12)
ax.set_ylabel(r'$v_{\perp}$ [km/s]', fontsize = 12)
ax.grid(c = 'silver', ls ='--', lw = 0.5)


save_fig_ang = '/Home/FO278650/Bureau/analyzed_shot/comparison/{}velocity_profile.pdf'
filename = save_fig_ang.format(shot)
#plt.savefig(filename)
plt.plot()
# %%
shot, isweep = 58333, 8
t_start = 14.2
dt = 0.6
from DBS.io.read import read_angle_pol
#%%
_time, _angle_pol_ver = read_angle_pol('west', shot, use_inclinometer= False)
_time, _angle_pol_inc = read_angle_pol('west', shot, use_inclinometer= True)


fig, ax = plt.subplots(figsize = (12, 4))
ax.set_title(f'SHOT {shot} -iteration {isweep} - dt {dt}', fontsize = 16)
ax.plot(_time, _angle_pol_inc, 'silver',  label = 'inclinometre')
ax.plot(_time, _angle_pol_ver, 'lightblue',  label = 'verin')


t_inc, ang_inc = get_angle_pol('west', shot, t_start, t_start + dt, True, return_val = False)
t_ver, ang_ver = get_angle_pol('west', shot, t_start, t_start + dt, False, return_val = False)

mean_ang_inc = get_angle_pol('west', shot,t_start, t_start + dt, True)
mean_ang_ver = get_angle_pol('west', shot,t_start, t_start + dt, False)


ax.plot(t_inc, ang_inc, 'grey', label = 'mean angle incl %.2f' %mean_ang_inc)
ax.plot(t_ver, ang_ver, 'b', label = 'mean angle verin %.2f' %mean_ang_ver)
ax.set_xlabel('time [s]', fontsize = 12)
ax.set_ylabel('angle [°]', fontsize = 12)
ax.legend(loc = 'center left', bbox_to_anchor = (0, 0.2))

ax.set_xlim(13.5, 15.5)
save_fig_ang = '/Home/FO278650/Bureau/analyzed_shot/comparison/angle_comp_shot{}_it{}_dt6_zoom.pdf'
filename = save_fig_ang.format(shot, isweep)
plt.savefig(filename)
plt.show()



# %%
shot, xmode, isweep, rshift, flag = 57558, 1, 19, 1, 2

#file0 = '/Home/difdop/WEST_traceray_2022/57558/57558WEST_1_reflecto_19_freq_all_dr0_old2.mat'
#out_LV_old = read_file(shot, xmode, isweep, 'LV', ch = 2, shift = 0.0, filename_btra = file0)

#file_git = '/Home/difdop/WEST_traceray_2022/57558/57558WEST_1_reflecto_19_freq_all_dr0_git.mat'
#out_git = loadmat(file_git)['out']['shot57558'][0][0]['d19'][0][0]

file_new = '/Home/difdop/WEST_traceray_2022/57558/57558WEST_1_reflecto_19_freq_all_dr0.mat'
out_LV_new = loadmat(file_new)['out']['shot57558'][0][0]['d19'][0][0]

out_FO = read_file(shot, xmode, isweep, 'FO', ch = 2, shift = 0.0)

file1 = "/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch2_modex1_isweep19_ver_dr0_1.mat"
out_FO_1 = read_file(shot, xmode, isweep, 'FO', ch = 2, shift = 0.0, filename_btra = file1)

#file2 = "/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch2_modex1_isweep19_ver_dr0_LVne.mat"
#out_FO_2 = read_file(shot, xmode, isweep, 'FO', ch = 2, shift = 0.0, filename_btra = file2)

#file3 = "/Home/difdop/DBSdata/processed/beamtracing/west/57558_FO/ch2_modex1_isweep19_inc_dr0_LVne.mat"
#out_FO_3 = read_file(shot, xmode, isweep, 'FO', ch = 2, shift = 0.0, filename_btra = file3)
# %%
fig, ax = plt.subplots(1, 2, figsize = (10, 4))

ax[0].plot(out_LV_new['rho'][0][0][0], out_LV_new['ne'][0][0][0]*1e-19, label = 'LV-new')
#ax[0].plot(out_LV_old['rho_psi'], out_LV_old['ne']*1e-19, 'k--', label = 'LV-old')
ax[0].plot(out_FO['rho_psi'], out_FO['ne']*1e-19, '-', label = 'FO')
#ax[0].plot(out_FO_2['rho_psi'], out_FO_2['ne']*1e-19, 'k-.', lw = 1, label = 'FO-LV')
ax[0].grid(c = 'silver', lw = 0.5)
ax[0].set_xlabel(r'$\rho_\psi$')
ax[0].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
ax[0].legend()

ax[1].plot(out_LV_new['rhodif'][0][0][0], out_LV_new['k_perp'][0][0][0]*2, '.', label = 'LV-new')
#ax[0].plot(out_LV_old['rho_psi'], out_LV_old['ne']*1e-19, 'k--', label = 'LV-old')
ax[1].plot(out_FO['rho'], out_FO['k_perp'], '.', label = 'FO_0')
ax[1].plot(out_FO_1['rho'], out_FO_1['k_perp'], '.', label = 'FO_1')
#ax[0].plot(out_FO_2['rho_psi'], out_FO_2['ne']*1e-19, 'k-.', lw = 1, label = 'FO-LV')
ax[1].grid(c = 'silver', lw = 0.5)
ax[1].set_xlabel(r'$\rho_\psi$')
ax[1].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
ax[1].legend()


fig, ax = plt.subplots()
ax.plot(out_LV_new['rhodif'][0][0][0], out_LV_new['k_perp'][0][0][0]*2, '.', label = 'LV-new')
#ax.plot(out_LV_old['rho'], out_LV_old['k_perp']*2, '.', label = 'LV-old')
ax.plot(out_FO['rho'], out_FO['k_perp'], '.', label = 'FO')
#ax[1].plot(out_FO_2['rho'], out_FO_2['k_perp'], '.', label = 'FO-LV')
ax.grid(c = 'silver', lw = 0.5)
ax.set_xlabel(r'$\rho_\psi$')
ax.set_ylabel(r'$k_\perp$ [$cm^{-1}$]')
ax.legend()

#%%
fig, ax = plt.subplots()

ax.plot(freqGHz_LV, k_perp_LV*2, '.', label = 'LV-git')
ax.plot(out_FO_1['freqGHz'], out_FO_1['k_perp'], '.', label = 'FO')
#ax.plot(out_FO_2['freqGHz'], out_FO_2['k_perp'], '.', label = 'FO-LV')
ax.grid(c = 'silver', lw = 0.5)
ax.set_xlabel(r'$f_{prob}$ [GHz]')
ax.set_ylabel(r'$k_\perp$ [$cm^{-1}$]')
ax.legend()

# %%
k_perp_0 = out_FO_1['outp'][0]['k_perp'][0][0]
dif = out_FO_1['outp'][0]['dif'][0][0]
# %%
drho = 0.05
rhoc = 0.97

#18 
combined = zip(out_18x['rho'][out_18x['val'] == 1], out_18x['v_perp'][out_18x['val'] == 1])
sort_comb = sorted(combined, key = lambda pair: pair[0])
rho18, v18 = zip(*sort_comb)
ind18 = (np.array(rho18) < rhoc + drho) & (np.array(rho18) > rhoc - drho)


ind = (out_9x['rho'] < 1 + drho) & (out_9x['rho'] > 1 - drho)
cond9 = ind & (out_9x['val'] == 1)
v9 = out_9x['v_perp'][cond18]
rho9 = np.linspace(rhoc - drho, rhoc + drho, v9.size)


plt.plot(rho9, v9*1e-19, 'r.')
plt.plot(np.array(rho18)[ind18], np.array(v18)[ind18]*1e-19, 'b.')

# %%

fig, ax = plt.subplots(figsize = (5, 3))
ax.set_title(r'$k_{\perp}$ vs $\rho_\psi$', fontsize = 14)
ax.plot(out_9x['rho'][out_9x['val'] == 1], out_9x['k_perp'][out_9x['val'] == 1], 'rx', label = r'9-Xmode')
ax.plot(out_9o['rho'][out_9o['val'] == 1], out_9o['k_perp'][out_9o['val'] == 1], 'r.', label = r'9-Omode')
ax.plot(out_18x['rho'][out_18x['val'] == 1], out_18x['k_perp'][out_18x['val'] == 1], 'bx', label = r'18-Xmode')
ax.plot(out_18o['rho'][out_18o['val'] == 1], out_18o['k_perp'][out_18o['val'] == 1], 'b.', label = r'18-Omode')
ax.legend()
ax.set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax.set_ylabel(r'$k_{\perp}$ $[cm^{-1}]$', fontsize = 12)
ax.grid(c = 'silver', ls ='--', lw = 0.5)


# %%
drho = [0.975, 1]

#sorting based on rho 
combined = zip(out_18x['rho'], out_18x['v_perp'], out_18x['val'], out_18x['k_perp'])

combine_sorted = sorted(combined,  key = lambda pair: pair[0])
rho, v_perp, val, k_perp = zip(*combine_sorted)
rho, v_perp, val, k_perp = np.array(rho), np.array(v_perp), np.array(val), np.array(k_perp)

ind = (rho[val == 1] > 0.975) & (rho[val == 1] < 1)

plt.plot(rho[ind], k_perp[ind])


#%%
shot, xmode, isweep, rshift, flag = 58333, 1, 9, 1, 2

out_8x = read_file(shot, 1, 8, 'FO', 'ver',ch = 2, shift = -0.035)
#out_8o = read_file(shot, 0, 8, 'FO', 'ver',ch = 1, shift = -0.035)
out_24x = read_file(shot, 1, 24, 'FO', 'ver',ch = 2, shift = -0.035)
#out_24o = read_file(shot, 0, 24, 'FO', 'ver',ch = 1, shift = -0.025)


fig, ax = plt.subplots(2, 2, figsize = (15, 11))

fig.suptitle('SHOT %d ' %(shot), size = 18)
ax[0,0].set_title('Density Profile', fontsize = 14)
ax[0,0].plot(out_9x['rho_psi'], out_9x['ne']*1e-19, 'r', label = 'it: 9')
ax[0,0].plot(out_18x['rho_psi'], out_18x['ne']*1e-19, 'b', label = 'it: 18')
ax[0,0].plot(out_12x['rho_psi'], out_12x['ne']*1e-19, 'g', label = 'it: 12')
ax[0,0].legend()
ax[0,0].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[0,0].set_ylabel(r'$n_e$ [$10^{-19}$ $m^{-3}$]', fontsize = 12)
ax[0,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[0,1].set_title(r'$k_{\perp}$ vs $\rho_\psi$', fontsize = 14)
ax[0,1].plot(out_9x['rho'][out_9x['val'] == 1], out_9x['k_perp'][out_9x['val'] == 1], 'rx', label = r'9-Xmode')
ax[0,1].plot(out_9o['rho'][out_9o['val'] == 1], out_9o['k_perp'][out_9o['val'] == 1], 'r.', label = r'9-Omode')
ax[0,1].plot(out_18x['rho'][out_18x['val'] == 1], out_18x['k_perp'][out_18x['val'] == 1], 'bx', label = r'18-Xmode')
ax[0,1].plot(out_18o['rho'][out_18o['val'] == 1], out_18o['k_perp'][out_18o['val'] == 1], 'b.', label = r'18-Omode')
ax[0,1].plot(out_12x['rho'][out_12x['val'] == 1], out_12x['k_perp'][out_12x['val'] == 1], 'gx', label = r'12-Xmode')
ax[0,1].plot(out_12o['rho'][out_12o['val'] == 1], out_12o['k_perp'][out_12o['val'] == 1], 'g.', label = r'12-Omode')
ax[0,1].legend()
ax[0,1].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[0,1].set_ylabel(r'$k_{\perp}$ $[cm^{-1}]$', fontsize = 12)
ax[0,1].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,0].set_title(r'velocity profile', fontsize = 14)
ax[1,0].plot(out_9x['rho'][out_9x['val'] == 1], out_9x['v_perp'][out_9x['val'] == 1]*1e-3, 'r.', label = 'it: 9')
ax[1,0].plot(out_9o['rho'][out_9o['val'] == 1], out_9o['v_perp'][out_9o['val'] == 1]*1e-3, 'r.')
ax[1,0].plot(out_18x['rho'][out_18x['val'] == 1], out_18x['v_perp'][out_18x['val'] == 1]*1e-3, 'b.', label = 'it: 18')
ax[1,0].plot(out_18o['rho'][out_18o['val'] == 1], out_18o['v_perp'][out_18o['val'] == 1]*1e-3, 'b.')
ax[1,0].plot(out_12x['rho'][out_12x['val'] == 1], out_12x['v_perp'][out_12x['val'] == 1]*1e-3, 'g.', label = 'it: 12')
ax[1,0].plot(out_12o['rho'][out_12o['val'] == 1], out_12o['v_perp'][out_12o['val'] == 1]*1e-3, 'g.')

ax[1,0].set_xlabel(r'$\rho_\psi$', fontsize = 12)
ax[1,0].set_ylabel(r'$v_{\perp}$ $[km/s]$', fontsize = 12)
ax[1,0].legend(loc = 4)
ax[1,0].grid(c = 'silver', ls ='--', lw = 0.5)

ax[1,1].set_title(r'$f_{Dop}$ estimation', fontsize = 14)
ax[1,1].plot(out_12x['freqGHz'][out_12x['val'] == 1], out_12x['fDop'][out_12x['val'] == 1]*1e-3, 'rx', label = 'Xmode')
ax[1,1].plot(out_9o['freqGHz'][out_9o['val'] == 1], out_9o['fDop'][out_9o['val'] == 1]*1e-3, 'r.', label = 'Omode')
ax[1,1].set_xlabel(r'$f_{prob}$ [GHz]', fontsize = 12)
ax[1,1].set_ylabel(r'$f_{Dop}$ [kHz]', fontsize = 12)
ax[1,1].legend(loc = 1)
ax[1,1].grid(c = 'silver', ls ='--', lw = 0.5)
#
save_fig = save_fig_gnr.format(shot, isweep, xmode)
fig.savefig(save_fig)
plt.show()