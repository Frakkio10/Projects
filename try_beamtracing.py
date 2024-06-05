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
from DBS.beamtracing.src.visualize import plot_beam, plot_gola


#fro
from WEST.custom_beamtracing import custom_beam3d
from WEST.reflectometers_v1 import DREFRAP
#%%
machine, shot, twindow, modex, channelval = 'west', 58333, [23.8, 24.6], 1, 2
flag = 0
isweep = 24
output_9, interface = custom_beam3d(machine, shot, isweep , twindow, modex, channelval, flag, Rshift = -0.041,  verbose = True)

#%%
machine, shot, twindow, modex, channelval = 'west', 58333, [14.2, 14.8], 1, 2
flag = 0
isweep = 8
output_3, interface = custom_beam3d(machine, shot, isweep , twindow, modex, channelval, flag, Rshift = -0.041, verbose = True)

#%%

fig, ax = plt.subplots(figsize = (8,8))
time = 7
lab = 0
polview(shot, time=7, ax=ax, colors='coral')
for ifreq, freq in enumerate(output_3.freqGHz[0:20:2]):
    beam = output_3.beam[ifreq]
    dif = output_3.dif[ifreq]
    beami = output_3.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = output_3.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_gola(beami, ax=ax, color='k', lw=1)
    if lab == 0:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='red', label = r'$\theta$ = 3.5°')
        lab = 1
    else:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='red')
    #plot_beam(beam, dif, beami, ax=ax, central_ray=False, other_rays=True, color='silver', alpha=0.5, lw=0.5)
lab = 0

for ifreq, freq in enumerate(output_9.freqGHz[0:20:2]):
    beam = output_9.beam[ifreq]
    dif = output_9.dif[ifreq]
    beami = output_9.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = output_9.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_gola(beami, ax=ax, color='k', lw=2)
    if lab == 0:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue', label = r'$\theta$ = 9.1°')
        lab = 1
    else:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue')    
ax.legend(loc = 4)
ax.set_xlabel('R [m]', fontsize = 12)
ax.set_ylabel('z [m]', fontsize = 12)
ax.set_xlim(2.9, 3.2)
ax.set_ylim(-0.1, 0.25)
plt.savefig('/Home/FO278650/Bureau/Analysis/beamtracing/FO58333_zoom.pdf')


    
#%%
fig, ax = plt.subplots()
interface.equilibrium.plot(ax = ax)
polview(shot, time=7, ax=ax, colors='k')

# %%
 
twindow =  [7, 7.2]
cond = ((time > twindow[0]) & (time < twindow[1]))
indices = np.arange(len(time))[cond]
# %%
import statistics_west as sw
from scipy.io import loadmat
from pywest import shot_to_campaign

gnr_path_DREFRAP = '/Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_{}_prof.mat'
filemat_ne = gnr_path_DREFRAP.format(shot)
data = loadmat(filemat_ne)

t_ignitron = sw.read_signal_shot(shot, 't_ignitron',shot_to_campaign(shot)).values[0]
time = np.transpose(data['tX'][0])

cond = (time > twindow[0]) & (time < twindow[1])
indices = np.arange(len(time))
ind = indices[cond]
# %%
cond
# %%
twindow[-1]
# %%
interface.densityprofile.plot()
# %%
plt.plot(output.plasma.neprofilerho, output.plasma.neprofile, 'o-')
# %%
fig, ax = plt.subplots()
ax.set_title('shot #60269')
nprof = DREFRAP(60269, 1, [8.6, 8.8], verbose = False) #if data from DREFRAP --> add the shot
data = nprof.get_ne(verbose = False, averaged=True)
nprof.plot_ne(verbose = False, averaged=True, ax = ax)
plt.savefig('/Home/FO278650/Bureau/Analysis/Figures/Densityprof/shot60269.pdf')
# %%
fig, ax = plt.subplots()
ax.set_title('shot #60270')
nprof = DREFRAP(60270, 1, [8.6, 8.8], verbose = False) #if data from DREFRAP --> add the shot
data = nprof.get_ne(verbose = False, averaged=True)
nprof.plot_ne(verbose = False, averaged=True, ax = ax)
plt.savefig('/Home/FO278650/Bureau/Analysis/Figures/Densityprof/shot60270.pdf')

#%%
gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}.mat"
from scipy.io import loadmat

shot, isweep, channelval, modex = 58333, 8, 2, 1
filename = gnr_path_btr.format(shot, channelval, modex, isweep)

beamtr = loadmat(filename)['outp']

#%%
ifreqs = beamtr['freqGHz'][0][0][0][0:20:2]
beam = beamtr.beam[ifreq]
dif = beamtr.dif[ifreq]
beami = beamtr.beami[ifreq]
beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
integopt = beamtr.integopt[ifreq]
# %%
machine, shot, twindow, modex, channelval = 'west', 60269, [5.2, 5.4], 1, 2
flag = 1
isweep = 12
output_12, interface = custom_beam3d(machine, shot, isweep , twindow, modex, channelval, flag, Rshift = -0.041,  verbose = True)
# %%
