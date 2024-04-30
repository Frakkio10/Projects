"""
Example of how to customize a beamtracing call (the stzeps that are hidden during the DBSbeam call). This is useful for debugging or if you want to run the beamtracing in a different way than the default way, or for fictituous runs (e.g. for scenario development).
"""
#%%
import numpy as np
from DBS import definitions as defs
from DBS.beamtracing import Equilibrium2d, DensityProf1d, Beam3dInterface, DIFDOP, TCVDUALV 
from fro.reflectometers_v1 import DREFRAP
import matplotlib.pyplot as plt
from DBS.io.interface import DataInterface
from DBS.beamtracing.src.west.io import retreive_west_density_data, retrieve_west_eq
from DBS.io.interface import get_angle_pol
from DBS.beamtracing.interface import Beam3dInterface
from DBS.beamtracing.src.visualize import plot_beam, plot_gola
from pywest import polview


#%%
def custom_beam3d(machine, shot, twindow, modex, channelval, flag,  verbose = True):
    
    #get density profiles and prepare the construct
    nprof = DREFRAP(shot, flag, twindow, verbose = False) #if data from DREFRAP --> add the shot
    data = nprof.get_ne(verbose, averaged=True)
    densityprof = DensityProf1d(data['rho_psi'], data['ne'], 'west', shot, twindow, description = data['header'])

    #retrieve or extract equi profiles  and build the construct 
    equil = retrieve_west_eq(shot, twindow, override=False, verbose=False)
    time = np.mean(twindow)
    equilibrium = Equilibrium2d.from_shot(machine, shot, time)
    
    #obtain info from difdop on the angle
    dataI = DataInterface.from_time(shot, time, channelval, machine)
    t0s = dataI.get_start_time_of_freq_step()
    dtStep = np.diff(t0s).mean()
    FreqGHz = dataI.params.F
    N = len(FreqGHz)

    anglepol = [get_angle_pol(machine,shot, t0s[i], t0s[i] + dtStep) for i in range(N)]
    if machine == 'west':
        anglepol = [np.mean(anglepol)] * N
        print('warning: using the mean anglepol over all frequencies on WEST') #can be changed
    
    #prepare the launcher 
    LauncherCls = DIFDOP if machine == 'west' else TCVDUALV
    launcher  = [LauncherCls(freqGHz=FreqGHz[i], anglepol=anglepol[i], modex=modex) for i in range(len(FreqGHz))]

    #create the path to save the date
    outpath = defs.get_path_for('beam3d_output', machine=machine,
                            shot=shot, xmode=modex,channelval=channelval, isweep=dataI.isweep)

    if not outpath.parent.exists():
        outpath.parent.mkdir(parents=True)

    #running and fetching the results 
    interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
    interface.run_beam3d()
    outp = interface.fetch_result()
    
    return outp , interface



#%%
machine, shot, twindow, modex, channelval = 'west', 58333, [14.2, 14.8], 1, 2
flag = 0
output_3, interface_3 = custom_beam3d(machine, shot, twindow, modex, channelval, flag , verbose = True)

#%%
machine, shot, twindow, modex, channelval = 'west', 58333, [24.4, 25], 1, 2
flag = 0
output_8, interface_8 = custom_beam3d(machine, shot, twindow, modex, channelval, flag , verbose = True)

output

#%%
fig, ax = plt.subplots(figsize = (8,8))
time = 7
lab = 0
polview(shot, time=7, ax=ax, colors='r')
for ifreq, freq in enumerate(output.freqGHz[0:20:2]):
    beam = output_3.beam[ifreq]
    dif = output_3.dif[ifreq]
    beami = output_3.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = output_3.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_gola(beami, ax=ax, color='k', lw=2)
    if lab == 0:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue', label = r'$\theta$ = 3°')
        lab = 1
    else:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue')
    #plot_beam(beam, dif, beami, ax=ax, central_ray=False, other_rays=True, color='silver', alpha=0.5, lw=0.5)
lab = 0
for ifreq, freq in enumerate(output.freqGHz[0:20:2]):
    beam = output_8.beam[ifreq]
    dif = output_8.dif[ifreq]
    beami = output_8.beami[ifreq]
    beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
    integopt = output_8.integopt[ifreq]
    # fig, ax = plt.subplots()
    plot_gola(beami, ax=ax, color='k', lw=2)
    if lab == 0:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='green', label = r'$\theta$ = 8°')
        lab = 1
    else:
        plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='green')    
ax.legend(loc = 4)
ax.set_xlabel('R [m]')
ax.set_ylabel('z [m]')
ax.set_xlim(2.9, 3.2)
ax.set_ylim(-0.1, 0.25)

    
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
# %%
