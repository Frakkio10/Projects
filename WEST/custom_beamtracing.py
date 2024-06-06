#%%
import numpy as np
from DBS import definitions as defs
from DBS.beamtracing import Equilibrium2d, DensityProf1d, Beam3dInterface, DIFDOP, TCVDUALV 
from WEST.reflectometers_v1 import DREFRAP
import matplotlib.pyplot as plt
from DBS.io.interface import DataInterface
from DBS.beamtracing.src.west.io import retreive_west_density_data, retrieve_west_eq
from DBS.io.interface import get_angle_pol
from DBS.beamtracing.interface import Beam3dInterface
from DBS.beamtracing.src.visualize import plot_beam, plot_gola
from pywest import polview
from pathlib import Path
from scipy.io import loadmat

gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_{}/ch{}_modex{}_isweep{}_{}_dr{}.mat"

def plot_density(nprof, data):    
    fig, ax = plt.subplots(1, 2, figsize = (10, 4))
    ax[0].plot(nprof.R + nprof.drho, nprof.ne*1e-19, 'silver', lw = 0.5)
    ax[0].plot(data['R'] + nprof.drho, data['ne_int']*1e-19, 'r')
    ax[0].set_xlabel('R [m]')
    ax[0].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
    ax[1].plot(data['rho_psi'], data['ne']*1e-19, 'r')
    ax[1].set_xlabel(r'$\rho_\psi$')
    ax[1].set_ylabel(r'$n_e$ [$10^{19}$ $m^{-3}$]')
    plt.show()

def custom_beam3d(machine, shot, isweep, drho, twindow, modex, channelval, flag, angle_choice = 'ver', user = 'FO', verbose = True):
    
    outpath = gnr_path_btr.format(shot, user, channelval, modex, isweep, angle_choice,  str(drho).split('.')[1])

    #get density profiles and prepare the construct --> density profile averaged on the time interval of the iteration 
    print('Retrieving Density Data from DREFRAP')
    nprof = DREFRAP(shot, flag, drho, twindow = twindow, verbose = True)
    data = nprof.get_ne(verbose=False, averaged=True)
    densityprof = DensityProf1d(data['rho_psi'], data['ne'], data['drho'], data['ne_max'], data['R_max'], 'west', shot, twindow, description = data['header'])
    plot_density(nprof, data)

    #retrieve or extract equi profiles and build the construct using the mean time from the time interval 
    print('Retrieving Equilibrium')
    equilibrium = retrieve_west_eq(shot, twindow, override=False, verbose=False) #using matlab script "get_equilibrium"
    time = np.mean(twindow)
    equilibrium = Equilibrium2d.from_shot(machine, shot, time) #reading from a file 
    
    #obtain info from difdop on the angle
    print('Retrieving DIFDOP data')
    dataI = DataInterface.from_time(shot, time, channelval, machine)
    t0s = dataI.get_start_time_of_freq_step()
    
    set_dt = input('Do you want to select a different dt for the angle? (y/n) Default dt = %.2f s' %(twindow[1] - twindow[0]))
    if set_dt == 'y':
        dtStep = float(input('dt to select (in s) : '))
    elif set_dt == 'n':
        dtStep = twindow[1] - twindow[0]
    else:
        print('!!Interrupting!!')
        return '', ''
        
    print('using dt = %.2f s for the mean angle' %dtStep)
    print(30*'-')
    FreqGHz = dataI.params.F
    N = len(FreqGHz)
    
    if angle_choice == 'ver':
        use_inclinometer = False
    else:
        use_inclinometer = True
        
    anglepol = [get_angle_pol(machine,shot, t0s[i], t0s[i] + dtStep, use_inclinometer) for i in range(N)]
    if machine == 'west':
        anglepol = [np.mean(anglepol)] * N
        print(f'warning: using the mean anglepol over all frequencies on WEST. dt for the average {dtStep}' ) #can be changed
    print('mean angle used = %.2f' %anglepol[0])
    
    #prepare the launcher 
    LauncherCls = DIFDOP if machine == 'west' else TCVDUALV
    launcher  = [LauncherCls(freqGHz=FreqGHz[i], anglepol=anglepol[i], modex=modex) for i in range(N)]
    
    #running and fetching the results 
    if Path(outpath).is_file():
        run_already = input("!! The beam tracing has been already performed!! Do you want to run already? y/n")
        if run_already == 'n':
            print('Fetching the results')
            interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
            outp = interface.fetch_result()
            #return outp, interface
        elif run_already == 'y':
            print('Running the beam tracing')
            interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
            interface.run_beam3d(outpath = outpath)
            outp = interface.fetch_result()
        else:
            print('!!Interrupting!!')
            return '', ''
    else:
        print('Running the beam tracing')
        interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
        interface.run_beam3d(outpath = outpath)
        outp = interface.fetch_result()
        

    return outp , interface



#%% Try 

if __name__ == '__main__':
    t_start, dt = 7.6, 0.2 #10.8
    machine, shot, twindow, modex, channelval = 'west', 57558, np.array([t_start, t_start + dt]), 1, 2
    sweep, drho = 19, 0.0
    flag = 0
    angle_choice = 'ver'
    output_3, interface_3= custom_beam3d(machine, shot, sweep, drho, twindow, modex, channelval, flag , angle_choice, verbose = True)





#%%
if __name__ == '__main__':
    machine, shot, twindow, modex, channelval = 'west', 57558, np.array([6.2, 6.4]), 1, 2
    sweep, drho = 12, -0.037
    flag = 0
    output_shift, interface_shift = custom_beam3d(machine, shot, sweep, drho, twindow, modex, channelval, flag , verbose = True)

#%% try plot 

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize = (8,8))
    time = int(np.mean(twindow))
    lab = 0
    polview(shot, time, ax=ax, colors='r')
    for ifreq, freq in enumerate(output_shift.freqGHz[0:20:2]):
        beam = output_shift.beam[ifreq]
        dif = output_shift.dif[ifreq]
        beami = output_shift.beami[ifreq]
        beami.diagnm = '' # 'tcvvx' or 'difdop', but not sure this is actually needed
        integopt = output_shift.integopt[ifreq]
        # fig, ax = plt.subplots()
        plot_gola(beami, ax=ax, color='k', lw=2)
        if lab == 0:
            plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue')
            lab = 1
        else:
            plot_beam(beam, dif, beami, ax=ax, central_ray=True, other_rays=False, color='blue')
        #plot_beam(beam, dif, beami, ax=ax, central_ray=False, other_rays=True, color='silver', alpha=0.5, lw=0.5)
    lab = 0
    ax.legend(loc = 4)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('z [m]')
    ax.set_xlim(2.8, 3.1)
    ax.set_ylim(-0.1, 0.25)

# %%
if __name__ == '__main__':
    fig, ax = plt.subplots(figsize = (8,8))
    time = int(np.mean(twindow))

    polview(shot, time, ax=ax, colors='r')

# %%
from DBS.io.read import read_angle_pol

#%%
if __name__ == '__main__':
    fig, ax = plt.subplots()
    for use_inclinometer, lab in zip([True, False], ['inclinometer', 'verin']):    
        t,y = read_angle_pol('west', 58333, use_inclinometer=use_inclinometer)
        ax.plot(t,y, label=lab); ax.legend(); ax.set_xlabel('time [s]'); ax.set_ylabel('poloidal angle [deg]')
    

# %%
