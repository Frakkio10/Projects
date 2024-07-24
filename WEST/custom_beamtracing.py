#%%
import numpy as np
from DBS import definitions as defs
from DBS.beamtracing import Equilibrium2d, DensityProf1d, Beam3dInterface, DIFDOP, TCVDUALV 
from DBS.beamtracing.custom.reflectometers_v1 import DREFRAP
from DBS.beamtracing.interface import Beam3dInterface
import matplotlib.pyplot as plt
from DBS.io.interface import DataInterface
from DBS.beamtracing.src.west.io import retreive_west_density_data, retrieve_west_eq
from DBS.io.interface import get_angle_pol
from DBS.beamtracing.interface import Beam3dInterface
from DBS.beamtracing.src.visualize import plot_beam, plot_gola
from pywest import polview
from pathlib import Path
from scipy.io import loadmat
from DBS.beamtracing.src.visualize import plot_beam, plot_gola

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

def get_twindow(shot, isweep, channelval):
    dataI = DataInterface(shot, isweep, channelval, machine='west')
    t0s = dataI.get_start_time_of_freq_step()
    dtStep = np.round(t0s[-1] - t0s[0], 1)
    twindow = twindow = [t0s[0].round(2), t0s[0].round(2) + dtStep]
    
    return np.mean(twindow), twindow

def get_isweep(shot, twindow, channelval):
    time = np.mean(twindow)
    dataI = DataInterface.from_time(shot=shot, time=time, channelval=channelval, machine='west', verbose=False)
    isweep = dataI.isweep
    
    return isweep
    
def custom_beam3d(machine, shot, drho, modex, channelval, isweep = None, twindow = None, flag = 0, angle_choice = 'ver', user = 'fe', verbose = True):

    if twindow is None and isweep is None:
        raise NotImplementedError('!! Either isweep or twindow must be specified, but none was given !!')
    elif twindow is None:
        dataI = DataInterface(shot, isweep, channelval, machine='west')
        twindow = dataI.get_sweep_twindow().round(2)
    elif isweep is None:
        isweep = get_isweep(shot, twindow, channelval)
    
    outpath = gnr_path_btr.format(shot, user, channelval, modex, isweep, angle_choice,  str(float(drho)).split('.')[1])
    if not Path(outpath).parent.exists():
        Path(outpath).parent.mkdir(parents=True)
        
    #get density profiles and prepare the construct --> density profile averaged on the time interval of the iteration 
    print('Retrieving Density Data from DREFRAP', flush = True)
    nprof = DREFRAP(shot, flag = flag, drho = drho, twindow = twindow, verbose = True)
    data = nprof.get_ne(verbose=False, averaged=True)
    densityprof = DensityProf1d(data['rho_psi'], data['ne'], 'west', shot, twindow, description = data['header'])
    plot_density(nprof, data)

    #retrieve or extract equi profiles and build the construct using the mean time from the time interval 
    print('Retrieving Equilibrium', flush = True)
    equilibrium = retrieve_west_eq(shot, twindow, override=False, verbose=False) #using matlab script "get_equilibrium"
    time = np.mean(twindow)
    equilibrium = Equilibrium2d.from_shot(machine, shot, time) #reading from a file 
    
    #obtain info from difdop on the angle and asking the user which dt he wants to use for the angle 
    print('Retrieving DIFDOP data', flush = True)
    dataI = DataInterface.from_time(shot, time, channelval, machine)
    t0s = dataI.get_start_time_of_freq_step()
    dtStep = 0.2
    
    dtStep = input('Select a different dt for the angle? Default dt = %.2f s' %(twindow[1] - twindow[0]))
    if dtStep == '':
        dtStep = 0.2
    elif float(dtStep):
        dtStep = float(dtStep)
    else:
        print('!!Interrupting!!')
        return '', ''
        
    print('using dt = %.2f s for the mean angle' %dtStep, flush = True)
    print(30*'-')
    FreqGHz = dataI.params.F
    N = len(FreqGHz)
    
    if angle_choice == 'ver':
        use_inclinometer = False
    else:
        use_inclinometer = True
        
    anglepol = get_angle_pol(machine, shot, twindow[0], twindow[0] + dtStep, return_val='mean', use_inclinometer = use_inclinometer,)
    anglepol = [anglepol] * N

    '''dtStep = np.diff(t0s).mean()
    anglepol = [get_angle_pol(machine,shot, t0s[i], t0s[i] + dtStep) for i in range(N)]
    print(anglepol)
    anglepol = [np.mean(anglepol)] * N'''
    print('warning: using the mean anglepol over all frequencies on WEST. dt for the average %.3f' %dtStep) #can be changed
    print('mean angle used = %.4f' %anglepol[0], flush = True)
    
    #prepare the launcher 
    LauncherCls = DIFDOP if machine == 'west' else TCVDUALV
    launcher  = [LauncherCls(freqGHz=FreqGHz[i], anglepol=anglepol[i], modex=modex) for i in range(N)]
    
    print(outpath)
    #running and fetching the results 
    if Path(outpath).is_file():
        run_already = input("!! The beam tracing has been already performed!! Do you want to run already? y/n")
        if run_already == 'n':
            print('Fetching the results', flush = True)
            interface = Beam3dInterface(shot, launcher, densityprof, equilibrium ,outpath=outpath)
            outp = interface.fetch_result()
            #return outp, interface
        elif run_already == 'y':
            print('Running the beam tracing', flush = True)
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


# %%
if __name__ == '__main__':
    machine, shot, sweep, drho,  modex, channelval = 'west', 57558, 19, 0, 1, 2
    angle_choice, user = 'ver', 'FO'  #these variables are optional
    output_3, interface_3 = custom_beam3d(machine, shot, drho, modex, channelval, isweep = None, twindow = [7.6, 7.8], flag = 0, user = user, verbose = True)



# %%
