#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import imas_west
import pywed as pw 
from scipy.io import loadmat
import equimap

gnr_path_DREFRAP = '/Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_{}_prof.mat'

# %%
shot = 58333
twindow = [14, 15]
data = imas_west.get(shot, 'reflectometer_profile')
t_ignitron = pw.tsbase(shot, 'RIGNITRON')[0][0]

time1 = data.time - t_ignitron
ne1 = data.channel[0].n_e.data
R1 = data.channel[0].position.r
Z1 = data.channel[0].position.z
Phi1 = data.channel[0].position.phi


cond = (time1 > twindow[0]) & (time1 < twindow[1])
indices = np.arange(len(time1))[cond]

#%%

def get_flattened_data(*quantities):
    
    return [q.reshape(-1,1, order='F').flatten() for q in quantities]

def interpolate(nmesh=50):
        """ interpolates the profile on a regular grid """

        R_mesh, R, ne = DREFRAP.select_on_common_R_mesh(self.R, self.ne, nmesh=nmesh, ascending_R=True)

        from scipy.interpolate import pchip_interpolate #, interp1d


def pchip_interpolation(X, Y, X_new):
    """ wrapper for pchip interpolation along first axis of the 2d arrays X and Y onto the 1d array X_new """
    Y_new = np.zeros((X_new.shape[0], Y.shape[1]))
    for i in range(X.shape[1]):
        Y_new[:,i] = pchip_interpolate(X[:,i], Y[:,i], X_new, axis=0)
    return Y_new
            
    ne_new = pchip_interpolation(R, ne, R_mesh)
        
    return R_mesh, ne_new
# %%
shot_1 = 60269

filemat_ne = gnr_path_DREFRAP.format(shot_1)
data = loadmat(filemat_ne)
t_ignitron = pw.tsbase(shot, 'RIGNITRON')[0][0]

time = data['tX'][0]
ne = np.transpose(data['NEX'])
R =  np.transpose(data['RX'])
a, b = np.shape(ne)
Phi = 2.4435*np.ones([a,b]) 
Z = np.zeros([a,b])

cond = (time > twindow[0]) & (time < twindow[1])
indices = np.arange(len(time))[cond]

#%%

time, R, Phi, Z = get_flattened_data(time, R, Phi, Z)


rho_psi = equimap.get(shot,time.mean() + t_ignitron, R, Phi, Z,'rho_pol_norm')
rho_psi = rho_psi.reshape(R.shape, order='F')
#%%
time, Phi, Z = get_flattened_data(time, Phi, Z)

R, ne_interp = interpolate(nmesh=50)

# apply horizontal shift for adjustment:
# R += Rshift
            
Phi = Phi.mean() * np.ones_like(R)
Z = Z.mean() * np.ones_like(R)


rho_psi = equimap.get(shot,time.mean() + t_ignitron, R, Phi, Z,' rho_pol_norm')
# ax.errorbar(R_mesh, ne_mean, yerr=ne_std, fmt='o', color='black', markersize=3, capsize=3, capthick=1, zorder=10)
            
ne = ne_interp.mean(axis=1)
ne_std = ne_interp.std(axis=1)
            
if twindow is None:
    twindow = [t.min(), t.max()]
                




# %%
