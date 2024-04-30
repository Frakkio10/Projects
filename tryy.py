#%% 
import equimap 
import numpy as np
import pandas as pd
import imas_west
import pywed as pw
from scipy.io import loadmat

#%%
shot, time, R, 
out = equimap.get(shot, time, T, Phi, Z, quantity, no_ripples = False, run = 0, occ = 0, user = 'imas_public', machine = 'west')


#%%
shot = 58333
t = 5.5
t_ignitron = pw.tsbase(shot, 'RIGNITRON')[0][0,0]
equi = imas_west.get(58333, 'equilibrium', imas_run = 0, imas_occurrence = 1, imas_user = 'imas_public', imas_machine = 'west')
time = equi.time - t_ignitron
idx_efit = np.where(np.abs(time - t) == min(np.abs(time - t)))

#%%
plasma = pd.DataFrame()

plasma.time = equi.time[idx_efit]
plasma.Rmag = equi.mag_ax_R[idx_efit]
plasma.Zmag = equi.mag_ax_Z[idx_efit]
plasma.rgrid = equi.RInterp[:,1]

# %%
equimap.get?
# %%
pw.tsmat?
# %%
mod_raf = pw.tsmat(shot, 'DREFRAP;BALAYAGE_V:mod_raf_V')
# %%
path = '/Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_60269_prof.mat'
data = loadmat(path)

tps_reflec = data['tX']
ne_reflec = data['NEX']
a, b = np.shape(ne_reflec)

# %%
