#%%
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd

from pywed import tsmat, tsbase
import statistics_west as sw
import equimap
import imas_west  # convenient package for WEST IMAS data
from pywest import shot_to_campaign
from DBS.io.utils import suppress_output
from scipy.io import loadmat 
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interpn
from scipy.interpolate import RegularGridInterpolator



# from abc import ABC, abstractmethod
from matlabtools import Struct
# %%

shot = 57558 
rho_bord = 1.3
t_ignitron = tsmat(shot, 'IGNITRON|1')[0]
t = 7.6

equi = imas_west.get(shot, 'equilibrium')
time = equi.time - t_ignitron
ind_efit = np.where(np.abs(time - t) == np.min(np.abs(time - t)))[0][0]
print(ind_efit)
#%%
time = time[ind_efit]
Rmag = equi.global_quantities.magnetic_axis.r[ind_efit]
Zmag = equi.global_quantities.magnetic_axis.z[ind_efit]
rgrid = equi.interp2D.r[:,0]
zgrid = equi.interp2D.z[0,:]
#%%
#magnetic field
Bpol = np.sqrt( np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])**2 )
Bpol[np.isnan(Bpol)] = 0
Btot = np.sqrt( np.squeeze(equi.interp2D.b_field_tor[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])**2 )
Btot[np.isnan(Btot)] = 0
Br = np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])
Br[np.isnan(Br)] = 0
Bz = np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])
Bz[np.isnan(Bz)] = 0
Btor = np.squeeze(equi.interp2D.b_field_tor[ind_efit, :, :])
Btor[np.isnan(Btor)] = 0

#%%
a = np.sign(equi.global_quantities.psi_axis[ind_efit] - equi.global_quantities.psi_boundary[ind_efit])
psiaxis = equi.global_quantities.psi_axis[ind_efit]*a
psibound = equi.global_quantities.psi_boundary[ind_efit]*a
psi = np.squeeze(equi.interp2D.psi[ind_efit, :, :])*a
 
# %%
i, j, k, l = 0, rgrid.size -1 , 0, zgrid.size -1
    
while (np.sum(np.isnan(psi[i:j+1, k:l+1])) != 0):
    i, j, k, l = i+1, j-1, k+1, l-1
        
ind1, ind2 = np.arange(i, j+1), np.arange(k, l+1)
# %%
rgrid2 = rgrid[ind1]
zgrid2 = zgrid[ind2]
psi2 = psi[(ind1[0]):(ind1[-1] + 1), (ind2[0]):(ind2[-1] + 1)]

#%%
psi2 = psi2.flatten('K')

# %%
Psi2 = np.reshape(psi2, [rgrid2.size, zgrid2.size]).T
Rgrid2, Zgrid2 = np.meshgrid(rgrid2, zgrid2, indexing='ij')
Rgrid3, Zgrid3 = np.meshgrid(rgrid, zgrid, indexing = 'ij')

# %%
points = (rgrid2, zgrid2)
psi = np.zeros([rgrid.size, zgrid.size])
interp = RegularGridInterpolator(points, Psi2, method = 'linear', bounds_error = False, fill_value = None)

for i in range(0, rgrid.size):
    for j in range(0, zgrid.size):  
        psi[i,j] = interp(np.array([Rgrid3[i,0], Zgrid3[j,0]]))

#%%
psi = psi.flatten('K')
# %%
psi = psi.reshape(rgrid.size, zgrid.size)
# %%
# 
rgrid2, zgrid2 = np.meshgrid(rgrid, zgrid)
Rgrid2fin, Zgrid2fin = np.meshgrid(np.linspace(rgrid.min(), rgrid.max(), 2000), np.linspace(zgrid.min(), zgrid.max(), 2000))
# %%
points = (rgrid2[0,:], zgrid2[:,0])
psifin = np.zeros([rgrid.size, zgrid.size])
for i in range(0, rgrid.size):
    for j in range(0, zgrid.size):  
        psifin[i,j] = interpn(points, psi, np.array([Rgrid2fin[0,i], Zgrid2fin[j,0]]), 'splinef2d')
# %%
psi0 = n