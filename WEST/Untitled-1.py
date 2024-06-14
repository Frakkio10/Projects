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

from matlabtools import Struct
from WEST.get_equilibrium import Equilibrium
#%%

wall = imas_west.get(57558, 'wall')
r_wall = wall.description_2d[0].limiter.unit[0].outline.r
z_wall = wall.description_2d[0].limiter.unit[0].outline.z

r_wall = np.append(r_wall, r_wall[0])
z_wall = np.append(z_wall, z_wall[0])

# %%
equilb = Equilibrium(57558, 7.7)
equi = equilb.get_equi()
# %%

levels_req = np.linspace(np.nanmin(equi['psi']), np.nanmax(equi['psi']), 30)
#levels_req = [-0.8, -0.6, -0.4, -0.2]
levels_req = np.array([0.2,0.4,0.6,0.8,1.0, 1.1, 1.2])*-1
fig, ax =plt.subplots(figsize = (5, 5))
ax.plot(r_wall, z_wall, 'k', lw = 2)
ax.contour(equi['rgrid'], equi['zgrid'], equi['psi'].T, levels = np.flip(levels_req), colors = 'silver', ls = '-', lw = 0.5)
ax.plot(equi['Rsep'], equi['Zsep'], c = 'r')
ax.plot(equi['Rmag'], equi['Zmag'], ls='', marker='+', color='k', lw=2)

# %%
Phi = 2.4435*np.ones([equi['rgrid'].size,equi['zgrid'].size])
rho_psi = equimap.get(57558,equi['time'] + equi['t_ignitron'], equi['rgrid'], Phi, equi['zgrid'], 'rho_pol_norm') 

# %%
