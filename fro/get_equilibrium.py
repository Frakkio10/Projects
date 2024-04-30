#%%
import numpy as np 
import imas_west
import pywed as pw
import pandas as pd 
from matlabtools import Struct

class data_source():
    
    def __init__(self, shot, time = None, run = 0, occ = 1, user = 'imas_public', machine = 'west', **kwargs):
        
        self.shot, self.time = shot, time
        self.run, self.occ = run, occ
        self.user, self.machine = user, machine 
        t_ignitron = pw.tsmat(shot, 'IGNITRON|1')[0]
        self.t_ignitron = t_ignitron
        
        equi = imas_west.get(shot, 'equilibrium', imas_run = 0, imas_occurrence = 1, imas_user = 'imas_public', imas_machine = 'west')
        
        mask_eq = np.asarray(equi.code.output_flag) >= 0 
        
        if time is None:
            #choose the mit time 
            time = 0.5 * (equi.time[mask_eq].max() + equi.time[mask_eq].min())
        
        idx_time = np.argmin(np.abs(equi.time[mask_eq] - t_ignitron - time)) 
        
        r = equi.interp2D.r
        z = equi.interp2D.z
                
        psi = np.squeeze(equi.interp2D.psi[mask_eq][idx_time, ...])
        psi_0 = equi.global_quantities.psi_axis[idx_time]
        psi_LCFS = equi.global_quantities.psi_boundary[idx_time]
        psi_norm = (psi - psi_0) / (psi_LCFS - psi_0)
        rho_psi = np.sqrt( psi_norm )
        
        r_axis = equi.global_quantities.magnetic_axis.r[idx_time]
        z_axis = equi.global_quantities.magnetic_axis.z[idx_time]
        
        self.time = equi.time[mask_eq][idx_time]
        self.equi = equi 
        self.r, self.z = r, z
        self.psi, self.psi_0, self.psi_LCFS, self.psi_norm = psi, psi_0, psi_LCFS, psi_norm
        self.rho_psi, self.r_axis, self.z_axis = rho_psi, r_axis, z_axis
        

    
    
        
        
#%% 
#shot = 58333
#equi = data_source(shot, time = 5.0)


# %%
def West_chamber(self, get_wall_from_imas = False):
        
    wall = imas_west.get(shot, 'wall')
        
    r_wall = wall_descrition_2d[0].limiter.unit[0].outline.r
    z_wall = wall_descrition_2d[0].limiter.unit[0].outline.z
        
    self.r_wall, self.z_wall = r_wall, z_wall
        
        
#%%
shot = 60270
equi = data_source(shot, time = 5.0)

# %%
