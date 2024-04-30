import numpy as np 
import pandas as pd 
from scipy.io import loadmat
import imas_west
import pywed as pw
import equimap

gnr_path_DREFRAP = 'Home/FC139710/WEST/drefrap/Profil/data_prof/WEST_{}_prof.mat'

def get_ne_profile(shot, t, t_ign, plasma, rhobord, dt_ne, ne_max_limit, R_max_limit):
    prof_reflecto = 3
    mod_raf = pw.tsmat(shot, 'DREFRAP;BALAYAGE_V:mod_raf_V')
    
    if prof_reflect == 1:
        density = imas_west.get(shot, 'reflectometer_profile')
        tps_reflec = density.time - t_ign
        ne_reflec = density.channel[0].n_e.data
        phi_reflec = density.channel[0],position.phi.data
        r_reflec = density.channel[0],position.r.data
        z_reflec = density.channel[0],position.z.data
        
    elif prof_reflecto == 3:
        filemat_ne = gnr_path_DREFRAP.format(shot)
        data = loadmat(filemat_ne)
        tps_reflec = data['tX']
        ne_reflec = data['NEX']
        a, b = np.shape(ne_reflec)
        phi_reflec = 2.4435*np.ones([a,b]) #??
        r_reflec = data['RX']
        z_reflec = np.zeros([a,b])
        
    if mod_raf == 1:
        ind_raf = np.where(np.diff(tps_reflec) > 1e-3)[0][0]
        traf = tps_reflec[ind_raf]
        traf = np.concatenate((traf, [tps_reflec[-1]]))
        ind_raf = np.concatenate((ind_raf, [len(tps_reflec)]))
        dd_ref = np.round(traf, t)
        ind_reflec1 = ind_raf[dd_ref] - ind_raf[0]
        ind_reflec2 = ind_rafe[dd_ref]
        ind_reflec = np.round((ind_reflec2 + ind_reflec1) / 2)
        
        R = r_reflec[:, ind_reflec1:ind_reflec2]
        ne_reflec_d = ne_reflec[:, ind_reflec1:ind_reflec2]
        sel = 0 
        
        if sel == 0 :
            
            if (ne_max_limit is None) or (R_max_limit is None):
                plt.plot(R, ne_refelc_d)
                ne_max_limit = input('Maximum density to select (in 1e19) : ')
                R_max_limit = input('Profiles have to reach at least this R value (m) : ')
            
            ind_sel = np.where((np.max(ne_reflec_d) > ne_max_limit*1e19) & (np.max(R) > R_max_limit))[0][0]
            r_sel = R[:, ind_sel]
            ne_sel = ne_reflec_d[:, ind_sel]
        
        else:
            r_sel = R
            ne_sel = ne_reflec_d
        
        a, b = np,shape(r_sel)
        NF = np.linspace(0, np.max((np.max(ne_sel))), 100)
        
        for i in range (0, b):
            RF[i,:] = np.interpN(NF, ne_sel[0:-2, i], r_sel[0:-1, i])
        
        if b == 1:
            r_reflec_m = RF
        else:
            r_reflec_m = np.mean(RF)
        
        ne_reflec_m = NF
        
        rho_pol_reflec_m = equimap.get(shot, t + t_ign, r_reflec_m, phi_reflec[:, ind_reflec], z_reflec[:, ind_reflec], 'rho_pol_norm')
        ne_reflec_m = ne_reflec_m[np.isnan(rho_pol_reflec_m)]
        R = r_reflec_m[np.isnan(rho_pol_reflec_m)]
        rho_pol_reflec_m = rho_pol_reflec_m[np.isnan(rho_pol_reflec_m)]
        
    else:
        
        dt = 0.4
        ind_reflec = np.where(abs(tps_reflec - t) <= dt)[0][0]
        tps_reflec = density.time - t_ign
        
        ind_reflec1 = np.round(tps_reflec, t)
        ind_reflec2 = np.round(tps_reflec, t + dt_ne)
        R = r_reflec[:, ind_reflec1:ind_reflec2]
        ne = ne_reflec[:, ind_reflec1:ind_reflec2]
        
        if (ne_max_limit is None) or (R_max_limit is None):
            plt.plot(R, ne, 'k')
            ne_max_limit = input('Maximum density to select (in 1e19) : ')
            R_max_limit = input('Profiles have to reach at least this R value (m) : ')
            
        ind_sel = np.where((np.max(ne_reflec_d) > ne_max_limit*1e19) & (np.max(R) > R_max_limit))[0][0]
        r_sel = R[:, ind_sel]
        ne_sel = ne_reflec_d[:, ind_sel]
    
    c, d = np.shape(r_sel)
    Rmin = np.max(np.min(r_sel), axis = 0)
    Rmax = np.max(np,max(r_sel), axis = 0)
    R_int = np.linspace(R_min, R_max, 400)
    
    for i in range (0, d):
        ne_int[:, i] = np.interp(R_int, r_sel[:, i], ne_sel[:, i], left = np.nan, right = np.nan)
    
    ne = np.mean(ne_int, 1)
    rho = equimap.get(shot, t + t_ign, R_int, np.zeros(R_int.size), np.zeros(R_int.size), 'rho_pol_norm')

    psi = np.reshape(plasma.psi, [plasma.rgrid.size, plasma.zgrid.size])
    rho_map = np.sqrt((psi - plasma-psiaxe) / (plasma.psibnd - plasma.psiaxe))
    ind_R = np.where(np.abs(plasma.zgrid) == np.min(np.abs(plasma.zgrid)))
    rho_R = rho_map[:, ind_R]
    rho = np.interp(R_int, plasma.rgrid, rho_R)

    return rho, ne
