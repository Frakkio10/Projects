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


# from abc import ABC, abstractmethod
from matlabtools import Struct

class WEST_wall(Struct):
    
    def __init__(self, shot):
        wall = imas_west.get(shot, 'wall')
        r_wall = wall.description_2d[0].limiter.unit[0].outline.r
        z_wall = wall.description_2d[0].limiter.unit[0].outline.z

        r_wall = np.append(r_wall, r_wall[0])
        z_wall = np.append(z_wall, z_wall[0])

        self.shot = shot
        self.r    = r_wall
        self.z    = z_wall
        
    def plot(self, ax = None, label = False):
        
        if ax is None:
            fig, ax = plt.subplots(figsize = (5, 5))
        ax.plot(self.r, self.z, c = 'k', lw = 3)
        if label:
            ax.set_title(f'West discharge #{self.shot}')
            ax.set_xlabel('R [m]', fontsize = 12)
            ax.set_ylabel('Z [m]', fontsize = 12)
        
    
    
class Equilibrium(Struct):
    
    def __init__(self, shot, t_start):
        t_ignitron = tsmat(shot, 'IGNITRON|1')[0]
        equi = imas_west.get(shot, 'equilibrium')
        time = equi.time - t_ignitron
        ind_efit = np.where(np.abs(time - t_start) == np.min(np.abs(time - t_start)))[0][0]
        time = time[ind_efit]
        rho_bord = 1.3
        
        #read parameters for the time given
        Rmag = equi.global_quantities.magnetic_axis.r[ind_efit]
        Zmag = equi.global_quantities.magnetic_axis.z[ind_efit]
        rgrid = equi.interp2D.r[:,0]
        zgrid = equi.interp2D.z[0,:]
        Rsep = equi.boundary.outline.r[ind_efit][0]
        Zsep = equi.boundary.outline.z[ind_efit][0]
        rstrike = equi.boundary.strike_point.r[ind_efit]
        zstrike = equi.boundary.strike_point.z[ind_efit]
        rxpoint = equi.boundary.x_point.r[ind_efit]
        zxpoint = equi.boundary.x_point.z[ind_efit]
        rmagaxi = equi.boundary.geometric_axis.r[ind_efit]
        zmagaxi = equi.boundary.geometric_axis.z[ind_efit]
        R0      = equi.vacuum_toroidal_field.r0
        B0      = equi.vacuum_toroidal_field.b0[ind_efit]
        
        #magnetic field 
        Bpol = np.sqrt( np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])**2 )
        Btot = np.sqrt( np.squeeze(equi.interp2D.b_field_tor[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])**2 + np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])**2 )
        Br   = np.squeeze(equi.interp2D.b_field_r[ind_efit, :, :])
        Bz   = np.squeeze(equi.interp2D.b_field_z[ind_efit, :, :])
        Btor = np.squeeze(equi.interp2D.b_field_tor[ind_efit, :, :])

        Bpol[np.isnan(Bpol)] = 0
        Btot[np.isnan(Btot)] = 0
        Br[np.isnan(Br)]     = 0
        Bz[np.isnan(Bz)]     = 0
        Btor[np.isnan(Btor)] = 0
        
        #fist psi 
        a = np.sign(equi.global_quantities.psi_axis[ind_efit] - equi.global_quantities.psi_boundary[ind_efit])
        psiaxis = equi.global_quantities.psi_axis[ind_efit]*a
        psibound = equi.global_quantities.psi_boundary[ind_efit]*a
        psi = np.squeeze(equi.interp2D.psi[ind_efit, :, :])*a
        
        self.description = f'Equilibrium for shot {shot} at time {t_start}'
        self.shot        = shot
        self.t_start     = t_start
        self.time        = time 
        self.ind_efit    = ind_efit
        self.t_ignitron  = t_ignitron
        self.Rmag        = Rmag
        self.Zmag        = Zmag
        self.rgrid       = rgrid
        self.zgrid       = zgrid
        self.Rsep        = Rsep
        self.Zsep        = Zsep
        self.Bpol        = Bpol
        self.Btor        = Btor
        self.Br          = Br
        self.Bz          = Bz
        self.Btot        = Btot
        self.psiaxis     = psiaxis
        self.psibound    = psibound
        self.psi         = psi
        self.rho_bord    = rho_bord
        self.rstrike     = rstrike
        self.zstrike     = zstrike        
        self.rxpoint     = rxpoint
        self.zxpoint     = zxpoint        
        self.rmagaxi     = rmagaxi
        self.zmagaxi     = zmagaxi
        self.R0          = R0
        self.B0          = B0

    @classmethod
    def get_equi(cls, shot, t_start):
        
        plasma = Equilibrium( shot, t_start)
        
        #remove nan from psi 
        i, j, k, l = 0, plasma.rgrid.size -1 , 0, plasma.zgrid.size -1
        while (np.sum(np.isnan(plasma.psi[i:j+1, k:l+1])) != 0):
            i, j, k, l = i+1, j-1, k+1, l-1
        ind1, ind2 = np.arange(i, j+1), np.arange(k, l+1)
        plasma.rgrid2 = plasma.rgrid[ind1]
        plasma.zgrid2 = plasma.zgrid[ind2]
        psi2 = plasma.psi[(ind1[0]):(ind1[-1] + 1), (ind2[0]):(ind2[-1] + 1)]
        plasma.psi2 = psi2.flatten('K')
        
        #creating the mesh for the interpolation
        Psi2 = np.reshape(plasma.psi2, [plasma.rgrid2.size, plasma.zgrid2.size]).T
        Rgrid2, Zgrid2 = np.meshgrid(plasma.rgrid2, plasma.zgrid2, indexing='ij')
        Rgrid3, Zgrid3 = np.meshgrid(plasma.rgrid, plasma.zgrid, indexing = 'ij')
        
        #interpolation 
        points = (plasma.rgrid2, plasma.zgrid2)
        psi = np.zeros([plasma.rgrid.size, plasma.zgrid.size])
        for i in range(0, plasma.rgrid.size):
            for j in range(0, plasma.zgrid.size):  
                psi[i,j] = interpn(points, Psi2, np.array([Rgrid3[i,0], Zgrid3[0,j]]), method = 'linear', bounds_error = False, fill_value = None)

        psi = psi.flatten('K')
        plasma.psi = psi.reshape(plasma.rgrid.size, plasma.zgrid.size)
        
        rgrid2, zgrid2 = np.meshgrid(plasma.rgrid, plasma.zgrid)
        rgrid2fin, zgrid2fin = np.meshgrid(np.linspace(plasma.rgrid.min(), plasma.rgrid.max(), 200), np.linspace(plasma.zgrid.min(), plasma.zgrid.max(), 200))

        psifin = np.zeros([rgrid2fin.shape[0], zgrid2fin.shape[0]])
        points2 = (rgrid2[0,:], zgrid2[:,0])

        for i in range(0, rgrid2fin.shape[0]):
            for j in range(0, zgrid2fin.shape[0]):
                psifin[i,j] = interpn(points2, plasma.psi, np.array([rgrid2fin[0,i], zgrid2fin[j,0]]), method = 'splinef2d')
                
                
        rmax = plasma.Rsep.max()
        rmin = plasma.Rsep.min()
        zmax = plasma.Zsep.max()
        zmin = plasma.Zsep.min()
        rr0 = (rmin+rmax)/2*100
        a03 = (rmax-rmin)/2*100
        zpos = (zmin + zmax)/2*100
        elong = (zmax - zmin)/(rmax - rmin)
        psi0 = (psifin.max()).max()
        psi1 = plasma.psibound
        dsha0 = rgrid2fin[psifin >= psi0].mean()*100 - rr0
        rhod = np.arange(0.1, 1, 0.1)

        dsha = np.zeros(rhod.size)
        for ir in range(0, rhod.size):
            dsha[ir] = np.mean(rgrid2fin[psifin >= ((1-rhod[ir]**2)*psi0 + rhod[ir]**2*psi1)])*100 - rr0
    
        dshap = np.real(np.mean(np.log(1-dsha/dsha0)/np.log(rhod)))
        
        plasma.rr0    = rr0
        plasma.a03    = a03
        plasma.zpos   = zpos
        plasma.elong  = elong
        plasma.dsha0  = dsha0
        plasma.dsha  = dsha
        plasma.dshap = dshap
        plasma.psifin = psifin
        
        return plasma
    
    def plot(self, shot, t_start, ax = None):
        
        
        wall = WEST_wall(shot)
        plasma = self.get_equi(shot, t_start)
        rho_psi = np.sqrt((plasma.psi - plasma.psiaxis) / (plasma.psibound - plasma.psiaxis))
        x,y = plasma.rgrid, plasma.zgrid
        
        if ax is None:
            fig, ax = plt.subplots(figsize = (7, 8))
        
        ax.set_title(r't = %d s -- $B_0$ = %.2f T -- $R_0$ = %.2f m' %(plasma.t_start, plasma.B0, plasma.R0), fontsize = 18)
        im = ax.contour(x, y, rho_psi.T, colors = 'k', levels = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1], linestyles = '-', origin='lower')
        im = ax.contour(x, y, rho_psi.T, colors = 'k', levels = [1.2], linestyles = '--', origin='lower')
        contour = ax.contour(x, y, rho_psi.T, colors = 'k', levels = [1.1], linestyles = '--', origin='lower')
        ax.clabel(contour, inline=True, fontsize=10, fmt="%.1f")

        contour = ax.contour(x, y, rho_psi.T, colors = 'k', levels = [0.5], linestyles = '-', origin='lower')
        ax.clabel(contour, inline=True, fontsize=10, fmt="%.1f")

        if plasma.zxpoint > 0:
            idx = np.argmax(plasma.Zsep)
        else:
            idx = np.argmin(plasma.Zsep)
            
        ax.plot(plasma.Rsep, plasma.Zsep, '--r', lw = 2)
        plt.plot([plasma.Rsep[idx], plasma.rstrike[0]], [plasma.Zsep[idx], plasma.zstrike[0]], 'r--', lw = 2)
        plt.plot([plasma.Rsep[idx], plasma.rstrike[1]], [plasma.Zsep[idx], plasma.zstrike[1]], 'r--', lw = 2)

        ax.plot(wall.r, wall.z, 'k', lw = 3)
        ax.plot(plasma.rstrike, plasma.zstrike, 'xr', markersize = 20)
        ax.plot(plasma.rxpoint, plasma.zxpoint, 'xk', markersize = 20)
       # ax.plot(plasma.rmagaxi, plasma.zmagaxi, '+k', markersize = 20)
        
        ax.set_xlabel('R [m]', fontsize = 12)
        ax.set_ylabel('Z [m]', fontsize = 12)
        fig.show()
        #return rho_psi
# %% 
if __name__ == '__main__':
    plasma = Equilibrium.get_equi(57558, 7.7)

# %%
if __name__ == '__main__':
    shots = [58333, 57558, 58108, 60269]
    t_start = [20, 6, 10, 7]
    
    for shot, ti in zip(shots, t_start):
        equi = Equilibrium(shot, ti)
        equi.plot(shot, ti)
        plt.show()
# %%
if __name__ == '__main__':
    wall = WEST_wall(57558)
# %%
