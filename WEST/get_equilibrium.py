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
        
# %%
if __name__ == '__main__':
    
    plasma = Equilibrium.get_equi(57558, 7.7)


# %%
