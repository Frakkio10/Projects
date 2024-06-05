#%%
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.io import loadmat 
import os 
import matlab.engine 

#WEST 
import imas_west
import pywest 
import pywed as pw
from pywed import tsbase
from pywed import tsmat
from matlabtools import Struct

#TORESUPRA
from TS_script.summary_TS import TIMETRACE as TT_ts, DIFDOP

#%%
class PROFILES():
    
    def __init__(self, shot):
        self.shot = shot
        R0, time_R0 = tsbase(shot, 'SRMAJ', nargout = 2)
        a, time_a = tsbase(shot, 'SAEQA', nargout = 2)
        am, time_am = tsbase(shot, 'SAMIN', nargout = 2)
        Z0, time_Z0 = tsbase(shot, 'SZPOS', nargout = 2)
        d0, time_d0 = tsbase(shot, 'SD0MAG', nargout = 2)
        piqd, time_piqd = tsbase(shot, 'SSMAG', nargout = 2)
        q_fit, time_qfit = tsbase(shot, 'GEFQ', nargout = 2)
        zeff, time_zeff = tsbase(shot, 'szfbrm', nargout = 2)
        rho_fit = tsbase(shot, 'SEFPSIN', nargout = 1)
        header_txt = 'Plasma information'
        plasma = dict({'header':header_txt, 'R0':R0[:,0], 'time_R0':time_R0[:,0], 'a':a[:,0], 'time_a':time_a[:,0], 'am':am[:,0], 'time_am':time_am[:,0], 'Z0':Z0[:,0], 'time_Z0':time_Z0[:,0], 'd0':d0[:,0], 'time_d0':time_d0[:,0], 'piqd':piqd[:,0], 'time_piqd':time_piqd[:,0], 'q_fit':q_fit[:,0], 'time_qfit':time_qfit[:,0], 'zeff':zeff[:,0], 'time_zeff':time_zeff[:,0], 'rho_fit':rho_fit[:,0]})
        self.plasma = plasma 
        
        #initializing matlab
        os.system('module load tools_dc')
        os.system('matlab -nodesktop')
        eng = matlab.engine.start_matlab()
        eng.addpath('/Home/FO278650/Bureau/Analysis/TS_script')
        self.eng = eng
    
    def rz2rho(self, r, z, a, r0, z0, d0, piqd = 0, el = 0, epl = 0):
        eng = self.eng
        r = matlab.double(r.tolist())
        z = matlab.double(z.tolist())
        if (piqd == 0) and (el == 0) and (epl == 0):
            rho, error = eng.rz2rho(r, z, matlab.double([a]), matlab.double([r0]), matlab.double([z0]), matlab.double([d0]), nargout = 2)
        else: 
            rho, error = eng.rz2rho(r, z, matlab.double([a]), matlab.double([r0]), matlab.double([z0]), matlab.double([d0]), matlab.double([piqd]), matlab.double([el]), matlab.double([epl]), nargout = 2)
        
        return rho, error
    #mean plasma paramters
    def mean_plasma(self, time_1, time_2):
        #mean between t1\   
        i1 = np.argmin(np.abs(self.plasma['time_a'] - time_1))
        i2 = np.argmin(np.abs(self.plasma['time_a'] - time_2))
        ie1 = np.argmin(np.abs(self.plasma['time_qfit'] - time_1))
        ie2 = np.argmin(np.abs(self.plasma['time_qfit'] - time_2))
        
        a_mean    = self.plasma['a'][i1:i2].mean()
        Z0_mean   = self.plasma['Z0'][i1:i2].mean()
        R0_mean   = self.plasma['R0'][i1:i2].mean()
        d0_mean   = self.plasma['d0'][i1:i2].mean()
        piqd_mean = self.plasma['piqd'][i1:i2].mean()
        q_mean    = self.plasma['q_fit'][ie1:ie2].mean()
        
        header_txt = f'Mean plasma parameters in [{time_1},{time_2}]'
        plasma_mean = dict({'header':header_txt, 'a':a_mean, 'Z0':Z0_mean,'R0':R0_mean,'d0':d0_mean,'piqd':piqd_mean,'q':q_mean})
        self.plasma_mean = plasma_mean
        return plasma_mean
    #Ece data 
    def ECE(self):
        Te_ece, time_ece, poubelle, poubelle, R_ece = tsbase(self.shot, 'GSHTE', 'GSHR', nargout = 5)
        if Te_ece.shape[0] == 0:
            Te_ece, time_ece, poubelle, poubelle, R_ece = tsbase(self.shot, 'GSHTENV','GSHR', nargout = 5)
        header_txt = 'ECE data'
        ece = dict({'header':header_txt, 'shot':self.shot, 'Te':Te_ece, 'time':time_ece[:,0], 'R':R_ece})
        self.ece = ece
        return ece
    #Thomson data
    def THOMSON(self):
        ne_thom, time_thom, z_thom, c_thom = tsbase(self.shot,'GNETHOM', nargout = 4)
        if time_thom.shape[0] == 0:
            print('!!NO Thomson data!!')
            return
        header_txt = 'Thomson data'
        thomson = dict({'header':header_txt, 'shot':self.shot, 'ne':ne_thom, 'time':time_thom[:,0], 'z':z_thom, 'c':c_thom})
        self.thomson = thomson
        return thomson
    #Interferometer data
    def INTERFEROMETER(self):
        interf, time_interf, rho_interf = tsbase(self.shot, 'GNE', nargout = 3)
        if time_interf.shape[0] == 0:
            print('!!NO Interferometer data!!')
            return
        header_txt = 'Interferometer data'
        interf = dict({'header':header_txt, 'shot':self.shot, 'interf':interf, 'time':time_interf[:,0], 'rho':rho_interf})
        self.interf = interf
        return interf
    #DREFLEC (reflectometer) data
    def DREFLEC(self):
        ne_dreflec, time_dreflec= tsbase(shot, 'GREFNEX', nargout = 2)
        r_dreflec, time_dreflec = tsbase(shot, 'GREFRX', nargout = 2)
        if time_dreflec.shape[0] == 0:
            print('!!NO DREFLEC data!!')
            check = input("Do you want to check in the local directory /reflec/TREFLEX_profils/ ? (1 or 0)")
            if check == 0:
                return
            elif check == 1:
                filename = '/home/sccp/gttm/reflec/TREFLEX_profils/choc{}_prof.mat'
                filename = filename.format(self.shot)
                data = loadmat(filename)
                ne_dreflec = data['NEX']
                time_dreflec = data['tps']
                r_dreflec = data['RX'] 
        header_txt = 'DREFLEC data'
        dreflec = dict({'header':header_txt, 'shot':self.shot, 'ne':ne_dreflec, 'time':time_dreflec[:,0], 'r':r_dreflec})
        self.dreflec = dreflec
        return dreflec
    #DREFLUC (reflectometer) data
    def DREFLUC(self):
        ne_drefluc_ref, time_drefluc_ref, ii = tsbase(self.shot, 'GNEFLUCREF', nargout = 3)
        r_drefluc_ref, time_drefluc_ref, ii = tsbase(self.shot, 'GRFLUCREF', nargout = 3)
        if time_drefluc_ref.shape[0] == 0:
            print('!!NO DREFLUC data!!')
        header_txt = 'DREFLUC data'
        drefluc = dict({'header':header_txt, 'shot':self.shot, 'ne':ne_drefluc_ref,  'time':time_drefluc_ref[:,0],  'r':r_drefluc_ref})
        self.drefluc = drefluc
        return drefluc

    def correction_ECE(self, time_1, time_2, cor):
        ece = self.ece
        shot = self.shot
        plasma_mean = self.mean_plasma(time_1, time_2)
        #check shot and calculate the mean 
        if ece['Te'].shape[0] != 0:
            ind_ece1 = np.argmin(np.abs(ece['time'] - time_1))
            ind_ece2 = np.argmin(np.abs(ece['time'] - time_2))
        if cor == 1:
            if (shot == 43076) or (shot == 43085) or (shot == 43082) or \
            (shot == 43083) or (shot == 43080) or (shot == 43087):
                Te_ece = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28])
                r_ece = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28])
            if (shot > 45475) and (shot < 45513):
                Te_ece = ece['Te'][ind_ece1:ind_ece2, 0:28]
                r_ece = ece['R'][ind_ece1:ind_ece1, 0:28] 
            elif (shot == 43089) or (shot == 43090) or (shot == 43100) or \
            (shot == 43092) or (shot == 43096) or (shot == 43097):
                Te_ece = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28], ece['Te'][ind_ece1:ind_ece2, 30:31])
                r_ece = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28], ece['R'][ind_ece1:ind_ece2, 30:31])
            elif (shot == 43081) or (shot == 43072)  or \
            (shot == 43069) or (shot == 43094):
                Te_ece = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28])
                r_ece = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28])
            elif shot == 43379:
                Te_ece = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:5], ece['Te'][ind_ece1:ind_ece2, 6:31])
                r_ece = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:5], ece['R'][ind_ece1:ind_ece2, 6:31])
            Te_ece_mean = np.mean(Te_ece, axis = 0)
            r_ece_mean = np.mean(r_ece, axis = 0)
        else:
            Te_ece_mean = np.mean(ece['Te'][ind_ece1:ind_ece2, :], axis = 0)
            r_ece_mean = np.mean(ece['R'][ind_ece1:ind_ece2, :], axis = 0)
        
        print(r_ece_mean)
        z_ece = np.zeros([1, r_ece_mean.size])
        rho_ece, error = self.rz2rho(r_ece_mean, z_ece, plasma_mean['a'], plasma_mean['R0'], plasma_mean['Z0'], plasma_mean['d0'], plasma_mean['piqd'])
        return rho_ece, error

# %%
if __name__ == '__main__':
    shot = 45511
    fig, ax = plt.subplots(figsize = (15,4))
    tt = TT_ts(shot)
    tt.plot(ax = ax)
    difdop = DIFDOP(shot)
    difdop.plot(ax = ax)
    #time_trace.plot_timetrace(shot)
    #plasma = tt.plasma(shot)

    
#%%
if __name__ == '__main__':
    shot = 45511
    prof = PROFILES(shot)
    ece_data = prof.ECE()
    rho_ece, error = prof.correction_ECE(6, 8, 0)

# %%
