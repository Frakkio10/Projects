#%%
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.io import loadmat 

#WEST 
import imas_west
import pywest 
import pywed as pw
from pywed import tsbase
from pywed import tsmat
from matlabtools import Struct


class TIMETRACE():
    
    def __init__(self, shot):
        self.shot = shot
        #reading plasma parameters
        self.plasma_dur = tsbase(shot, 'RDUREE', nargout = 1)[0][0]
        self.B0 = tsbase(shot, 'RBTOR', nargout = 1)[0][0] 
        self.Ip, time_Ip = tsbase(shot, 'SIMAG', nargout = 2) 
        #ref time
        if time_Ip.shape[0] != 0:
            tref = time_Ip[:,0]
            l = tref.shape[0]
        else:
            print('!!No Ip data!!')
    
        def check_data(X_ref, X, Y):
            if Y.shape[0] != 0:
                Y_new = np.interp(X_ref, X[:,0], Y[:, 0])
            else:
                Y_new = np.zeros(X_ref.shape[0]) 
            return Y_new
        #zeff, nl, Pohm, ICRH
        zeff, time_zeff = tsbase(shot, 'szfbrm', nargout = 2)
        self.zeff = check_data(tref, time_zeff, zeff)
        nl, time_nl = tsbase(shot,'gnl%4', nargout = 2)
        self.nl = check_data(tref, time_nl, nl)
        Pohm, time_Pohm = tsbase(shot, 'spohm', nargout = 2)
        self.Pohm = check_data(tref, time_Pohm, Pohm)
        #Picrh, time_Picrh = tsbase(shot, 'GBILAN%3', nargout = 2)
        Picrh, time_Picrh = tsbase(shot,'GPUIFCI%4', nargout = 2)
        self.Picrh = check_data(tref, time_Picrh, Picrh)
        self.tref = tref

    def plasma(self, shot):
        R0, time_R0 = tsbase(shot, 'SRMAJ', nargout = 2)
        a, time_a = tsbase(shot, 'SAEQA', nargout = 2)
        am, time_am = tsbase(shot, 'SAMIN', nargout = 2)
        Z0, time_Z0 = tsbase(shot, 'SZPOS', nargout = 2)
        d0, time_d0 = tsbase(shot, 'SD0MAG', nargout = 2)
        piqd, time_piqd = tsbase(shot, 'SSMAG', nargout = 2)
        q_fit, time_qfit = tsbase(shot, 'GEFQ', nargout = 2)
        rho_fit = tsbase(shot, 'SEFPSIN', nargout = 1)
        header_txt = 'Plasma information'
        plasma = dict({'header':header_txt, 'R0':R0[:,0], 'time_R0':time_R0[:,0], 'a':a[:,0], 'time_a':time_a[:,0], 'am':am[:,0], 'time_am':time_am[:,0], 'Z0':Z0[:,0], 'time_Z0':time_Z0[:,0], 'd0':d0[:,0], 'time_d0':time_d0[:,0], 'piqd':piqd[:,0], 'time_piqd':time_piqd[:,0], 'q_fit':q_fit[:,0], 'time_qfit':time_qfit[:,0], 'rho_fit':rho_fit[:,0]})
        self.plasma = plasma 
        return plasma 
    
    def plot_timetrace(self, shot, ax = None):
        
        if ax is None:
            fig, ax = plt.subplots(figsize = (15, 4))
        ax.set_title('TORESUPRE SHOT #%d, B = %.1f T' %(shot, self.B0))
        ax.plot(self.tref, self.Pohm, 'dodgerblue', label = r'$P_{ohm}$')
        ax.plot(self.tref, self.Ip, 'magenta', label = r'$I_p$')

        ax.plot(self.tref, self.Picrh, 'r', label = r'$P_{ICRH}$')
        ax.plot(self.tref, self.nl * 1e-19, 'b', label = r'$n_e$')
        #ax.plot(self.tref, self.zeff, 'dodgerblue', label = r'$z_{eff}$')
        ax.set_xlim(0, self.plasma_dur + 2)
        ax.set_ylim(0, 10)
        ax.grid(c = 'silver', ls ='--', lw = 0.5)
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_xlabel('time [s]')
        
class DIFDOP():
    
    def __init__(self, shot):
        plasma_dur = tsbase(shot, 'RDUREE', nargout = 1)[0][0]
        svideo, tlent, srms, tlent, sposverin, tlent, sfreq, tlent = tsmat(shot, 'DIFDOP-SRMSVIDEO', 'DIFDOP-SRMSCOS', 'DIFDOP-SPOSVERIN', 'DIFDOP-SFREQ', nargout = 8)
        #conversion in volts 
        sposverin = sposverin[:,0] * 10 / 2048
        srms = srms[:,0] * 10 / 2048
        svideo = svideo[:,0] * 0.014 - 8.3  # x 10 dBm calibration AT
        sfreq = np.arcsin( sfreq[:,0] / 2048) * 180 / np.pi - 7.5 #inclinometre since 2006
        #calibration of the verrain signal 
        if shot > 39100: #campagne 2007 
            sangle = (sposverin**2) * 0.00785 + sposverin * 1.09 + 4.64 # 2007
            if shot > 40050:
                sangle = sangle + 7
        elif shot > 34300: #campagne 2005
            sangle = (sposverin**2)  *0.0115 + sposverin * 1.1156 + 6.7398
            print(3)
        elif shot > 32040: #campagne 2003
            sangle = (sposverin + 2.3427) / 0.3365 - 7.5
            print(4)
        else:
            sangle=[] #any calibration on the motor
        header_txt = 'DIFDOP data'
        difdop = dict({'header':header_txt, 'shot':shot, 'plasma_dur':plasma_dur, 'time':tlent[:,0], 'svideo':svideo, 'srms':srms, 'sposverin':sposverin, 'sfreq':sfreq, 'sangle':sangle})
        self.difdop = difdop
        
    def get_difdop(self, shot):
        return self.difdop
    
    def plot(self, ax = None):
        if ax is None:
            fig, ax = plt.subplots(figsize = (15, 4))
            
        #ax.plot(self.difdop['time'], self.difdop['svideo'])
        #ax.plot(self.difdop['time'], self.difdop['srms'])
        ax.plot(self.difdop['time'], self.difdop['sfreq'], 'silver', label = 'inclinometer')
        #ax.plot(self.difdop['time'], self.difdop['sangle'], label = 'angle')
        ax.set_xlim(0, self.difdop['plasma_dur'] + 2)
        ax.grid(c = 'silver', ls ='--', lw = 0.5)
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_xlabel('time [s]')


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
    
    #ECE data 
    def ECE(self):
        Te_ece, time_ece, poubelle, poubelle, R_ece = tsbase(self.shot, 'GSHTE', 'GSHR', nargout = 5)
        if Te_ece.shape[0] == 0:
            Te_ece, time_ece, poubelle, poubelle, R_ece = tsbase(self.shot, 'GSHTENV','GSHR', nargout = 5)
        header_txt = 'ECE data'
        ece = dict({'header':header_txt, 'shot':self.shot, 'Te':Te_ece, 'time':time_ece[:,0], 'R':R_ece})
        self.ece = ece
        return ece
    
    def THOMSON(self):
        ne_thom, time_thom, z_thom, c_thom = tsbase(self.shot,'GNETHOM', nargout = 4)
        if time_thom.shape[0] == 0:
            print('!!NO Thomson data!!')
            return
        header_txt = 'Thomson data'
        thomson = dict({'header':header_txt, 'shot':self.shot, 'ne':ne_thom, 'time':time_thom[:,0], 'z':z_thom, 'c':c_thom})
        self.thomson = thomson
        return thomson

    def INTERFEROMETER(self):
        interf, time_interf, rho_interf = tsbase(self.shot, 'GNE', nargout = 3)
        if time_interf.shape[0] == 0:
            print('!!NO Interferometer data!!')
            return
        header_txt = 'Interferometer data'
        interf = dict({'header':header_txt, 'shot':self.shot, 'interf':interf, 'time':time_interf[:,0], 'rho':rho_interf})
        self.interf = interf
        return interf
    
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

    def DREFLUC(self):
        ne_drefluc_ref, time_drefluc_ref, ii = tsbase(self.shot, 'GNEFLUCREF', nargout = 3)
        r_drefluc_ref, time_drefluc_ref, ii = tsbase(self.shot, 'GRFLUCREF', nargout = 3)
        if time_drefluc_ref.shape[0] == 0:
            print('!!NO DREFLUC data!!')
        header_txt = 'DREFLUC data'
        drefluc = dict({'header':header_txt, 'shot':self.shot, 'ne':ne_drefluc_ref,  'time':time_drefluc_ref[:,0],  'r':r_drefluc_ref})
        self.drefluc = drefluc
        return drefluc
    
    def rz2rho(self, r, z, a, r0, z0, d0, piqd = 0, el = 0, epl = 0):
        #chack the input arguments
        mode2 = 1
        elip = 0
        if np.array(piqd).size == 0:
            elip = 0
            mode2 = 1
        elif (np.array(el).size == 0) or (np.array(epl).size == 0):
            elip = 0
            mode2 = 1
        else:
            elip = 1
            mode2 = 0
        #check the input size
        if (not np.array_equal(a.size, r0.size)) or (not np.array_equal(a.size, z0.size)) or (not np.array_equal(a.size, d0.size)) or \
            ((mode2 == 0) and (not np.array_equal(a.size, np.array(piqd).size))) or ((elip == 1) and (not np.array_equal(a.size, np.array(el).size))) or \
            ((elip == 1) and (not np.array_equal(a.size, np.array(epl).size))):
            
            return print('!!ERROR the plasma geometry data must have all the same dimension.')
        
        if (not np.array_equal(r.size, z.size)):
            return print('!!ERROR r and z must have all the same dimension.')
        
        if (r.size == 1):
            r = r * np.ones(a.size)
            z = z * np.ones(a.size)
        elif (r.size == 1):
            r = np.ones(a.size) * r
            z = np.ones(a.size) * z
        elif (r.size != a.size ):
            return print('!!ERROR the dimensione 1 of r (or z) must be the same of a, r0 or z0')
        
        if (a.size == 1):
            a = a * np.ones([1, r.size])
            r0 = r0 * np.ones([1, r.size])
            z0 = z0 * np.ones([1, r.size])
            d0 = d0 * np.ones([1, r.size])
            if np.array(piqd).size != 0:
                piqd = piqd * np.ones([1, r.size])
            if np.array(el).size != 0:
                el = el * np.ones([1, r.size])
                epl = epl * np.ones([1, r.size])
        elif (a.size != r.size):
            return print('!!ERROR the dimensione 2 of either r0, z0, d0 must be either 1 or the same as r')
        
        if (np.array(el).size == 0) or (elip == 0) or (mode2 == 1):
            epl = np.zeros(a.shape)
            el = np.ones(a.shape)
        
        if (np.array(piqd).size == 0) or (mode2 == 1):
            piqd = 2 * np.ones(a.shape)
        #end of checks 
        rho = np.nan * np.ones(r.shape)
        error = rho
        test = not (np.isnan(r) or np.isnan(z))
        print(test)
        if np.any(test):
            rho[test], error[test] = self.rz2rho(r[test, :], z[test, :], a[test, :], r0[test], z0[test], d0[test], piqd[test], el[test], epl[test])
        
        matnan = rho < 0 
        cpt = np.nan * np.ones([np.sum(np.sum(matnan)), 1])
        rho[matnan] = cpt[:, 0]
        
        return rho, error
        
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
                Te_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28])
                r_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28])
            if (shot > 45475) and (shot < 45513):
                Te_ece[ind_ece1:ind_ece2, :] = ece['Te'][ind_ece1:ind_ece2, 0:28]
                r_ece[ind_ece1:ind_ece2, :] = ece['R'][ind_ece1:ind_ece1, 0:28] 
            elif (shot == 43089) or (shot == 43090) or (shot == 43100) or \
            (shot == 43092) or (shot == 43096) or (shot == 43097):
                Te_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28], ece['Te'][ind_ece1:ind_ece2, 30:31])
                r_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28], ece['R'][ind_ece1:ind_ece2, 30:31])
            elif (shot == 43081) or (shot == 43072)  or \
            (shot == 43069) or (shot == 43094):
                Te_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:26], ece['Te'][ind_ece1:ind_ece2, 28])
                r_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:26], ece['R'][ind_ece1:ind_ece2, 28])
            elif shot == 43379:
                Te_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['Te'][ind_ece1:ind_ece2, 0:5], ece['Te'][ind_ece1:ind_ece2, 6:31])
                r_ece[ind_ece1:ind_ece2, :] = np.concatenate(ece['R'][ind_ece1:ind_ece1, 0:5], ece['R'][ind_ece1:ind_ece2, 6:31])
            te_ece_mean = np.mean(Te_ece[ind_ece1:ind_ece2, :])
            r_ece_mean = np.mean(r_ece[ind_ece1:ind_ece2, :])
        else:
            te_ece_mean = np.mean(ece['Te'][ind_ece1:ind_ece2, :])
            r_ece_mean = np.mean(ece['R'][ind_ece1:ind_ece2, :])
        
        z_ece = np.zeros([1, r_ece_mean.size])
        rho_ece, error = self.rz2rho(r_ece_mean, z_ece, plasma_mean['a'], plasma_mean['R0'], plasma_mean['Z0'], plasma_mean['d0'], plasma_mean['piqd'])
        return rho_ece, error 
                  
# %%
if __name__ == '__main__':
    shot = 45511
    #tt = TIMETRACE(shot)
    prof = PROFILES(shot)
    #time_trace.plot_timetrace(shot)
    #plasma = tt.plasma(shot)
    plasma_mean = prof.plasma_mean(4, 6)

#%%
    difdop = DIFDOP(shot)
    data = difdop.get_difdop(shot)
    difdop.plot()
    
    fig, ax = plt.subplots(figsize = (15, 4))
    time_trace.plot_timetrace(shot,ax = ax)
    difdop.plot(ax = ax)
    
#%%
if __name__ == '__main__':
    shot = 45511
    prof = PROFILES(shot)
    ece_data = prof.ECE()
    rho_ece, error = prof.correction_ECE(5, 6, 0)
# %%
