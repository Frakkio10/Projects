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
# %%

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

        #zeff, nl, Pohm, ICRH
        zeff, time_zeff = tsbase(shot, 'szfbrm', nargout = 2)
        self.zeff = self.check_data(tref, time_zeff, zeff)
        nl, time_nl = tsbase(shot,'gnl%4', nargout = 2)
        self.nl = self.check_data(tref, time_nl, nl)
        Pohm, time_Pohm = tsbase(shot, 'spohm', nargout = 2)
        self.Pohm = self.check_data(tref, time_Pohm, Pohm)
        #Picrh, time_Picrh = tsbase(shot, 'GBILAN%3', nargout = 2)
        Picrh, time_Picrh = tsbase(shot,'GPUIFCI%4', nargout = 2)
        self.Picrh = self.check_data(tref, time_Picrh, Picrh)
        self.tref = tref
        
    def check_data(self, X_ref, X, Y):
        if Y.shape[0] != 0:
            Y_new = np.interp(X_ref, X[:,0], Y[:, 0])
        else:
            Y_new = np.zeros(X_ref.shape[0]) 
        return Y_new

    def plasma(self):
        shot = self.shot
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
    
    def plot(self, ax = None):
        shot = self.shot
        if ax is None:
            fig, ax = plt.subplots(figsize = (15, 4))
        ax.set_title('TORESUPRE SHOT #%d, B = %.1f T' %(shot, self.B0))
        ax.plot(self.tref, self.Pohm, 'dodgerblue', label = r'$P_{ohm}$ [MW]')
        ax.plot(self.tref, self.Ip, 'magenta', label = r'$I_p$ [100 kA]')

        ax.plot(self.tref, self.Picrh, 'r', label = r'$P_{ICRH}$ [MW]')
        ax.plot(self.tref, self.nl * 1e-19, 'b', label = r'$n_e$ [$10^{19}$ $m^{-3}$]')
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
        ax.plot(self.difdop['time'], self.difdop['sangle'], label = 'angle')
        ax.set_xlim(0, self.difdop['plasma_dur'] + 2)
        ax.grid(c = 'silver', ls ='--', lw = 0.5)
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        ax.set_xlabel('time [s]')
        
# %%
if __name__ == '__main__':
    shot = 45511
    fig, ax = plt.subplots(figsize = (15,4))
    tt = TIMETRACE(shot)
    tt.plot(ax = ax)
    difdop = DIFDOP(shot)
    difdop.plot(ax = ax)
    #time_trace.plot_timetrace(shot)
    #plasma = tt.plasma(shot)
# %%
