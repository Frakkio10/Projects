#%%
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

from pywed import tsbase, tsmat
# %%
class TS_profiles:

    def __init__(self, shot):
        self.shot = shot 
        #self.rho_efit = tsbase(shot, 'SEFPSIN')
        return None 


    def DIFDOP_data(self, shot):
        DIFDOP_SRMSVIDEO = tsmat(shot, 'DIFDOP-SRMSVIDEO')
        svideo, tlent = DIFDOP_SRMSVIDEO[0], DIFDOP_SRMSVIDEO[1]
        DIFDOP_SRMSCOS = tsmat(shot, 'DIFDOP-SRMSCOS')
        srms, tlent = DIFDOP_SRMSCOS[0], DIFDOP_SRMSCOS[1]
        DIFDOP_SPOSVERIN = tsmat(shot, 'DIFDOP-SPOSVERIN')
        sposverin, tlent = DIFDOP_SPOSVERIN[0], DIFDOP_SPOSVERIN[1]
        DIFDOP_SFREQ = tsmat(shot, 'DIFDOP-SFREQ')
        sfreq, tlent = DIFDOP_SFREQ[0], DIFDOP_SFREQ[1]

        #conversuin in volt
        sposverin = np.array(sposverin) * 10 / 2048
        srms = np.array(srms) * 10 / 2048
        svideo = np.array(svideo) * .014 - 8.3

        sfreq = np.arcsin(np.array(sfreq) / 2048) * 180 / np.pi - 7.5

        if shot > 39100: #2007 campaign
            sangle = (sposverin ** 2) * 0.00785 + sposverin * 1.09 + 4.64
            if shot > 40050:
                sangle = sangle + 7
        elif shot > 34300: #2005 campaign
            sangle = (sposverin ** 2) * 0.0115 + sposverin * 1.1156 + 6.7398
        elif shot > 32040: #2003 campaign 
            sangle = (sposverin + 2.3427) / 3365 - 7.5
        else:
            sangle = []

        DIFDOP = pd.DataFrame()
        #DIFDOP['tlent'] = np.array(tlent)
        DIFDOP['svideo'] = np.array(svideo)
        DIFDOP['srms'] = np.array(srms)
        DIFDOP['sfreq'] = sfreq
        DIFDOP['sangle'] = sangle

        return DIFDOP

# %%
shot = 45511
profiles = TS_profiles(shot)
DIFDOP = profiles.DIFDOP_data(shot)
# %%
