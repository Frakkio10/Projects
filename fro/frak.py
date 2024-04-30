import imas_west
import numpy as np 
import pandas as pd 
import pywed as pw

class get_data():
    
    def __init__(self, shot):
        
        self.shot = shot 
        self.tignitron = pw.tsbase(shot, 'RIGNITRON')[0][0,0]

    def summary(self, shot ):
        summ = imas_west.get(shot, 'summary')
        
        summary_df = pd.DataFrame()
        summary_df = pd.DataFrame()
        summary_df['time'] = summ.time - t_ignitron 
        summary_df['ip'] = summ.global_quantities.ip.value 
        summary_df['p_ohm'] = summ.global_quantities.power_ohm.value 
        summary_df['n_e'] = summ.line_average.n_e.value 
        summary_df['p_ic'] = summ.heating_current_drive.power_ic.value 
        summary_df['p_lh'] = summ.heating_current_drive.power_lh.value 
        
        return summary_df
    


#%%
get_data.summary(58333)
# %%
