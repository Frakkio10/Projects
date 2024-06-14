#%%
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.io import loadmat 
import imas_west

from DBS.io.interface import get_angle_pol

from DBS.io.interface import extract_DBS_files

# %%
shot_i = 58000
shot_f = 58100

shots = np.arange(shot_i, shot_f, 1)
# %%
with open('/Home/FO278650/Bureau/Analysis/database.txt', 'w') as file:
    for shot in shots:

        _, succ = extract_DBS_files(int(shot), machine='west',  extract_all=True)
        
        if succ == 1:
            t, ang = get_angle_pol('west', int(shot), 0, 60, t_averaged = False, use_inclinometer = True)
            m = (ang[-1] - ang[0])/(t[-1] - t[0])
        
            if m > 0.01 or m < -0.01:
                print(shot)
                file.write(f'{shot}\n')
            


# %%
