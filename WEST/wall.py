#%%
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import numpy as np
import imas_west  # convenient package for WEST IMAS data

# from abc import ABC, abstractmethod
from matlabtools import Struct

class WEST_wall(Struct):
    
    #this class reads the files from imas_west for the walls and it saves them in a tmp 
    #folder for a given shot. Then, if the same shot is required another time it will use 
    #read the data from the file instead of imas_west to make it faster
    
    def __init__(self, shot, **kwargs):
        
        machine = kwargs.get('machine', 'west')
        p = Path('/Home/FO278650/tmp_FO/wall/wall_{}_{}.pkl'.format(machine, shot))
        
        if p.exists():
            with open(p, 'rb') as f:
                wall = pickle.load(f)
                
        else:
            
            wall = imas_west.get(shot, 'wall')
            try:
                with open(p, 'wb') as f:
                    pickle.dump(wall, f)
            except:
                print('Could not save equilibrium data to {}. Check if the directory exists and has writing permissions.'.format(p))

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
# %%
if __name__ == '__main__':
    
    wall = WEST_wall(shot = 58333)
    wall.plot()
# %%
