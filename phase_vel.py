#%%
import numpy as np
import matplotlib.pyplot as plt
import imas_west
import scipy.io as sc
import pywest
import pywed as pw
import pandas as pd
from DBS.io.interface import get_angle_pol
from pywest import polview
from DBS.analysis import DBS_Profile
from DBS.beamtracing import DBSbeam
import DBS
from scipy.io import loadmat
from DBS.io.interface import DataInterface

gnr_path_fdop = "/Home/difdop/DBSdata/processed/fDop_estimation/west/{}_FO/_ch{}_isweep{}.mat"
gnr_path_btr = "/Home/difdop/DBSdata/processed/beamtracing/west/{}_FO/ch{}_modex{}_isweep{}.mat"
gnr_path_savefig = '/Home/FO278650/Bureau/Analysis/fig_analised/FO{}_{}.pdf'

# %%
import shot_profile
# %%
