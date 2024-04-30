#%%
import matplotlib.pyplot as plt
import numpy as np
from pywest.shot_analysis.src.visualization import plot_shot_summary
shot = 59824

fig, ax = plt.subplots()
plot_shot_summary(shot, ax=ax)
ax.grid(True)
plt.tight_layout()

#%%

