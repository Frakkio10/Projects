# %%
#%matplotlib widget
import matplotlib.pyplot as plt
import numpy as np
from pywest.shot_analysis.src.visualization import plot_shot_summary

shots = [59824] + [59785,59796]

nrow=len(shots)
ncol=1

fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True, figsize=(10.0,1.8*nrow), gridspec_kw={'hspace': 0.0})
for i,shot in enumerate(shots):
    
    if len(shots)>1: 
        ax = axs[i]
    else:
        ax = axs
    
    try:
        plot_shot_summary(shot, ax=ax, legend=False)
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"Failed to plot shot {shot}")
        
    if i==0:
        ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', borderaxespad=0., ncol=6, handlelength=0.7)

    if i==nrow-1:
        ax.set_xlabel('time [s]')

    ax.grid(True)
    
plt.tight_layout()


# %%
from pywest import polview

shots = [55562, 59824, 59785]
figure, ax = plt.subplots()
for shot, col in zip(shots, ['r', 'b','g']):
    polview(shot, time=6.5, ax=ax, colors=col)
ax.legend(shots)

