#%%

import matplotlib.pyplot as plt
import seaborn as sns
import statistics_west as sw

stats = sw.read_stats('C7')

#%%
stats_filt = stats[(stats['P_LH_mean_plto'] > 5e6) &
                   (stats['Ip_mean_plto'] > 400)
                   ]


# %%
stats_filt['shot']
# %%
