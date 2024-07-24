#%%

from DBS.io.DBSsync import push_processed_data
# %%
shots = [58333, 57558, 60269, 58108, 60270]

for shot in shots:
    push_processed_data('west', shot, data = 'beamtracing', dest_hostname = 'altair1')
    push_processed_data('west', shot, data = 'fDop_estimation', dest_hostname = 'altair1')
    
# %%
