#%%
from WEST.custom_beamtracing import custom_beam3d

#%%
t_start, dt = 10.8, 0.2
machine, shot, twindow, modex, channelval = 'west', 60269, [t_start, t_start + dt], 1,2
flag = 1 #flag = 0 for data from imas_west
sweep = 40
angle_choice = 'ver'
Rshift = 0.0
outp1, interface1 = custom_beam3d(machine, shot, sweep, Rshift, twindow, modex, channelval, flag, angle_choice, verbose = True)
#done
# %%
t_start, dt = 4.6, 0.2
machine, shot, twindow, modex, channelval = 'west', 57558, [t_start, t_start + dt], 0,1
flag = 0 #flag = 0 for data from imas_west
sweep = 4
angle_choice = 'ver'
Rshift = -0.0308
outp1, interface1 = custom_beam3d(machine, shot, sweep, Rshift, twindow, modex, channelval, flag, angle_choice, verbose = True)
#done

# %%
