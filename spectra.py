#%% Retrieve the spectral data exported by DBSgui / fitspec.py:
import matplotlib.pyplot as plt
import numpy as np
from DBS.io.interface import DataInterface
from DBS.io.interface import get_angle_pol
from DBS.processing.sigprocessing import init_specobjs, show_spec, make_title, OutputWrapper
# a sweep that has been treated already:
machine, shot, channelval, isweep = 'west', 60269, 2, 12

data_interface = DataInterface(shot, isweep, channelval, machine=machine)
output_wrapper = OutputWrapper(data_interface) # automatically loads the existing data from file if any
ifreqs = data_interface._get_freq_choice('all')


# initialize the specobjs:
specobjs = init_specobjs(data_interface, ifreqs)
#%%
# check if some of the specobjs already have been treated (validated or rejected), in which case they are overwritten:
ifreqs_treated = output_wrapper.header.ifreqs_treated
if len(ifreqs_treated) > 0:
    print(f'Found {len(output_wrapper.header.ifreqs_treated)} specobjs already treated.')
        
    for i in ifreqs_treated:
        if i in ifreqs:
            k = (np.where(i==ifreqs))[0][0]
            specobjs[k] = output_wrapper.merged_specobj.specobjs[i-1]

output_wrapper.update_specobjs(specobjs)
#%%
for s in specobjs[:1]:
    fig, ax = plt.subplots()
    plot_dict = show_spec(s, ax=ax, verbose=True)
    make_title(ax, data_interface, s.header.ifreq)
# %% Retrieve the fit results (see also the function `show_spec` used above)
#fig, ax = plt.subplots()
from DBS.processing.fit_utils import fitfuncs


i = 0
color_fit = ['r', 'b', 'g', 'darkturquoise', 'magenta']
color_masked = ['indianred', 'royalblue']
color_raw = ['lightcoral', 'cornflowerblue', 'palegreen', 'paleturquoise', 'plum']
#%%
fig, ax = plt.subplots(figsize = (12, 6))
idx = [0]
ax.set_title(f'west, #{shot}, ch {channelval}, sweep {isweep}')
for s in  specobjs[0:8:6]:
    include_mask = s.include_mask
    xscale = s.xscale
    yscale = s.yscale
    xdata = s.f[include_mask] / xscale
    ydata = s.P[include_mask] / yscale

    curve_type = 'taylor'
    curve_func   = fitfuncs[curve_type]
    curve_params = s.fit_params[curve_type]
    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
    yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)

    x = xfit * xscale
    y = yfit * yscale + s.P_noise
    l = ax.plot(xdata * xscale / 1e6, 10 * np.log10(ydata * yscale), c = color_raw[i])
    l = ax.plot(x / 1e6, 10 * np.log10(y), c = color_fit[i], label= r'$f_{prob}$ = %.2f' %s.header.freqGHz)
    i += 1

ax.set_xlabel('f [MHz]')
ax.set_ylabel('PSD [dB]')
ax.legend(loc = 'center left', bbox_to_anchor = (0.8, 0.5))
ax.grid(c = 'silver', ls ='--', lw = 0.5)
#plt.savefig('/Home/FO278650/Bureau/Analysis/spectra/FO60269_12.pdf')
#plt.savefig('/Home/FO278650/Bureau/Analysis/spectra/FO60269_12.jpg')
# %%

def spectra(machine, shot, channelval, isweep):
    
    data_interface = DataInterface(shot, isweep, channelval, machine=machine)
    output_wrapper = OutputWrapper(data_interface) # automatically loads the existing data from file if any
    ifreqs = data_interface._get_freq_choice('all')
    specobjs = init_specobjs(data_interface, ifreqs)
    
    ifreqs_treated = output_wrapper.header.ifreqs_treated
    if len(ifreqs_treated) > 0:
        print(f'Found {len(output_wrapper.header.ifreqs_treated)} specobjs already treated.')
            
        for i in ifreqs_treated:
            if i in ifreqs:
                k = (np.where(i==ifreqs))[0][0]
                specobjs[k] = output_wrapper.merged_specobj.specobjs[i-1]

    output_wrapper.update_specobjs(specobjs)
    return specobjs

#%%
machine, shot, channelval, isweep = 'west', 60269, 2, 12
specobjs_12 = spectra(machine, shot, channelval, isweep)
machine, shot, channelval, isweep = 'west', 60269, 2, 29
specobjs_29 = spectra(machine, shot, channelval, isweep)
machine, shot, channelval, isweep = 'west', 60269, 2, 40
specobjs_40= spectra(machine, shot, channelval, isweep)
#%%
for specobjs in [specobjs_12, specobjs_29,specobjs_40]:
    ifreqs_treated = output_wrapper.header.ifreqs_treated
    if len(ifreqs_treated) > 0:
        print(f'Found {len(output_wrapper.header.ifreqs_treated)} specobjs already treated.')
            
        for i in ifreqs_treated:
            if i in ifreqs:
                k = (np.where(i==ifreqs))[0][0]
                specobjs[k] = output_wrapper.merged_specobj.specobjs[i-1]

    output_wrapper.update_specobjs(specobjs)

#%%
color_fit = ['r', 'b', 'g', 'darkturquoise', 'magenta']
color_masked = ['indianred', 'royalblue']
color_raw = ['lightcoral', 'cornflowerblue', 'palegreen', 'paleturquoise', 'plum']
dt = 0.2
ti = 5.2
t, anglepol = get_angle_pol(machine, shot, ti, ti + dt, return_val='array')
theta_12 = np.mean(anglepol)   
ti = 10.8
t, anglepol = get_angle_pol(machine, shot, ti, ti + dt, return_val='array')
theta_40 = np.mean(anglepol)
i = 0
specobjs = [specobjs_40[7], specobjs_12[7]]
theta = [theta_40, theta_12]

fig, ax = plt.subplots(figsize = (6, 6))
ax.set_title(r'SHOT #%d, $f_{prob} = %.2f$ GHz' %(shot, specobjs_29[7].header.freqGHz), fontsize = 14)
for s, th in  zip(specobjs, theta):
    include_mask = s.include_mask
    xscale = s.xscale
    yscale = s.yscale
    #xdata = s.f[include_mask] / xscale
    #ydata = s.P[include_mask] / yscale
    xdata = s.f / xscale
    ydata = s.P / yscale    

    x_masked = s.f[include_mask] / xscale
    y_masked = s.P[include_mask] / yscale
    
    curve_type = 'taylor'
    curve_func   = fitfuncs[curve_type]
    curve_params = s.fit_params[curve_type]
    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
    yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)

    x = xfit * xscale
    y = yfit * yscale + s.P_noise
    l = ax.plot(xdata * xscale / 1e6, 10 * np.log10(ydata * yscale), alpha = 0.5, c = color_raw[i])
    l = ax.plot(x_masked * xscale / 1e6, 10 * np.log10(y_masked * yscale), lw = 2, alpha = 0.8,   c = color_raw[i])
    #l = ax.plot(x / 1e6, 10 * np.log10(y), c = color_fit[i], lw =3,  label= r'$\theta$ = %.1f' %th)
    i += 1

i = 0
for s, th in  zip(specobjs, theta):
    include_mask = s.include_mask
    xscale = s.xscale
    yscale = s.yscale
    #xdata = s.f[include_mask] / xscale
    #ydata = s.P[include_mask] / yscale
    xdata = s.f / xscale
    ydata = s.P / yscale    

    x_masked = s.f[include_mask] / xscale
    y_masked = s.P[include_mask] / yscale
    
    curve_type = 'taylor'
    curve_func   = fitfuncs[curve_type]
    curve_params = s.fit_params[curve_type]
    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
    yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)

    x = xfit * xscale
    y = yfit * yscale + s.P_noise
    #l = ax.plot(xdata * xscale / 1e6, 10 * np.log10(ydata * yscale), alpha = 0.5, c = color_raw[i])
    #l = ax.plot(x_masked * xscale / 1e6, 10 * np.log10(y_masked * yscale), lw = 2, alpha = 0.8,   c = color_raw[i])
    #l = ax.plot(x / 1e6, 10 * np.log10(y), c = color_fit[i], lw =3,  label= r'$\theta$ = %.1f' %th)
    i += 1 

i = 0
for s, th in  zip(specobjs, theta):
    include_mask = s.include_mask
    xscale = s.xscale
    yscale = s.yscale
    #xdata = s.f[include_mask] / xscale
    #ydata = s.P[include_mask] / yscale
    xdata = s.f / xscale
    ydata = s.P / yscale    

    x_masked = s.f[include_mask] / xscale
    y_masked = s.P[include_mask] / yscale
    
    curve_type = 'taylor'
    curve_func   = fitfuncs[curve_type]
    curve_params = s.fit_params[curve_type]
    xfit = np.linspace(np.min(xdata), np.max(xdata), 1000)
    yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)

    x = xfit * xscale
    y = yfit * yscale + s.P_noise
    l = ax.plot(x / 1e6, 10 * np.log10(y), c = color_fit[i], lw =3,  label= r'$\theta$ = %.1f' %th)
    i += 1
    
ax.set_xlabel('f [MHz]', fontsize = 12)
ax.set_ylabel('PSD [dB]', fontsize = 12)
ax.legend(loc = 'center left', bbox_to_anchor = (0, 0.9))
ax.grid(c = 'silver', ls ='--', lw = 0.5)
ax.set_xlim(-4, 4)
plt.savefig('/Home/FO278650/Bureau/Analysis/spectra/FO60269_angle_fit.pdf')
plt.show()
# %%
specobjs_12.include_mask
# %%
