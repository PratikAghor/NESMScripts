import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
from matplotlib.ticker import FuncFormatter, MaxNLocator, ScalarFormatter
import matplotlib as mpl

##########################################
font = {'family': 'normal',
        'weight': 'bold',
        'size': 30}

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})

font = {'family': 'monospace',
        'weight': 'bold',
        'size': 20}

rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})
mpl.rcParams['hatch.linewidth'] = 3
##########################################

#-------------------------------------------------
indxRange = '952_3877'  # entire year
variables = ['kmke', 'vebf']
#-------------------------------------------------
configs = {
    'NESM_1km': {
        'prefix': 'nesm_2019_2020',
        'box': 'NE_1km',
    },
    'NESM_5km': {
        'prefix': 'nesm_2019_2020_5km',
        'box': 'NE_5km',
    },
    'SM3_1km': {
        'prefix': 'sm3_2019_2020',
        'box': 'NE_1km',
    },
}



colors = {
    'NESM_1km': 'black',
    'NESM_5km': 'tab:cyan',
    'SM3_1km': 'tab:red',
}

linewidths = {
    'NESM_1km': 8,
    'NESM_5km': 8,
    'SM3_1km':8,
}


linestyles = {
    'NESM_1km': '-',
    'NESM_5km': '-',
    'SM3_1km':'-',
}

markers = {
    'NESM_1km': '^',
    'NESM_5km': '^',
    'SM3_1km':'o',
}
markerskip=50
markersize=20

labels = {
    'NESM_1km': r'NESM 1km',
    'NESM_5km': r'NESM 5km',
    'SM3_1km': r'SM3 1km',
}

#-------------------------------------------------
# formatter for scientific notation
def scientific_formatter(x, pos):
    if x == 0:
        return '0'
    exponent = int(np.floor(np.log10(np.abs(x))))
    coeff = x / 10**exponent
    
    # Better handling of coefficients
    if abs(coeff - 1.0) < 0.01:
        return r'$10^{{{:d}}}$'.format(exponent)
    elif abs(coeff + 1.0) < 0.01:
        return r'$-10^{{{:d}}}$'.format(exponent)
    else:
        return r'${:.1f} \times 10^{{{:d}}}$'.format(coeff, exponent)
#-------------------------------------------------
#-------------------------------------------------
# #-------------------------------------------------
# only KmKe vert profile, but horizontal aspect ratio
var = 'kmke'
fig, ax = plt.subplots(1, 1, figsize=(24, 14))
    
box = configs['NESM_1km']['box'] 

for config_name, info in configs.items():
    prefix = info['prefix']
    box = info['box']
    
    filename = f"{var}/{prefix}_{var}_z_{box}_nt_{indxRange}.asc"

    data = np.loadtxt(filename)
    z = data[:, 0]
    val = data[:, 1]

    mask = ~np.isnan(val)
    z = z[mask]
    val = val[mask]

    # Interpolation
    z_fine = np.linspace(z.min(), z.max(), 1000)
    interp_func = interp1d(z, val, kind='cubic', bounds_error=False, fill_value=0)
    val_smooth = interp_func(z_fine)

    ax.plot(z_fine, val_smooth,
            color=colors[config_name],
            linewidth=linewidths[config_name],
            linestyle=linestyles[config_name],
            marker=markers[config_name],
            markersize=markersize,
            markevery=50,
            alpha=0.8,
            label=labels[config_name])

ax.axhline(y=0, color='black', linewidth=4, linestyle='--')
ax.tick_params(axis='x', labelsize=48)
ax.tick_params(axis='y', labelsize=48)
ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(48)
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))


if(var == 'kmke'):
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0), useMathText=True)
    ax.set_xticks([-5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500, -1000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda xval, p: f'${-xval/1000:.1f}$'))
    ax.text(0.93, 0.06, r'$\times 10^{3}$', transform=ax.transAxes, 
        fontsize=48, verticalalignment='top')
    # ax.yaxis.get_offset_text().set_fontsize(16)
else:
    ax.set_yticklabels([])
    ax.legend(loc=0, fontsize=48)

# ax.set_title(f"{var.upper()} vertical profile")
ax.set_ylabel(rf'$\langle KmKe \rangle$', fontsize=48)
# ax.set_xlabel('Depth [m]')
ax.set_xlim(-5000, -1500)
ax.invert_xaxis()
if(var == 'kmke'):
    ax.set_ylim(-5e-8, 1.8e-7)
elif(var == 'vebf'):
    ax.set_ylim(-1e-7, 1e-7)
ax.grid(True)
# ax.legend()

plt.tight_layout()
outname = f'profiles_{var}_box_{box}_nt_{indxRange}.png'
plt.savefig(outname, dpi=300)
plt.close()
print(f"Saved figure: {outname}")
# #-------------------------------------------------
# only vebf vert profile, but horizontal aspect ratio
var = 'vebf'
fig, ax = plt.subplots(1, 1, figsize=(24, 14))
    
box = configs['NESM_1km']['box'] 

for config_name, info in configs.items():
    prefix = info['prefix']
    box = info['box']
    
    filename = f"{var}/{prefix}_{var}_z_{box}_nt_{indxRange}.asc"

    data = np.loadtxt(filename)
    z = data[:, 0]
    val = data[:, 1]

    mask = ~np.isnan(val)
    z = z[mask]
    val = val[mask]

    # Interpolation
    z_fine = np.linspace(z.min(), z.max(), 1000)
    interp_func = interp1d(z, val, kind='cubic', bounds_error=False, fill_value=0)
    val_smooth = interp_func(z_fine)

    ax.plot(z_fine, val_smooth,
            color=colors[config_name],
            linewidth=linewidths[config_name],
            linestyle=linestyles[config_name],
            marker=markers[config_name],
            markersize=markersize,
            markevery=50,
            alpha=0.8,
            label=labels[config_name])

ax.axhline(y=0, color='black', linewidth=4, linestyle='--')
ax.tick_params(axis='x', labelsize=48)
ax.tick_params(axis='y', labelsize=48)
ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(48)
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))


ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0), useMathText=True)
ax.set_xticks([-5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500, -1000])
ax.xaxis.set_major_formatter(FuncFormatter(lambda xval, p: f'${-xval/1000:.1f}$'))
ax.text(0.93, 0.06, r'$\times 10^{3}$', transform=ax.transAxes, 
    fontsize=48, verticalalignment='top')
# ax.yaxis.get_offset_text().set_fontsize(16)

# ax.set_title(f"{var.upper()} vertical profile")
ax.set_ylabel(rf'$\langle VEBF \rangle$', fontsize=48)
# ax.set_xlabel('Depth [m]')
ax.set_xlim(-5000, -1500)
ax.invert_xaxis()
if(var == 'kmke'):
    ax.set_ylim(-5e-8, 1.8e-7)
elif(var == 'vebf'):
    ax.set_ylim(-1e-7, 1e-7)
ax.grid(True)
# ax.legend()

plt.tight_layout()
outname = f'profiles_{var}_box_{box}_nt_{indxRange}.pdf'
plt.savefig(outname, dpi=300)
plt.close()
print(f"Saved figure: {outname}")
# #-------------------------------------------------
