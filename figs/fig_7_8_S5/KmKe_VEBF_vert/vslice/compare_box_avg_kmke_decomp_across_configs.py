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
components = ['kmke', 'hrs', 'vrs']
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
    }
}



labels = {
    'NESM_1km': r'NESM (1km)',
    'NESM_5km': r'NESM (5km)',
    'SM3_1km': r'SM3 (1km)'
}

component_labels = {
    'kmke': 'KmKe',
    'hrs': 'HRS',
    'vrs': 'VRS'
}

component_colors = {
    'kmke': 'black',
    'hrs': 'tab:blue',
    'vrs': 'tab:green'
}

component_markers = {
    'kmke': None,
    'hrs': 's',
    'vrs': '^'
}
markerskip = 50
markersize= 6
#-------------------------------------------------
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
# Plot 
fig, axes = plt.subplots(1, 3, figsize=(8, 5), sharey=True)
axes = axes.flatten()
for i, (config_name, info) in enumerate(configs.items()):
    prefix = info['prefix']
    box = info['box']
    ax = axes[i]
    for var in components:
        if(var == 'kmke'):
            linewidth=4
        else:
            linewidth=3
        
        filename = f"{var}/{prefix}_{var}_z_{box}_nt_{indxRange}.asc"
        print(f"filename = {filename}")
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
        # Plot with component-specific colors and markers
        ax.plot(val_smooth, z_fine,
                color=component_colors[var],
                linewidth=linewidth,
                marker=component_markers[var],
                markevery=markerskip,
                markersize=markersize,
                alpha=0.8,
                label=component_labels[var])
    #-------------------------------------------------
    # Add vertical line at x=0
    ax.axvline(x=0, color='black', linewidth=1, linestyle='--')
    ax.set_ylim(-5000, -1000)
    ax.set_xlim(-2e-7, 2e-7)
    # Set subplot title
    # ax.set_title(labels[config_name])
    # ax.set_xlabel('Value')
    ax.grid(True)
    if(i == 0):
        ax.legend(loc=1)
    
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    # ax.set_yticks([-5000, -4000, -3000, -2000, -1000])
    # ax.set_yticklabels([5000, 4000, 3000, 2000, 1000])
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0), useMathText=True)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
    ax.xaxis.get_offset_text().set_fontsize(16)

    if(i==0):
        ax.set_yticks([-5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500, -1000])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda yval, p: f'${-yval/1000:.1f}$'))
        ax.text(-0.15, 1.1, r'$\times 10^{3}$', transform=ax.transAxes, 
            fontsize=16, verticalalignment='top')

plt.tight_layout()
plt.savefig(f'kmke_decomp_box_{box}_nt_{indxRange}.png', dpi=300, bbox_inches='tight')
plt.show()

#-------------------------------------------------
