import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from numpy import sqrt
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
#-------------------------------------------

# Paths and settings
plots_path = './'
boxname = 'NE' # for saving figs
vlevel=-4000
w_vlevel = vlevel
indxRange = '952_3877'
w_indxRange = '0_3877'
Gamma = 0.2 # mixing efficiency

n_files_per_day = 8
if (indxRange == '0_3877'):
    start_day_for_averaging = 120
elif(indxRange == '952_3877'):
    start_day_for_averaging = 0

base_date = datetime(2019, 5, 1)
t0 = start_day_for_averaging
start_day_for_averaging_speed = 0
t0_speed = start_day_for_averaging_speed
markevery, markersize = 10, 10

# Define dataset paths and filename prefixes
sim_paths = {
    'nesm_1km': {
        'kmke': './',
        'prefix': 'nesm_2019_2020_',
        'box':'NE_1km'
    },
    'nesm_5km': {
        'kmke': './',
        'prefix': 'nesm_2019_2020_5km_',
        'box':'NE_5km'
    },
    'sm3_1km': {
        'kmke': './',
        'prefix': 'sm3_2019_2020_',
        'box':'NE_1km'
    },
}

labels = {
    'nesm_1km': r'NESM (1km)',
    'nesm_5km': r'NESM (5km)',
    'sm3_1km': r'SM3 (1km)',
}

def load_kmke(path, box, var, prefix, vlevel=vlevel, indxRange='952_3877'):
    fname = f"{prefix}{box}_box_avg_{var}_t_vlevel_{vlevel}_nt_{indxRange}.nc"
    return xr.open_dataset(path + fname)

def daily_avg(data, n):
    Nt_avg = len(data) // n
    return np.array([np.nanmean(data[i*n:(i+1)*n]) for i in range(Nt_avg)])

def compute_stats(x, label):
    mean, std = np.nanmean(x), np.nanstd(x)
    print(f'{label:25s}  mean = {mean:.4e}, std = {std:.4e}')
    return mean, std
#-----------------------------------------------------------------------------------------------
# Load and process datasets
datasets = {}
daily_avgs = {}
for label, paths in sim_paths.items():
    ds_kmke = load_kmke(paths['kmke'], paths['box'], 'KmKe', paths['prefix'], vlevel=vlevel, indxRange=indxRange)
    print(f"\nLoading data for: {label}")

    kmke = daily_avg(ds_kmke['KmKe_t'].values, n_files_per_day)
    hrs = daily_avg(ds_kmke['HRS_t'].values, n_files_per_day)
    vrs = daily_avg(ds_kmke['VRS_t'].values, n_files_per_day)

    # Time array (only once)
    if 't_arr' not in locals() and 't_arr' in ds_kmke:
        t_arr = ds_kmke['t_arr'].values

    datasets[label] = {
        'hrs': hrs,
        'vrs': vrs,
        'kmke': kmke,
    }
#-----------------------------------------------------------------------------------------------
# Date range
date_range = [base_date + timedelta(days=i) for i in range(len(next(iter(datasets.values()))['kmke']))]

# Plotting
colors = ['black', 'tab:cyan', 'tab:red']
styles = ['-', '-', '-']
linewidths = [4, 2, 2]
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
labels_list = []
std_kmke_list = []
bar_colors = ['black', 'tab:cyan', 'tab:red']

for label in sim_paths.keys():
    kmke = datasets[label]['kmke'][t0:]
    std_kmke = np.nanstd(kmke)
    std_kmke_list.append(std_kmke)
    labels_list.append(labels[label])
#------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(24, 14))

for idx, (config, data) in enumerate(datasets.items()):
    kmke = data['kmke'][t0:]
    vrs = data['vrs'][t0:]
    hrs = data['hrs'][t0:]
    # Make sure lengths align
    min_len = min(len(kmke), len(hrs))
    dates = date_range[t0:][:min_len]

    std = std_kmke_list[idx]
    s = f"{std:.2e}"
    base, exp = s.split("e")
    text = f"${base} \\times 10^{{{int(exp)}}}$"

    ax.plot(dates, kmke, color=colors[idx], linestyle=styles[idx], linewidth=linewidths[idx],
            label=rf'{labels.get(config, config)}, Std. = {text} m$^2$/s$^3$')

ax.set_ylabel(r'$\langle KmKe \rangle$', color='k', fontsize=48)
ax.tick_params(axis='y', labelcolor='k')
ax.grid(True)
ax.legend(fontsize=48, loc=2)

# Format x-axis
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right', fontsize=48)
ax.yaxis.get_offset_text().set_fontsize(48)
ax.tick_params(axis='y', labelcolor='k', labelsize=48)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(f"{plots_path}{boxname}_kmke_timeseries_daily_comparison_vlevel_{vlevel}.png", dpi=300)
#-----------------------------------------------------------------------------------------------
