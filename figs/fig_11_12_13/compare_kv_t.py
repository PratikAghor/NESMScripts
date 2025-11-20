import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from matplotlib.ticker import ScalarFormatter
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
# Paths
plots_path = './'
box = "box" # box or domain
w_vlevel = -4000
indxRange = '952_3877'
n_files_per_day = 8
start_day_for_averaging = 0

base_date = datetime(2019, 5, 1)
t0 = start_day_for_averaging
markevery, markersize = 10, 10

sim_paths = {
    'nesm_1km': {
        'kv': './',
        'prefix': 'nesm_2019_2020_'
    },
    'nesm_5km': {
        'kv': './',
        'prefix': 'nesm_2019_2020_5km_'
    },
    'sm3_1km': {
        'kv': './',
        'prefix': 'sm3_2019_2020_'
    }
}

labels = {
    'nesm_1km': r'NESM (1km)',
    'nesm_5km': r'NESM (5km)',
    'sm3_1km': r'SM3 (1km)'
}

def load_variable(path, var, prefix):
    fname = f"{prefix}model_{box}_avg_{var}_t_orig_vlevel_{w_vlevel}_nt_{indxRange}.nc"
    return xr.open_dataset(path + fname)

def daily_avg(data, n):
    Nt_avg = len(data) // n
    return np.array([np.nanmean(data[i*n:(i+1)*n]) for i in range(Nt_avg)])

def compute_stats(x, label):
    mean, std = np.nanmean(x), np.nanstd(x)
    print(f'{label:25s}  mean = {mean:.4e}, std = {std:.4e}')
    return mean, std
#-----------------------------------------------------------------------------------------------
datasets = {}
daily_avgs = {}
for label, paths in sim_paths.items():
    ds_AKv = load_variable(paths['kv'], 'kv', paths['prefix'])
    if 't_arr' not in locals():
        t_arr = ds_AKv['t_arr'].values
    kv_raw = ds_AKv['box_avg_kv_t'].values

    kv_capped = kv_raw; # no-cap
    datasets[label] = {
        'box_avg_kv_t': daily_avg(kv_capped, n_files_per_day),
        }


# Date range
date_range = [base_date + timedelta(days=i) for i in range(len(next(iter(datasets.values()))['box_avg_kv_t']))]

#-----------------------------------------------------------------------------------------------
sim_labels = list(sim_paths.keys())
def collect_stats(key, dataset_dict):
    means = []
    stds = []
    for label in sim_labels:
        mean, std = compute_stats(dataset_dict[label][key][t0:], f'{label} {key}')
        means.append(mean)
        stds.append(std)
    return np.array(means), np.array(stds)

#-----------------------------------------------------------------------------------------------
# Collect stats
kv_means, kv_stds = collect_stats('box_avg_kv_t', datasets)

#-----------------------------------------------------------------------------------------------
# Plotting
fig, ax = plt.subplots(figsize=(18, 8), sharex=True)
colors = ['black', 'tab:cyan', 'tab:red']
styles = ['-', '-', '-']
linewidths = [4, 2, 2]


for idx, (config, data) in enumerate(datasets.items()):
    mean = kv_means[idx]
    s = f"{mean:.2e}"
    base, exp = s.split("e")
    text = f"${base} \\times 10^{{{int(exp)}}}$"

    ax.plot(date_range[t0:], data['box_avg_kv_t'][t0:], color=colors[idx], linestyle=styles[idx], linewidth=linewidths[idx],
             label=rf'{labels.get(config, config)}, Mean = {text}')

ax.set_ylabel(r'$\langle k_v \rangle [m^2/s]$', color='k', fontsize=48)
ax.tick_params(axis='y', labelcolor='k', labelsize=48) 

ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(48)
ax.tick_params(axis='y', labelcolor='k', labelsize=48)

ax.grid(True)
ax.legend(fontsize=32)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right', fontsize=48)

plt.tight_layout(rect=[0, 0.05, 1, 1])
# plt.savefig(f"{plots_path}{box}_kv_timeseries_daily_comparison_vlevel_{w_vlevel}.png", dpi=300)
plt.savefig(f"fig_12.pdf", dpi=300)

#-----------------------------------------------------------------------------------------------
