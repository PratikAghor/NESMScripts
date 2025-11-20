import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl

from datetime import datetime, timedelta
import matplotlib.dates as mdates
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
event_dates = [datetime(2020, 1, 1), datetime(2020, 1, 15)]

#--------------------------------------------
plots_path = './'

data_path_ke_1km = './'
data_path_ke_5km = './'
data_path_ke_sm3_1km = './'

#--------------------------------------------
indxRange = '0_3877'
box = "box"
n_files_per_day = 8
start_day_for_averaging = 120
#--------------------------------------------
def daily_avg(timeseries, n_files_per_day=8):
    """
    Compute daily averages.
    """
    Nt = len(timeseries)
    Nt_avg = Nt // n_files_per_day

    daily_avg_series = np.zeros(Nt_avg)
    for ll in range(0, Nt_avg):
        daily_avg_series[ll] = np.nanmean(
            timeseries[ll * n_files_per_day : (ll + 1) * n_files_per_day]
        )
    
    return daily_avg_series

#--------------------------------------------
def load_and_process_ke(data_path, sim_name, w_vlevel, indxRange, box):
    """
    Load KE data, compute daily average, and convert to KE.
    Returns the daily averaged KE timeseries.
    """
    filename = f'{sim_name}_{box}_avg_speed_t_hslice_vlevel_{w_vlevel}_nt_{indxRange}.nc'
    ds = xr.open_dataset(data_path + filename)
    
    # Get the speed data
    box_avg_speed_t = ds['box_avg_speed_t'].values
    
    # Compute daily average
    box_avg_speed_t_daily = daily_avg(box_avg_speed_t, n_files_per_day=8)
    
    # Convert to KE (0.5 * v^2)
    box_avg_ke_t_daily = 0.5 * (box_avg_speed_t_daily)**2
    
    ds.close()
    return box_avg_ke_t_daily
#--------------------------------------------
def compute_stats(x, t0=120):
    """
    Compute mean and std from day t0 onwards.
    """
    mean = np.nanmean(x[t0:])
    std = np.nanstd(x[t0:])
    return mean, std
#--------------------------------------------
results = {}
timeseries_data = {}
# Depths
depths = [-500, -4000]

for w_vlevel in depths:
    print(f"\nProcessing depth: {w_vlevel} m")
    print("-"*60)
    
    ke_nesm_1km = load_and_process_ke(
        data_path_ke_1km, 
        'nesm_2019_2020', 
        w_vlevel, 
        indxRange, 
        box
    )
    mean_1km, std_1km = compute_stats(ke_nesm_1km, start_day_for_averaging)
    print(f"NESM (1km):  mean = {mean_1km:.4e}, std = {std_1km:.4e}")
    
    ke_nesm_5km = load_and_process_ke(
        data_path_ke_5km, 
        'nesm_2019_2020_5km', 
        w_vlevel, 
        indxRange, 
        box
    )
    mean_5km, std_5km = compute_stats(ke_nesm_5km, start_day_for_averaging)
    print(f"NESM (5km):  mean = {mean_5km:.4e}, std = {std_5km:.4e}")
    
    ke_sm3_1km = load_and_process_ke(
        data_path_ke_sm3_1km, 
        'sm3_2019_2020', 
        w_vlevel, 
        indxRange, 
        box
    )
    mean_sm3, std_sm3 = compute_stats(ke_sm3_1km, start_day_for_averaging)
    print(f"SM3 (1km):   mean = {mean_sm3:.4e}, std = {std_sm3:.4e}")
    
    results[w_vlevel] = {
        'nesm_1km': (mean_1km, std_1km),
        'nesm_5km': (mean_5km, std_5km),
        'sm3_1km': (mean_sm3, std_sm3)
    }

    timeseries_data[w_vlevel] = {
        'nesm_1km': ke_nesm_1km,
        'nesm_5km': ke_nesm_5km,
        'sm3_1km': ke_sm3_1km
    }

#-----------------------------------------------------
labels = [r'NESM (1km)', r'NESM (5km)', r'SM3 (1km)']
line_colors = ['black', 'tab:cyan', 'tab:red']
bar_colors = ['grey', 'tab:cyan', 'tab:red']
#-----------------------------------------------------
# time series

base_date = datetime(2019, 1, 1)
n_days = len(timeseries_data[-500]['nesm_1km'])
date_range = [base_date + timedelta(days=i) for i in range(n_days)]
t0 = 120 # discard transients

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12), sharex=True)

# Top subplot - 500 m depth
ax1.plot(date_range[t0:], timeseries_data[-500]['nesm_1km'][t0:], 
        color=line_colors[0], linewidth=2, linestyle='-', 
        label=r'NESM (1km)')
ax1.plot(date_range[t0:], timeseries_data[-500]['nesm_5km'][t0:], 
        color=line_colors[1], linewidth=2, linestyle='-', 
        label=r'NESM (5km)')
ax1.plot(date_range[t0:], timeseries_data[-500]['sm3_1km'][t0:], 
        color=line_colors[2], linewidth=2, linestyle='-', 
        label=r'SM3 (1km)')

ax1.set_ylabel(r'$\langle KE \rangle$ (500 m)', color='k', fontsize=48)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax1.yaxis.get_offset_text().set_fontsize(48)
ax1.tick_params(axis='y', labelcolor='k', labelsize=48)
ax1.legend(fontsize=32, loc='best')
ax1.grid(True)

# highlight events on top plot
for ed in event_dates:
    ax1.axvline(ed, color='black', linestyle='--', linewidth=3)
ax1.axvspan(event_dates[0], event_dates[1], color='grey', alpha=0.3)

# Bottom subplot - 4000 m depth
ax2.plot(date_range[t0:], timeseries_data[-4000]['nesm_1km'][t0:], 
        color=line_colors[0], linewidth=2, linestyle='--',  
        marker='o', markevery=8,
        label=r'NESM (1km)')
ax2.plot(date_range[t0:], timeseries_data[-4000]['nesm_5km'][t0:], 
        color=line_colors[1], linewidth=2, linestyle='--', 
        marker='o', markevery=8,
        label=r'NESM (5km)')
ax2.plot(date_range[t0:], timeseries_data[-4000]['sm3_1km'][t0:], 
        color=line_colors[2], linewidth=2, linestyle='--', 
        marker='o', markevery=8,
        label=r'SM3 (1km)')

ax2.set_ylabel(r'$\langle KE \rangle$ (4000 m)', color='k', fontsize=48)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax2.yaxis.get_offset_text().set_fontsize(48)
ax2.tick_params(axis='y', labelcolor='k', labelsize=48)
ax2.legend(fontsize=32, loc='best')
ax2.grid(True)

# highlight events on bottom plot
for ed in event_dates:
    ax2.axvline(ed, color='black', linestyle='--', linewidth=3)
ax2.axvspan(event_dates[0], event_dates[1], color='grey', alpha=0.3)

ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=30, ha='right')
ax2.tick_params(axis='x', labelsize=40)

plt.tight_layout()
output_filename = f'{plots_path}{box}_ke_timeseries_comparison.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
#-----------------------------------------------------
# bar plot
ke_means_500 = [
    results[-500]['nesm_1km'][0],
    results[-500]['nesm_5km'][0],
    results[-500]['sm3_1km'][0]
]
ke_stds_500 = [
    results[-500]['nesm_1km'][1],
    results[-500]['nesm_5km'][1],
    results[-500]['sm3_1km'][1]
]

ke_means_4000 = [
    results[-4000]['nesm_1km'][0],
    results[-4000]['nesm_5km'][0],
    results[-4000]['sm3_1km'][0]
]
ke_stds_4000 = [
    results[-4000]['nesm_1km'][1],
    results[-4000]['nesm_5km'][1],
    results[-4000]['sm3_1km'][1]
]

# Create grouped bar plot
x = np.arange(len(labels))
width = 0.25

fig, ax = plt.subplots(figsize=(14, 12))

bars_500 = ax.bar(
    x - width/2, 
    ke_means_500, 
    width, 
    label=r'500 m',
    color=bar_colors,
    edgecolor='black',
    linewidth=1.5
)

bars_4000 = ax.bar(
    x + width/2, 
    ke_means_4000, 
    width, 
    label=r'4000 m',
    color=bar_colors,
    edgecolor='black',
    linewidth=1.5,
    hatch='//'
    )

# Annotate
max_val = max(max(ke_means_500), max(ke_means_4000))
for i in range(len(labels)):
    # 500 m annotation
    mean_500 = ke_means_500[i]
    std_500 = ke_stds_500[i]
    s_500 = f"{mean_500:.2e}"
    base_500, exp_500 = s_500.split("e")
    text_500 = f"${base_500} \\times 10^{{{int(exp_500)}}}$"
    ax.text(i - width/2, mean_500 + 1e-3, 
            text_500, ha='center', va='bottom', fontsize=32)
    
    # 4000 m annotation
    mean_4000 = ke_means_4000[i]
    std_4000 = ke_stds_4000[i]
    s_4000 = f"{mean_4000:.2e}"
    base_4000, exp_4000 = s_4000.split("e")
    text_4000 = f"${base_4000} \\times 10^{{{int(exp_4000)}}}$"
    ax.text(i + width/2, mean_4000 + 1e-3, 
            text_4000, ha='center', va='bottom', fontsize=32)

ax.set_ylabel(r'Mean $\langle KE \rangle \quad [m^2/s^2]$', fontsize=40)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(48)
ax.tick_params(axis='y', labelsize=48)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=48)
ax.legend(fontsize=30, loc = 9)
ax.grid(True, axis='y', linestyle='--', alpha=0.6)

plt.tight_layout()
output_filename = f'{plots_path}{box}_ke_grouped_barplot_depths_comparison.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Figure saved: {output_filename}")
#-----------------------------------------------------