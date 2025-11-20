import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
from datetime import datetime, timedelta
from scipy.stats import gaussian_kde
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

# -------------------------------------------------------------------------------------
# configs and paths
sim_paths = {
    'nesm_1km': {
        'pv': './',
        'prefix': 'nesm_2019_2020_'
    },
    'nesm_5km': {
        'pv': './',
        'prefix': 'nesm_2019_2020_5km_'
    },
    'sm3_1km': {
        'pv': './',
        'prefix': 'sm3_2019_2020_'
    }
}

labels = {
    'nesm_1km': r'NESM (1km)',
    'nesm_5km': r'NESM (5km)',
    'sm3_1km': r'SM3 (1km)',
}

colors = ['black', 'tab:cyan', 'tab:red']
styles = ['-', '-', '-']
linewidths = [4, 2, 2]

# -------------------------------------------------------------------------------------
indxRange = '952_3877'
z1 = -4000
z2 = -1000
pv_thresh = 1e-9
n_files_per_day = 8
start_day = 0  # to skip spin-up
base_date = datetime(2019, 5, 1) 
# ----------
sci_notation = f"{pv_thresh:.0e}"
base, exponent = sci_notation.split('e')
exponent = int(exponent)
if base == '1':
    formatted_thresh = f"10^{{{exponent}}}"
else:
    formatted_thresh = f"{base} × 10^{{{exponent}}}"
# ----------
# -------------------------------------------------------------------------------------
def daily_avg(timeseries, n_files_per_day=8):
    Nt = len(timeseries)
    Nt_avg = Nt // n_files_per_day
    if Nt % n_files_per_day != 0:
        print("saving daily timeseries.")
    daily_avg_series = np.zeros(Nt_avg)
    for ll in range(Nt_avg):
        daily_avg_series[ll] = np.nanmean(
            timeseries[ll * n_files_per_day : (ll + 1) * n_files_per_day]
        )
    return daily_avg_series

# -------------------------------------------------------------------------------------
results = {}

for sim, info in sim_paths.items():
    prefix = info['prefix']

    path = info['pv']
    label = labels[sim]
    print("sim = ", sim)
    filename = f"{path}{prefix}IS_V_t_z1_{z1}_z2_{z2}_nt_{indxRange}.nc"
    
    
    ds = xr.open_dataset(filename)
    IS_V = ds['IS_V'].transpose('Nt').values
    fq   = ds['fq'].transpose('Nt').values
    IS_V_daily = daily_avg(IS_V, n_files_per_day=8)
    fq_daily = daily_avg(fq, n_files_per_day=8)

    date_range = [base_date + timedelta(days=i) for i in range(len(fq_daily))]

    results[sim] = {
        'date': date_range,
        'q_abs_lt_thresh_vol': fq_daily,
        'IS_V_daily': IS_V_daily,
    }

# -------------------------------------------------------------------------------------
# mean IS_V
means = []
stds = []
labels_list = []
for sim in results:
    vals = results[sim]['IS_V_daily'][start_day:] * 100
    means.append(np.nanmean(vals))
    stds.append(np.nanstd(vals))
    labels_list.append(labels[sim])
# -------------------------------------------------------------------------------------
# plot V_q = Vol fraction where |q| < threshold
fig3, ax3 = plt.subplots(figsize=(18, 12))
for i, sim in enumerate(results):
    IS_V_mean = means[i]
    s = f"{IS_V_mean:.2e}"
    base, exp = s.split("e")
    text = f"${base} \\times 10^{{{int(exp)}}}$"

    ax3.plot(results[sim]['date'][start_day:], results[sim]['q_abs_lt_thresh_vol'][start_day:] * 100,
             label=rf"{labels[sim]}, Mean IS = {text} $s^{-3}$", 
             color=colors[i], linestyle=styles[i], linewidth=linewidths[i])

ax3.set_ylabel(fr"$|q| <  {formatted_thresh}$ Volume ($\%$)", fontsize=48)
ax3.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=30, ha='right')
ax3.tick_params(axis='x', labelsize=48)

ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax3.yaxis.get_offset_text().set_fontsize(48)
ax3.tick_params(axis='y', labelcolor='k', labelsize=48)


ax3.grid(True)
ax3.legend(fontsize=24, loc=1)

plt.tight_layout()
fig3.savefig(f"compare_q_lt_thresh_vol_z1_{z1}_z2_{z2}_nt_{indxRange}.png", dpi=300)
fig3.savefig(f"IS_V.pdf", dpi=300)

# -------------------------------------------------------------------------------------
