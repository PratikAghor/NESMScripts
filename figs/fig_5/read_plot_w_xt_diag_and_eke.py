# --------------------
# Hovmoller for w_t_daily
# --------------------
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from datetime import datetime, timedelta

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})

# --------------------
# Inputs
# --------------------
indxRange = '952_3877'   # entire year
eke_indxRange = '0_3877'   # entire year

vlevel = -4000
base_date = datetime(2019, 5, 1)  # if indxRange == '952_3877'
eke_base_date = datetime(2019, 1, 1)  # if indxRange == '0_3877'

# start day
w_t0 = 0
eke_t0 = 120
               
n_files_per_day = 8  # data saved 3-hourly

# File locations and name prefixes
configs = {
    'nesm_1km': {
        'w': './',
        'eke_path':'./',
        'prefix': 'nesm_2019_2020_'
    },
    'nesm_5km': {
        'w': './',
        'eke_path':'./',
        'prefix': 'nesm_2019_2020_5km_'
    },
    'sm3_1km': {
        'w': './',
        'eke_path':'./',
        'prefix': 'sm3_2019_2020_'
    },
}

labels = {
    'nesm_1km': r'NESM (1km)',
    'nesm_5km': r'NESM (5km)',
    'sm3_1km': r'SM3 (1km)',
}

latvals = {
    'nesm_1km': 38.70,
    'nesm_5km': 38.68,
    'sm3_1km': 38.70,
}

# Plot settings
w_cmap = 'RdBu_r'
w_vmin, w_vmax = -1e-2, 1e-2

# --------------------
# daily_avg, 8 files a day
# w data is 2d [Nt x Nx], eke data is 1d [Nt]
# --------------------
def daily_avg(data, n):
    Nt, Nx = data.shape
    r = Nt % n
    if r != 0:
        data = np.vstack([data, np.repeat(data[-1:], n - r, axis=0)])
        Nt = data.shape[0]
    Nd = Nt // n
    return np.array([data[i*n:(i+1)*n].mean(axis=0) for i in range(Nd)])
# --------------------
def daily_avg_eke(timeseries, n):
    Nt = len(timeseries)
    Nt_avg = Nt // n_files_per_day

    print(f"Nt_avg = {Nt_avg}")

    daily_avg_series = np.zeros(Nt_avg)
    for ll in range(0, Nt_avg):
        daily_avg_series[ll] = np.nanmean(
            timeseries[ll * n_files_per_day : (ll + 1) * n_files_per_day]
        )
    
    return daily_avg_series

# --------------------
# load w daily data 
# --------------------
def load_w_daily(key):
    base_dir = configs[key]['w']
    prefix = configs[key]['prefix']
    latval = latvals[key]
    fname = f"{prefix}w_temp_salt_slice_const_lat_{latval:.2f}_const_z_{vlevel}_nt_{indxRange}.nc"
    path = os.path.join(base_dir, fname)
    
    ds = xr.open_dataset(path)
    w_t = ds['w_t_vslice_const_lat_const_z'].transpose('Nt', 'Nx').values
    lon = ds['lon_rho_vec'].transpose('Nx').values.astype(float)

    w_d = daily_avg(w_t, n_files_per_day)
    Nd = w_d.shape[0]
    dates = [base_date + timedelta(days=i) for i in range(Nd)]
    return w_d, lon, dates


def load_eke_daily(key):
    base_dir = configs[key]['eke_path']
    prefix = configs[key]['prefix']
    fname = f"{prefix}box_avg_eke_t_hslice_vlevel_{vlevel}_nt_{eke_indxRange}.nc"
    path = os.path.join(base_dir, fname)
    
    ds = xr.open_dataset(path)
    eke_t = ds['box_avg_eke_t'].transpose('Nt').values
    eke_t = np.append(eke_t, eke_t[-2:])

    eke_d = daily_avg_eke(eke_t, n_files_per_day)
    Nd = eke_d.shape[0]
    dates = [eke_base_date + timedelta(days=i) for i in range(Nd)]
    return eke_d, dates

# --------------------
# Load all
# --------------------
keys = ['nesm_1km', 'nesm_5km', 'sm3_1km']

W = {}
EKE = {}
for k in keys:
    w_d, lon, dates = load_w_daily(k)
    eke_d, eke_dates = load_eke_daily(k)

    # remove transients if any
    w_d = w_d[w_t0:-1]
    eke_d = eke_d[eke_t0:]
    dates = eke_dates[eke_t0:]

    W[k] = dict(w_d=w_d, lon=lon, dates=dates)
    EKE[k] = dict(eke_d=eke_d, dates=dates)
    print(f"{labels[k]}: w_d.shape={w_d.shape}, lon.size={lon.size}")
    print(f"{labels[k]}: eke_d.shape={eke_d.shape}")

# --------------------
# Hovmoller
# --------------------
latval = latvals['nesm_1km'] # for naming purposes

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(20, 12), sharex=True)  # shorter height than 18

for ax in axes:
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0 + 0.02, pos.width, pos.height * 0.75])  # reduce height

for ax, k in zip(axes, keys):
    dates = W[k]['dates']
    lon = W[k]['lon']
    Z = W[k]['w_d'].T
    X, Y = np.meshgrid(mdates.date2num(dates), lon)
    extent = [mdates.date2num(dates[0]),
          mdates.date2num(dates[-1]),
          lon[0],
          lon[-1]]

    pc = ax.imshow(Z, aspect='auto', origin='lower',
                   extent=extent,
                   cmap=w_cmap, vmin=w_vmin, vmax=w_vmax,
                   interpolation='spline36')

    ax.set_title(rf"{labels[k]}")
    ax.set_ylabel('Longitude', fontsize=32)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right')
    ax.tick_params(axis='x', labelsize=32)
    ax.set_ylim([-63.5, -62.5])
    yticks = [-63.5, -63.0, -62.5]
    ax.set_yticks(yticks)

    # superpose EKE
    ax2 = ax.twinx()
    eke_dates = EKE[k]['dates']
    eke_d = EKE[k]['eke_d']
    ax2.plot(eke_dates, eke_d, color='black', linewidth=4, alpha=0.7)

    ax2.set_ylabel(r'$\langle EKE \rangle$[m$^2$/s$^2$]', color='k', rotation=270, 
        labelpad=50,
        fontsize=28)
    ax2.set_ylim([0, 0.1])
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax2.yaxis.get_offset_text().set_fontsize(24)
    ax2.tick_params(axis='y', labelcolor='k', labelsize=24)

    ax2.grid(False)  # Turn off grid for secondary axis

cbar_ax = fig.add_axes([0.15, 0.97, 0.7, 0.01])
cbar = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label(r'$w [m/s]$')
cbar.set_ticks([w_vmin, w_vmax])

fig.subplots_adjust(top=0.88, hspace=0.25)
fig.savefig(f'compare_w_t_const_lat_{latval:.2f}_daily_hovmoller_with_eke.png', dpi=300, bbox_inches='tight')
fig.savefig(f'compare_w_t_const_lat_{latval:.2f}_daily_hovmoller_with_eke.pdf')
plt.close(fig)
# --------------------
