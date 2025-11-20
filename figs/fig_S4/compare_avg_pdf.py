import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from scipy.stats import gaussian_kde

#-------------------------------------------
# Paths and settings
plots_path = './'
vlevel_top = 1850
n_files_per_day = 8
start_day_for_averaging = 120
indxRange='3024_3152' # Jan 15 - Jan 30, 2025

# paths
sim_paths = {
    'nesm_1km': {
        'vort': './vort_pdf_data/',
        'prefix': 'nesm_2019_2020_'
    },
    'nesm_5km': {
        'vort': './vort_pdf_data/',
        'prefix': 'nesm_2019_2020_5km_'
    },
    'sm3_1km': {
        'vort': './vort_pdf_data/',
        'prefix': 'sm3_2019_2020_'
    }
}

labels = {
    'nesm_1km': r'NESM 1km',
    'nesm_5km': r'NESM 5km',
    'sm3_1km': r'SM3 1km'
}

colors = {
    'nesm_1km': 'k',
    'nesm_5km': 'tab:cyan',
    'sm3_1km': 'tab:red'
}

#-------------------------------------------
# plot pdfs
fig, ax = plt.subplots(figsize=(10, 6))

for sim in sim_paths.keys():
    nc_file = sim_paths[sim]['vort'] + sim_paths[sim]['prefix'] + 'dimless_vort_kde_avg_1850m_to_bottom.nc'

    ds = xr.open_dataset(nc_file)
    vort_bins = ds['vort_bins'].values
    mean_pdf = ds['mean_pdf'].values
    ax.plot(vort_bins, mean_pdf, label=labels[sim], linewidth=2, color=colors[sim])

    # Normalize
    pdf_norm = mean_pdf / np.trapz(mean_pdf, vort_bins)

    # Compute mean and std. dev.
    mu = np.trapz(pdf_norm * vort_bins, vort_bins)
    variance = np.trapz(pdf_norm * (vort_bins - mu) ** 2, vort_bins)
    sigma = np.sqrt(variance)

    print(f"{labels[sim]}: mean = {mu:.4f}, std. dev. = {sigma:.4f}")

# suit-putting
ax.set_yscale('log')
ax.set_ylim(1e-4, 2e1)
ax.set_xlim(-1, 1)

ax.set_xlabel(r'$\zeta/f$', fontsize=20)
ax.set_ylabel('Time-averaged PDF (log scale)', fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.legend(fontsize=16)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)

fig.tight_layout()
fig.savefig(plots_path + 'dimless_vort_kde_comparison_nt_' + indxRange + '.pdf', dpi=300)
#-------------------------------------------