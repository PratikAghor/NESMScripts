import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from netCDF4 import Dataset

vlevel = -4000
indxRange = '0_3877'
cases = [
    'nesm_2019_2020',
    'nesm_2019_2020_5km',
    'sm3_2019_2020'
]

# Mapping from case file prefix to plot label
labels = {
    'nesm_2019_2020': r'NESM (1km)',
    'nesm_2019_2020_5km': r'NESM (5km)',
    'sm3_2019_2020': r'SM3 (1km)',
}

colors = ['black', 'tab:cyan', 'tab:red']

# Convert frequencies to Hz (from radians/s)
f_coriolis = 0.909e-4 / (2 * np.pi)  # Hz
m2_freq = 1.40519e-4 / (2 * np.pi)   # Hz


plt.figure(figsize=(12, 6))

for i, (case, label) in enumerate(labels.items()):
    ncfile = f'{case}_box_avg_w_t_hslice_vlevel_{vlevel}_nt_{indxRange}.nc'
    with Dataset(ncfile, 'r') as nc:
        w = nc.variables['box_avg_w_t'][:]

    dt = 3 * 3600
    fs = 1 / dt
    frequencies, psd = welch(w, fs=fs, nperseg=512)

    plt.loglog(frequencies, psd, label=label, color=colors[i])

# Vertical lines for M2 and f
plt.axvline(m2_freq, color='r', linestyle='--', label='M2')
plt.axvline(f_coriolis, color='b', linestyle='--', label='f')

plt.xlabel(r'Frequency [Hz]', fontsize=24)
plt.ylabel(r'PSD $[(\mathrm{m}^2\,\mathrm{s}^{-2})/\mathrm{Hz}]$', fontsize=24)
# plt.title(f'Vertical Velocity PSD at z = {vlevel} m')
plt.xlim([1e-6, 4e-5])
plt.legend(fontsize=20, loc=3)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.grid(True, which='both', ls='--')
plt.tight_layout()
plt.savefig(f'w_psd_vlevel_{vlevel}.pdf')
# plt.show()
