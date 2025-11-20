import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from skimage.measure import regionprops, label
from math import radians, sin, cos, sqrt, atan2
import matplotlib.cm as cm
import imageio.v3 as iio
import matplotlib as mpl
import pandas as pd
import trackpy as tp
# ----------------------------------------
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
##########################################
indxRange = '952_3877'
vlevel = -4000
iso_depth = -vlevel # for bathy contours
f0 = 0.909e-4 # Coriolis freq.
fig_dir = './figs'
# -------------------------------
trackfile_path= './trackfiles'
sim_paths = {
    'nesm_1km': {
        'grid_file': '../fig_3/NESM_grd.nc',
        'prefix': 'nesm_2019_2020_'

    },
    'nesm_5km': {
        'grid_file': '../fig_3/NESM_grd_5km.nc',
        'prefix': 'nesm_2019_2020_5km_'
    },
    'sm3_1km': {
        'grid_file': '../fig_3/SM3_grd.nc',
        'prefix': 'sm3_2019_2020_'
    }
}

# ----------------------------------------
for config, paths in sim_paths.items():
    prefix = paths['prefix']
    tracks_file = f'{trackfile_path}/{prefix}vortex_tracks_vlevel_{vlevel}_nt_{indxRange}.nc'
    
    ds_tracks = xr.open_dataset(tracks_file)
    
    # Extract arrays
    track_ids = ds_tracks['track_id'].values
    track_lengths = ds_tracks['track_length'].values
    latitudes = ds_tracks['latitude'].values
    longitudes = ds_tracks['longitude'].values
    n_tracks = len(track_ids)
    
    print(f"Loaded {n_tracks} tracks from {tracks_file}")
    # -----------------------------------
    # load grid
    grid_file = paths['grid_file'] 
    ds = xr.open_dataset(grid_file)
    lat_rho = ds['lat_rho'].transpose("eta_rho", "xi_rho").values
    lon_rho = ds['lon_rho'].transpose("eta_rho", "xi_rho").values
    h = ds['h'].transpose("eta_rho", "xi_rho").values
    mask_rho = ds['mask_rho'].transpose("eta_rho", "xi_rho").values
    ds.close()

    # Convert to lat/lon bounds using edges of the box
    lat_rho_vec = lat_rho[:, 0]
    lon_rho_vec = lon_rho[0, :]
    # -----------------------------------
    # plot all tracks on bathymetry
    fname=f"{fig_dir}/{config}_all_tracks.png"
    highlight_lifetime_days = 7
    plt.figure(figsize=(10, 8))
     
    cs = plt.contour(lon_rho, lat_rho, h, levels=[iso_depth], colors='k', linewidths=1.2)
    plt.clabel(cs, fmt={iso_depth: f'{iso_depth:.0f} m'}, fontsize=9, colors='k')

    for i in range(n_tracks):
        track_length = track_lengths[i]
        
        lats = latitudes[i, :track_length]
        lons = longitudes[i, :track_length]
        lifetime_days = track_length * 3 / 24
        
        # highlight if lifetime_days > 7 days
        lw = 4 if lifetime_days > highlight_lifetime_days else 1
        
        plt.plot(lons, lats, '-', alpha=0.8, linewidth=lw)
        plt.plot(lons[0], lats[0], 'o', markersize=3)
        plt.plot(lons[-1], lats[-1], 'x', markersize=3)
    
    plt.xlabel('Longitude', fontsize=24, labelpad=10)
    plt.ylabel('Latitude', fontsize=24, labelpad=10)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)    
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close()

    print(f"done saving: {fname}")
    # -----------------------------------
    ds_tracks.close()




# ----------------------------------------
