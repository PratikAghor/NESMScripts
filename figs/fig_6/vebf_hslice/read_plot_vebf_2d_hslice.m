%-------------------------------------------------------
% use plot_vebf_hslice function to plot vebf_hslices
%-------------------------------------------------------
aghor_extras_path = '/media/aghor/aghor/prometheus_backup/Aghor/GT/GT_project/aghor_extras/';
addpath(genpath(fullfile(aghor_extras_path, 'export_fig')));
addpath(genpath(fullfile(aghor_extras_path, 'm_map')));
addpath(genpath(fullfile(aghor_extras_path, 'cmap_manual')));
addpath(genpath(fullfile(aghor_extras_path, 'cmap_manual/centered-colormap/')));
%-------------------------------------------------------
% Common params
plots_path = './';
vlevel = -4000;
lon_min = -64.5; lon_max = -62.8;
lat_min = 38; lat_max = 39.5;
cMax = 3e-7;
%-------------------------------------------------------
% NESM 1km
cbar = false; ref_arrow = false;
show_xlabel = true; show_ylabel = true;
show_xticks = true; show_yticks = true;
xyskip = 15;
lat_idx_arr = ([133]);
box = "NE_1km";

plot_vebf_hslice('nesm_2019_2020', vlevel, plots_path, ...
    lon_min, lon_max, lat_min, lat_max, lat_idx_arr, box, cbar, ref_arrow, cMax, xyskip, ...
    show_xlabel, show_ylabel, show_xticks, show_yticks);
%-------------------------------------------------------
% NESM 5km
cbar = false; ref_arrow = false;
show_xlabel = false; show_ylabel = false;
show_xticks = false; show_yticks = false;
xyskip = 3;
lat_idx_arr = ([27]);
box = "NE_5km";

plot_vebf_hslice('nesm_2019_2020_5km', vlevel, plots_path, ...
    lon_min, lon_max, lat_min, lat_max, lat_idx_arr, box, cbar, ref_arrow, cMax, xyskip, ...
    show_xlabel, show_ylabel, show_xticks, show_yticks);
%-------------------------------------------------------
%-------------------------------------------------------
% SM3 1km
cbar = true; ref_arrow = false;
show_xlabel = false; show_ylabel = false;
show_xticks = false; show_yticks = false;
xyskip = 15;
lat_idx_arr = ([133]);
box = "NE_1km";

plot_vebf_hslice('sm3_2019_2020', vlevel, plots_path, ...
    lon_min, lon_max, lat_min, lat_max, lat_idx_arr, box, cbar, ref_arrow, cMax, xyskip, ...
    show_xlabel, show_ylabel, show_xticks, show_yticks);
%-------------------------------------------------------
%-------------------------------------------------------
