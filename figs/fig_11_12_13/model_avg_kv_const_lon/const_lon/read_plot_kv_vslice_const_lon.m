% plot kv_vslice_const_lon based on plot_kv_vslice_const_lon.m
%--------------------------------------------------------------------------
aghor_extras_path = '../../../../../aghor_extras/';
addpath(fullfile(aghor_extras_path, 'export_fig'));
addpath(fullfile(aghor_extras_path, 'm_map'));
addpath(fullfile(aghor_extras_path, 'cmap_manual'));
%--------------------------------------------------------------------------

%-----------------------------
% Plot Kv vertical slice for NESM 1km config
case_name = 'nesm_2019_2020';
% lonval = -63.25;
lonval = -63.15;

% indxRange = 952:3877; % entire year
% indxRange=2176:2296; % Oct 1 - Oct 15, 2019
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
indxRange=3024:3152; % Jan 15 - Jan 30, 2020

plots_path = './';

lat_min = 38; lat_max = 39.3;
z_min = -5500; z_max = 0;
cbar = false;           % suppress colorbar for panel-style plots
lat_idx = [];           % no vertical line

show_xlabels = true;    % show only on bottom panel
show_ylabels = true;    % show only on leftmost panel

plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
    lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
    show_xlabels, show_ylabels);
%-----------------------------
%-----------------------------
% Plot Kv vertical slice for NESM 5km config
case_name = 'nesm_2019_2020_5km';
% lonval = -63.26;
lonval = -63.14;

plots_path = './';

cbar = false;           % no colorbar here
lat_idx = [];           % no vertical line

show_xlabels = false;   % no x labels here
show_ylabels = false;   % no y labels here

plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
    lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
    show_xlabels, show_ylabels);
%-----------------------------
%-----------------------------
% Plot Kv vertical slice for NESM 5km_to_1km config
% case_name = 'nesm_2019_2020_5km_to_1km';
% % lonval = -63.26;
% lonval = -63.16;
% 
% plots_path = './';
% 
% cbar = false;           % no colorbar here
% lat_idx = [];           % no vertical line
% 
% show_xlabels = false;   % no x labels here
% show_ylabels = false;   % no y labels here
% 
% plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
%     lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
%     show_xlabels, show_ylabels);
%-----------------------------
%-----------------------------
% Plot Kv vertical slice for SM3 config
case_name = 'sm3_2019_2020';
% lonval = -63.25;
lonval = -63.15;

plots_path = './';

cbar = true;            % show colorbar only here
lat_idx = [];           % no vertical line

show_xlabels = false;   % no x labels here
show_ylabels = false;   % no y labels here

plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
    lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
    show_xlabels, show_ylabels);
%-----------------------------
%-----------------------------
% Plot Kv vertical slice for SM1 config
% case_name = 'sm1_2019_2020';
% % lonval = -63.25;
% lonval = -63.15;
% 
% plots_path = './';
% 
% cbar = true;            % show colorbar only here
% lat_idx = [];           % no vertical line
% 
% show_xlabels = false;   % no x labels here
% show_ylabels = false;   % no y labels here
% 
% plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
%     lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
%     show_xlabels, show_ylabels);
%-----------------------------

